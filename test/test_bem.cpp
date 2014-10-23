#include "UnitTest++.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

#include "bem.h"
#include "numerics.h"
#include "mesh.h"
#include "util.h"
#include <iostream>

inline double one_kernel(double r2, std::array<double,3> delta,
                         std::array<double,3> n) {
    return 1.0;
}

inline double laplace_single(double r2,
                             std::array<double,3> delta,
                             std::array<double,3> nsrc) {
    return 1.0 / (4.0 * M_PI * std::sqrt(r2));
}

inline double laplace_double(double r2,
                             std::array<double,3> delta,
                             std::array<double,3> nsrc) {
    return (nsrc[0] * delta[0] + nsrc[1] * delta[1] + nsrc[2] * delta[2]) / 
           (4 * M_PI * pow(r2, 1.5));
}

struct IntegrationProb {
    IntegrationProb() {
        src_locs = {{{0,0,0},{2,0,0},{0,1,0}}};
        src_vals = {1.0, 1.0, 1.0};
        obs_loc = {2.0, 2.0, 2.0};
    }
    std::array<std::array<double,3>,3> src_locs;
    std::array<double,3> src_vals;
    std::array<double,3> obs_loc;
};

TEST_FIXTURE(IntegrationProb, IntegralOne) {
    auto q = tri_gauss(2);
    double abc = integrate(q, [](double x_hat, double y_hat){return 1.0;});
    CHECK_CLOSE(abc, 0.5, 1e-12);
    KernelFnc kernel = one_kernel;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc);
    CHECK_CLOSE(result, 1.0, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceSingle) {
    auto q = tri_gauss(2);
    KernelFnc kernel = laplace_single;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc);
    CHECK_CLOSE(result, 0.0269063, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceDouble) {
    auto q = tri_gauss(2);
    KernelFnc kernel = laplace_double;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc);
    CHECK_CLOSE(result, 0.00621003, 1e-5);
}

TEST(Richardson) {
    double result = richardson_step({0.5, 0.3, 0.2, 0.15});
    CHECK_CLOSE(result, 0.1, 1e-12);
}

TEST_FIXTURE(IntegrationProb, RichardsonIntegral) {
    auto q = tri_gauss(3);
    KernelFnc kernel = laplace_single;
    double offset = 0.5;
    std::vector<double> vals;
    for (int i = 0; i < 4; i++) {
        obs_loc = {2.0, 2.0, 2.0 + offset};
        vals.push_back(integral(q, kernel, src_locs, src_vals, obs_loc));
        offset /= 2;
    }
    double result = richardson_step(vals);
    CHECK_CLOSE(result, 0.0269063, 1e-6);
}

TEST(EvalIntegralEquationSphereSurfaceArea) {
    //TODO: Convert to autocheck
    auto sphere = sphere_mesh({0, 0, 0}, 3.0);
    int refine_level = 9;
    for (int i = 0; i < refine_level; i++) {
        sphere = refine_mesh(sphere);
    }
    auto q = tri_gauss(3);
    KernelFnc kernel = one_kernel;
    NearEval ne(3);
    auto obs_pt = random_pt();
    auto obs_normal = random_pt();
    std::vector<double> src_strength(sphere.vertices.size(), 1.0);
    double result = eval_integral_equation(sphere, q, kernel, ne, obs_pt,
                                           obs_normal, src_strength);
    double exact_surf_area = 4*M_PI*9;
    CHECK_CLOSE(result, exact_surf_area, 1e-3);
}

TEST(ConstantLaplace) {
    //TODO: Convert to autocheck
    auto sphere = sphere_mesh({0, 0, 0}, 4.0);
    int refine_level = 4;
    for (int i = 0; i < refine_level; i++) {
        sphere = refine_mesh(sphere);
    }
    auto q = tri_gauss(2);
    KernelFnc kernel = laplace_double;
    NearEval ne(6);
    auto obs_pt = random_pt();
    auto obs_normal = random_pt();
    std::vector<double> src_strength(sphere.vertices.size(), 1.0);
    double result = eval_integral_equation(sphere, q, kernel, ne, obs_pt,
                                           obs_normal, src_strength);
    CHECK_CLOSE(result, 1.0, 1e-3);
}


double harmonic_u(std::array<double,3> x) {
    return 1.0 / hypot(x);
}

double harmonic_dudn(std::array<double,3> x,
                     std::array<double,3> center) {
    double nx = x[0] - center[0];
    double ny = x[1] - center[1];
    double nz = x[2] - center[2];
    double n_mag = hypot(nx, ny, nz);
    nx /= n_mag; ny /= n_mag; nz /= n_mag;
    double r = hypot(x);
    return -(x[0] * nx + x[1] * ny + x[2] * nz) / (r * r * r);
}


//TODO: This could be the start of a Vec3 class. 
class Pt {
public:
    std::array<double,3> x;
    friend std::ostream& operator<<(std::ostream& os, const Pt& obj)
    {
        os << "(" << obj.x[0] << ", " << obj.x[1] << ", " << obj.x[2] << ")";
        return os;
    }
};

class SpherePtGenerator {
public:
    SpherePtGenerator(std::array<double,3> center, double r):
        center(center),
        r(r)
    {}

    typedef Pt result_type;

    result_type operator()(size_t size = 0) {
        return Pt{random_pt_sphere(center, r)};
    }
    const std::array<double,3> center;
    const double r;
};

TEST(LaplaceHarmonic) {
    //THIS IS HOT!
    const std::array<double,3> center = {5, 0, 0};
    double r = 3.0;
    int refine_level = 5;
    int near_field = 2;
    int far_gauss_pts = 2;

    auto sphere = sphere_mesh(center, r);
    for (int i = 0; i < refine_level; i++) {
        sphere = refine_mesh(sphere);
    }
    auto q = tri_gauss(far_gauss_pts);
    KernelFnc K = laplace_single;
    KernelFnc Kdn = laplace_double;
    NearEval ne(near_field);

    int n_verts = sphere.vertices.size();
    std::vector<double> u(n_verts);
    for (int i = 0; i < n_verts; i++) {
        u[i] = harmonic_u(sphere.vertices[i]);
    }
    std::vector<double> dudn(n_verts);
    for (int i = 0; i < n_verts; i++) {
        dudn[i] = harmonic_dudn(sphere.vertices[i], center);
    }

    auto arb = ac::make_arbitrary(SpherePtGenerator(center, r));
    ac::check<Pt>(
        [&](Pt p) {
            auto obs_pt = p.x;
            auto obs_normal = normalize(diff(center, obs_pt));

            // const double cdist = std::sqrt(dist2<3>(obs_pt, center));
            /* std::cout << cdist << " " << obs_pt[0] << " " << obs_pt[1] << " " << obs_pt[2] <<  " " << obs_normal[0] << " " << obs_normal[1] << " " << obs_normal[2] << std::endl; */

            double result = eval_integral_equation(sphere, q, Kdn, ne, obs_pt,
                                                   obs_normal, u);
            result += eval_integral_equation(sphere, q, K, ne, obs_pt,
                                             obs_normal, dudn);
            double exact = 1.0 / hypot(obs_pt);
            // std::cout << result << " " << exact << std::endl;
            return std::fabs(exact - result) < 1e-2;
        }, 100, arb);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
