#include "UnitTest++.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

#include "kernels.h"
#include "bem.h"
#include "numerics.h"
#include "quadrature.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "util.h"
#include "shared.h"

struct IntegrationProb {
    IntegrationProb() {
        src_locs = {{{0,0,0},{2,0,0},{0,1,0}}};
        src_vals = {1.0, 1.0, 1.0};
        obs_loc = {2.0, 2.0, 2.0};
    }
    std::array<Vec3<double>,3> src_locs;
    Vec3<double> src_vals;
    Vec3<double> obs_loc;
};

TEST_FIXTURE(IntegrationProb, IntegralOne) {
    auto q = tri_gauss(2);
    double abc = integrate(q, [](double x_hat, double y_hat){return 1.0;});
    CHECK_CLOSE(abc, 0.5, 1e-12);
    KernelFnc kernel = one;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc, {0,0,1});
    CHECK_CLOSE(result, 1.0, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceSingle) {
    auto q = tri_gauss(2);
    KernelFnc kernel = laplace_single;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc, {1, 0, 0});
    CHECK_CLOSE(result, 0.0269063, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceDouble) {
    auto q = tri_gauss(2);
    KernelFnc kernel = laplace_double;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc, {1,0,0});
    CHECK_CLOSE(result, 0.00621003, 1e-5);
}

//TODO: Test a kernel that depends on the observation normal vector -- adjoint 
//double layer

TEST(Richardson) {
    double result = richardson_step({0.5, 0.3, 0.2, 0.15});
    CHECK_CLOSE(result, 0.1, 1e-12);
}

TEST(RichardsonZeros) {
    double result = richardson_step({0.0, 0.0, 0.0, 0.0});
    CHECK_CLOSE(result, 0.0, 1e-12);
}

TEST_FIXTURE(IntegrationProb, RichardsonIntegral) {
    auto q = tri_gauss(3);
    KernelFnc kernel = laplace_single;
    double offset = 0.5;
    std::vector<double> vals;
    for (int i = 0; i < 4; i++) {
        obs_loc = {2.0, 2.0, 2.0 + offset};
        vals.push_back(integral(q, kernel, src_locs, src_vals, obs_loc, {1, 0, 0}));
        offset /= 2;
    }
    double result = richardson_step(vals);
    CHECK_CLOSE(result, 0.0269063, 1e-6);
}

struct EvalProb {
    EvalProb(int refine_level, int near_eval, int gauss_order, KernelFnc k,
             Vec3<double> center = Vec3<double>{0,0,0}, double r = 3.0):
        q(tri_gauss(gauss_order)),
        kernel(k),
        ne(near_eval)
    {
        sphere = sphere_mesh(center, r);
        for (int i = 0; i < refine_level; i++) {
            sphere = refine_mesh(sphere);
        }
        obs_pt = random_pt();
        obs_len_scale = get_len_scale(sphere, 0, gauss_order);
        obs_normal = random_pt();
        src_strength = std::vector<double>(sphere.vertices.size(), 1.0);
    }

    double go() {
        return eval_integral_equation(sphere, q, kernel, ne, obs_pt,
                                      obs_normal, obs_len_scale, src_strength);
    }
    Mesh sphere;
    QuadratureRule2D q;
    KernelFnc kernel;
    NearEval ne;
    Vec3<double> obs_pt;
    Vec3<double> obs_normal;
    double obs_len_scale;
    std::vector<double> src_strength;
}; 

TEST(EvalIntegralEquationSphereSurfaceArea) {
    EvalProb ep(9, 3, 2, one);
    double result = ep.go();
    double exact_surf_area = 4*M_PI*9;
    CHECK_CLOSE(result, exact_surf_area, 1e-3);
}

TEST(ConstantLaplace) {
    EvalProb ep(5, 3, 2, laplace_double);
    double result = ep.go();
    CHECK_CLOSE(result, 1.0, 1e-3);
}

TEST(ConstantLaplaceBoundary) {
    EvalProb ep(2, 3, 3, laplace_double);
    int n_verts = ep.sphere.vertices.size();
    // n_verts = 1;
    for (int i = 0; i < n_verts;i++) {
        ep.obs_pt = ep.sphere.vertices[i];
        ep.obs_normal = -normalized(ep.obs_pt);
        double result = ep.go();
        CHECK_CLOSE(result, 1.0, 1e-2);
    }
}

TEST(Whoa_baNANnas) {
    EvalProb ep(2, 3, 4, laplace_double, {5,0,0}, 3.0);
    ep.obs_pt = {6.5309310892394858, 1.8371173070873836, -1.5309310892394863};
    ep.obs_normal = {-0.5773502691896254, -0.57735026918962595, 0.57735026918962595};
    ep.obs_len_scale = 0.0094194506485734651;
    
    double result = ep.go();
    CHECK_CLOSE(result, 0.561132, 1e-3);
}

TEST(MassTerm) {
    Mesh sphere = clean_mesh(sphere_mesh({0,0,0}, 1.0));
    std::vector<double> str(sphere.vertices.size(), 1.0);
    str[0] = 0.0;
    auto res = mass_term(sphere, tri_gauss(2), str);
    CHECK_EQUAL(res.size(), sphere.vertices.size());
    double true_area = 0.0;
    for (auto f: sphere.faces) {
        true_area += tri_area({
                sphere.vertices[f[0]],
                sphere.vertices[f[1]],
                sphere.vertices[f[2]]
            });
    }
    double mass_area = 0.0;
    for (auto r: res) {
        mass_area += r;
    }   
    CHECK_CLOSE(mass_area, (5.0 / 6.0) * true_area, 1e-12);
}

TEST(DirectInteractOne) {
    Mesh sphere = clean_mesh(sphere_mesh({0,0,0}, 1.0));
    for (int i = 0; i < 4; i++) {
        sphere = refine_mesh(sphere);
    }
    sphere = clean_mesh(sphere);
    auto q = tri_gauss(2);
    std::vector<double> str(sphere.vertices.size(), 1.0);
    auto res = direct_interact(sphere, sphere, q, q, one, str, 2);
    double total = 0.0;
    for (auto r: res) {
        total += r;
    }
    double sa2 = pow(4 * M_PI, 2);
    double error = std::fabs((sa2 - total) / sa2);
    CHECK_CLOSE(error, 0.0, 6.4e-3);
}

TEST(DirectInteractConstantLaplace) {
    Mesh sphere = clean_mesh(sphere_mesh({0,0,0}, 1.0));
    for (int i = 0; i < 3; i++) {
        sphere = refine_mesh(sphere);
    }
    sphere = clean_mesh(sphere);
    auto q = tri_gauss(2);
    std::vector<double> str(sphere.vertices.size(), 1.0);
    auto res0 = direct_interact(sphere, sphere, q, q, laplace_double, str, 2);
    auto res1 = direct_interact(sphere, sphere, q, q, laplace_single, str, 2);
    auto res2 = mass_term(sphere, q, str);
    CHECK_ARRAY_CLOSE(res0, res2, sphere.vertices.size(), 1e-2);
    CHECK_ARRAY_CLOSE(res1, res2, sphere.vertices.size(), 1e-2);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
