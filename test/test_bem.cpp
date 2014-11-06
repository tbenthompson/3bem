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
#include "taylor.h"

struct IntegrationProb {
    IntegrationProb():
        src_locs({{{0,0,0},{2,0,0},{0,1,0}}}),
        src_vals({1.0, 1.0, 1.0}),
        obs_loc({2.0, 2.0, 2.0}),
        obs_n({1.0, 0.0, 0.0}),
        q(tri_gauss(2)),
        face(Mesh{src_locs, {{0,1,2}}}, 0)
    { }
    std::vector<Vec3<double>> src_locs;
    Vec3<double> src_vals;
    Vec3<double> obs_loc;
    Vec3<double> obs_n;
    QuadratureRule2D q;
    FaceInfo face;
};

TEST_FIXTURE(IntegrationProb, IntegralOne) {
    double abc = integrate(q, [](double x_hat, double y_hat){return 1.0;});
    CHECK_CLOSE(abc, 0.5, 1e-12);
    auto kernel = one<double>;
    double result = integral<double>(q, kernel, face, src_vals, obs_loc, obs_n);
    CHECK_CLOSE(result, 1.0, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceSingle) {
    auto kernel = laplace_single<double>;
    double result = integral<double>(q, kernel, face, src_vals, obs_loc, obs_n);
    CHECK_CLOSE(result, 0.0269063, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceDouble) {
    auto kernel = laplace_double<double>;
    double result = integral<double>(q, kernel, face, src_vals, obs_loc, obs_n);
    CHECK_CLOSE(result, 0.00621003, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralBasis) {
    auto kernel = laplace_double<double>;
    src_vals = {1.0, 0.7, 0.3};
    double result = integral<double>(q, kernel, face, src_vals, obs_loc, obs_n);
    Vec3<double> basis =
        basis_integrals<double>(q, kernel, face, obs_loc, obs_n);
    double result2 = dot(basis, src_vals);
    CHECK_CLOSE(result, result2, 1e-5);
}

TEST_FIXTURE(IntegrationProb, TaylorIntegral) {
    auto Tk = laplace_double<Td<taylor_degree>>;

    // The goal point
    obs_loc = {0.5, 0.2, 0.0};

    //Define an expansion point for the taylor expansion
    auto x = Td<taylor_degree>::var(-1.0);
    double s1 = -0.1;
    std::array<Td<taylor_degree>,3> t_exp_pt = {
        obs_loc[0], obs_loc[1], s1 * x + obs_loc[2]
    };

    // Taylor expand
    int high_order = 50;
    auto acc1 = integral<Td<taylor_degree>>(tri_gauss(high_order), Tk, face,
                                               src_vals, t_exp_pt, obs_n);
    for (int i = 20; i < 21; i+=2) {
        auto taylor_exp = integral<Td<taylor_degree>>(tri_gauss(i), Tk, face,
                                            src_vals, t_exp_pt, obs_n);
        double exp1 = taylor_exp.eval(1);
        std::cout << "Far point with Gauss order: " << i << " " 
                  << taylor_exp - acc1 << " " << exp1 << std::endl;

        std::array<double,taylor_degree+1> a = taylor_exp.c;
        std::array<double,taylor_degree+1> s;
        s[0] = a[0];
        s[1] = s[0] + a[1];
        std::cout << s[0] << std::endl;
        std::cout << s[1] << std::endl;
        for (int i = 2; i < taylor_degree + 1; i++) {
            s[i] = s[i - 1] + a[i];
            std::cout << s[i] << std::endl;
        }
        std::cout << taylor_exp << std::endl;
    }
}

TEST(Richardson) {
    std::vector<double> input = {0.5, 0.3, 0.2, 0.15};
    double result = richardson_step(input);
    CHECK_CLOSE(result, 0.1, 1e-12);
}

TEST(RichardsonZeros) {
    std::vector<double> input = {0.0, 0.0, 0.0, 0.0};
    double result = richardson_step(input);
    CHECK_CLOSE(result, 0.0, 1e-12);
}

TEST_FIXTURE(IntegrationProb, RichardsonIntegral) {
    auto q = tri_gauss(3);
    auto kernel = laplace_single<double>;
    double offset = 0.5;
    std::vector<double> vals;
    for (int i = 0; i < 4; i++) {
        obs_loc = {2.0, 2.0, 2.0 + offset};
        vals.push_back(integral<double>(q, kernel, face, src_vals, obs_loc, {1, 0, 0}));
        offset /= 2;
    }
    double result = richardson_step(vals);
    CHECK_CLOSE(result, 0.0269063, 1e-6);
}

struct EvalProb {
    EvalProb(int refine_level, int near_eval, int gauss_order, Kernel k,
             TaylorKernel Tk, Vec3<double> center = Vec3<double>{0,0,0},
             double r = 3.0):
        q(tri_gauss(gauss_order)),
        kernel(k),
        taylor_kernel(Tk),
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
        return eval_integral_equation(sphere, q, kernel, taylor_kernel, ne, obs_pt,
                                      obs_normal, obs_len_scale, src_strength);
    }
    double go_row() {
        auto row = integral_equation_vector(sphere, q, kernel, taylor_kernel, ne, obs_pt,
                                      obs_normal, obs_len_scale);
        double row_sum = 0.0;
        for(std::size_t i = 0; i < row.size(); i++) {
            row_sum += row[i] * src_strength[i];
        }
        return row_sum;
    }

    Mesh sphere;
    QuadratureRule2D q;
    Kernel kernel;
    TaylorKernel taylor_kernel;
    NearEval ne;
    Vec3<double> obs_pt;
    Vec3<double> obs_normal;
    double obs_len_scale;
    std::vector<double> src_strength;
}; 

TEST(EvalIntegralEquationSphereSurfaceArea) {
    EvalProb ep(5, 3, 2, one<double>, one<Td<taylor_degree>>);
    double result = ep.go();
    double result2 = ep.go_row();
    double exact_surf_area = 4*M_PI*9;
    CHECK_CLOSE(result, exact_surf_area, 1e-1);
    CHECK_CLOSE(result2, exact_surf_area, 1e-1);
}

TEST(ConstantLaplace) {
    EvalProb ep(5, 3, 2, laplace_double<double>,
                laplace_double<Td<taylor_degree>>);
    double result = ep.go();
    double result2 = ep.go_row();
    CHECK_CLOSE(result, 1.0, 1e-3);
    CHECK_CLOSE(result2, 1.0, 1e-3);
}

TEST(ConstantLaplaceBoundary) {
    EvalProb ep(2, 3, 3, laplace_double<double>,
                laplace_double<Td<taylor_degree>>);
    int n_verts = ep.sphere.vertices.size();
    // n_verts = 1;
    for (int i = 0; i < n_verts;i++) {
        ep.obs_pt = ep.sphere.vertices[i];
        ep.obs_normal = -normalized(ep.obs_pt);
        double result = ep.go();
        double result2 = ep.go();
        CHECK_CLOSE(result, 1.0, 1e-2);
        CHECK_CLOSE(result2, 1.0, 1e-2);
    }
}

TEST(MatrixRowVsEval) {
    EvalProb ep(5, 3, 2, laplace_double<double>,
                laplace_double<Td<taylor_degree>>);
    double result = ep.go();
    double result2 = ep.go_row();
    CHECK_CLOSE(result, 1.0, 1e-3);
    CHECK_CLOSE(result2, 1.0, 1e-3);
}

TEST(Whoa_baNANnas) {
    EvalProb ep(2, 3, 4, laplace_double<double>,
                laplace_double<Td<taylor_degree>>, {5,0,0}, 3.0);
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
    for (int i = 0; i < 3; i++) {
        sphere = refine_mesh(sphere);
    }
    sphere = clean_mesh(sphere);
    auto q = tri_gauss(2);
    int n_verts = sphere.vertices.size();
    std::vector<double> str(n_verts, 1.0);
    auto res = direct_interact(sphere, sphere, q, q,
                               one<double>, one<Td<taylor_degree>>, str, 2);
    auto matrix = interact_matrix(sphere, sphere, q, q,
                                  one<double>, one<Td<taylor_degree>>, 2);
    std::vector<double> res2(n_verts, 0.0);
    for (int i = 0; i < n_verts; i++) {
        for (int j = 0; j < n_verts; j++) {
            res2[i] += matrix[i][j] * str[j]; 
        }
    }
    double total = 0.0;
    for (auto r: res) {
        total += r;
    }
    double total2 = 0.0;
    for (auto r2: res2) {
        total2 += r2;
    }
    double sa2 = pow(4 * M_PI, 2);
    double error = std::fabs((sa2 - total) / sa2);
    double error2 = std::fabs((sa2 - total2) / sa2);
    CHECK_CLOSE(error, 0.0, 6.4e-2);
    CHECK_CLOSE(error2, 0.0, 6.4e-2);
}

//TODO: Fixture for this and the next one.
TEST(DirectInteractConstantLaplace) {
    Mesh sphere = clean_mesh(sphere_mesh({0,0,0}, 1.0));
    for (int i = 0; i < 3; i++) {
        sphere = refine_mesh(sphere);
    }
    sphere = clean_mesh(sphere);
    auto q = tri_gauss(2);
    std::vector<double> str(sphere.vertices.size(), 1.0);
    auto res0 = direct_interact(sphere, sphere, q, q,
                                laplace_double<double>, 
                                laplace_double<Td<taylor_degree>>, str, 2);
    auto res1 = direct_interact(sphere, sphere, q, q,
                                laplace_single<double>,
                                laplace_double<Td<taylor_degree>>, str, 2);
    auto res2 = mass_term(sphere, q, str);
    CHECK_ARRAY_CLOSE(res0, res2, sphere.vertices.size(), 1e-2);
    CHECK_ARRAY_CLOSE(res1, res2, sphere.vertices.size(), 1e-2);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
