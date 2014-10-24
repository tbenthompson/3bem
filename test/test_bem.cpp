#include "UnitTest++.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

#include "bem.h"
#include "numerics.h"
#include "mesh.h"
#include "util.h"
#include "shared.h"

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
    KernelFnc kernel = BEMone;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc);
    CHECK_CLOSE(result, 1.0, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceSingle) {
    auto q = tri_gauss(2);
    KernelFnc kernel = BEMlaplace_single;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc);
    CHECK_CLOSE(result, 0.0269063, 1e-5);
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceDouble) {
    auto q = tri_gauss(2);
    KernelFnc kernel = BEMlaplace_double;
    double result = integral(q, kernel, src_locs, src_vals, obs_loc);
    CHECK_CLOSE(result, 0.00621003, 1e-5);
}

TEST(Richardson) {
    double result = richardson_step({0.5, 0.3, 0.2, 0.15});
    CHECK_CLOSE(result, 0.1, 1e-12);
}

TEST_FIXTURE(IntegrationProb, RichardsonIntegral) {
    auto q = tri_gauss(3);
    KernelFnc kernel = BEMlaplace_single;
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
    KernelFnc kernel = BEMone;
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
    KernelFnc kernel = BEMlaplace_double;
    NearEval ne(6);
    auto obs_pt = random_pt();
    auto obs_normal = random_pt();
    std::vector<double> src_strength(sphere.vertices.size(), 1.0);
    double result = eval_integral_equation(sphere, q, kernel, ne, obs_pt,
                                           obs_normal, src_strength);
    CHECK_CLOSE(result, 1.0, 1e-3);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
