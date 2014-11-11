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
    IntegrationProb():
        src_locs({{{0,0,0},{2,0,0},{0,1,0}}}),
        src_vals({1.0, 1.0, 1.0}),
        obs_loc({2.0, 2.0, 2.0}),
        obs_n({1.0, 0.0, 0.0}),
        q(tri_gauss(2)),
        face(Mesh{src_locs, {{0,1,2}}}, 0)
    { }
    
    void go() {
        auto basis = integrate<Vec3<double>,2>(q, [&] (std::array<double,2> x_hat) {
                return eval_quad_pt(x_hat, kernel, face, obs_loc, obs_n);
            });
        result = dot(basis, src_vals);
    }

    void check() {
        CHECK_CLOSE(result, exact, 1e-5);
    }

    Kernel kernel;
    double result;
    double exact;
    std::vector<Vec3<double>> src_locs;
    Vec3<double> src_vals;
    Vec3<double> obs_loc;
    Vec3<double> obs_n;
    QuadRule2d q;
    FaceInfo face;
};

TEST_FIXTURE(IntegrationProb, IntegralOne) {
    double abc = integrate<double,2>(q, [](std::array<double,2> x_hat){return 1.0;});
    CHECK_CLOSE(abc, 0.5, 1e-12);
    kernel = one;
    exact = 1.0; go(); check();
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceSingle) {
    kernel = laplace_single;
    exact = -0.0269063; go(); check();
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceDouble) {
    kernel = laplace_double;
    exact = 0.00621003; go(); check();
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
    q = tri_gauss(3);
    kernel = laplace_single;
    double offset = 0.5;
    std::vector<double> vals;
    obs_n = {1,0,0};
    for (int i = 0; i < 5; i++) {
        obs_loc = {2.0, 2.0, 2.0 + offset};
        go();
        vals.push_back(result);
        offset /= 2;
    }
    double result = richardson_step(vals);
    CHECK_CLOSE(result, -0.0269063, 1e-6);
}

//TODO: This should be refactored a bit!
struct EvalProb {
    EvalProb(int refine_level, int near_eval, int gauss_order, const Kernel& k,
             Vec3<double> center = Vec3<double>{0,0,0},
             double r = 3.0):
        sphere(refine_clean(sphere_mesh(center,r), refine_level)),
        qs(gauss_order, gauss_order, gauss_order, near_eval, 2.0, 1e-2),
        kernel(k),
        obs_pt(random_pt()),
        obs_n(random_pt()),
        obs_length_scale(get_len_scale(sphere, 0, gauss_order)),
        src_strength(std::vector<double>(sphere.vertices.size(), 1.0))
    {}

    double go() {
        Problem p = {sphere, sphere, kernel, src_strength};

        return eval_integral_equation(p, qs, {obs_length_scale, obs_pt, obs_n});
    }
    double go_row() {
        Problem p = {sphere, sphere, kernel, src_strength};

        auto row = integral_equation_vector(p, qs, {obs_length_scale, obs_pt, obs_n});
        double row_sum = 0.0;
        for(std::size_t i = 0; i < row.size(); i++) {
            row_sum += row[i] * src_strength[i];
        }
        return row_sum;
    }

    Mesh sphere;
    QuadStrategy qs;
    Kernel kernel;
    Vec3<double> obs_pt;
    Vec3<double> obs_n;
    double obs_length_scale;
    std::vector<double> src_strength;
}; 

TEST(EvalIntegralEquationSphereSurfaceArea) {
    EvalProb ep(5, 3, 2, one);
    double result = ep.go();
    double result2 = ep.go_row();
    double exact_surf_area = 4*M_PI*9;
    CHECK_CLOSE(result, exact_surf_area, 1e-1);
    CHECK_CLOSE(result2, exact_surf_area, 1e-1);
}

TEST(ConstantLaplace) {
    EvalProb ep(5, 3, 2, laplace_double);
    double result = ep.go();
    double result2 = ep.go_row();
    CHECK_CLOSE(result, 1.0, 1e-3);
    CHECK_CLOSE(result2, 1.0, 1e-3);
}

TEST(ConstantLaplaceBoundary) {
    EvalProb ep(1, 3, 3, laplace_double);
    int n_verts = ep.sphere.vertices.size();
    // n_verts = 1;
    for (int i = 0; i < n_verts;i++) {
        ep.obs_pt = ep.sphere.vertices[i];
        ep.obs_n = -normalized(ep.obs_pt);
        double result = ep.go();
        double result2 = ep.go();
        CHECK_CLOSE(result, 1.0, 1e-2);
        CHECK_CLOSE(result2, 1.0, 1e-2);
    }
}

TEST(MatrixRowVsEval) {
    EvalProb ep(4, 3, 2, laplace_double);
    double result = ep.go();
    double result2 = ep.go_row();
    CHECK_CLOSE(result, 1.0, 1e-3);
    CHECK_CLOSE(result2, 1.0, 1e-3);
}

TEST(MassTerm) {
    Mesh sphere = clean_mesh(sphere_mesh({0,0,0}, 1.0));
    std::vector<double> str(sphere.vertices.size(), 1.0);
    str[0] = 0.0;
    Problem p = {sphere, sphere, one, str};
    QuadStrategy qs(2);
    auto res = mass_term(p, qs);
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
    Mesh sphere = refine_clean(sphere_mesh({0,0,0}, 1.0), 3);
    int n_verts = sphere.vertices.size();
    std::vector<double> str(n_verts, 1.0);

    QuadStrategy qs(2, 2, 3, 3, 3.0, 1e-2);
    Problem p = {sphere, sphere, one, str};
    auto res = direct_interact(p, qs);
    auto matrix = interact_matrix(p, qs);

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
    Mesh sphere = refine_clean(sphere_mesh({0,0,0}, 1.0), 2);
    std::vector<double> str(sphere.vertices.size(), 1.0);

    QuadStrategy qs(2, 2, 3, 3, 3.0, 1e-2);
    Problem p_double = {sphere, sphere, laplace_double, str};
    Problem p_single = {sphere, sphere, laplace_single, str};
    auto res0 = direct_interact(p_double, qs);
    auto res1 = direct_interact(p_single, qs);
    for (unsigned int i = 0; i < res1.size(); i++) {
        res1[i] = -res1[i];
    }

    Problem p_mass = {sphere, sphere, one, str};
    auto res2 = mass_term(p_mass, qs);
    CHECK_ARRAY_CLOSE(res0, res2, sphere.vertices.size(), 3e-2);
    CHECK_ARRAY_CLOSE(res1, res2, sphere.vertices.size(), 3e-2);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
