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
        face({src_locs[0], src_locs[1], src_locs[2]})
    { }
    
    void go() {
        auto basis = integrate<Vec3<double>,2>(q, [&] (std::array<double,2> x_hat) {
                return eval_quad_pt<3>(x_hat, kernel, FaceInfo<3>(face),
                                       obs_loc, obs_n);
            });
        result = dot(basis, src_vals);
    }

    void check() {
        CHECK_CLOSE(result, exact, 1e-5);
    }

    Kernel<3> kernel;
    double result;
    double exact;
    std::vector<Vec3<double>> src_locs;
    Vec3<double> src_vals;
    Vec3<double> obs_loc;
    Vec3<double> obs_n;
    QuadRule<2> q;
    Facet<3> face;
};

TEST_FIXTURE(IntegrationProb, IntegralOne) {
    double abc = integrate<double,2>(q, [](std::array<double,2> x_hat){return 1.0;});
    CHECK_CLOSE(abc, 0.5, 1e-12);
    kernel = one<3>;
    exact = 1.0; go(); check();
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceSingle) {
    kernel = laplace_single;
    exact = 0.0269063; go(); check();
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceDouble) {
    kernel = laplace_double;
    exact = -0.00621003; go(); check();
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
    CHECK_CLOSE(result, 0.0269063, 1e-6);
}

//TODO: This should be refactored a bit!
struct EvalProb {
    EvalProb(int refine_level, int near_eval, int gauss_order, const Kernel<3>& k,
             Vec3<double> center = Vec3<double>{0,0,0},
             double r = 3.0):
        sphere(sphere_mesh(center,r).refine_repeatedly(refine_level)),
        qs(gauss_order, gauss_order, gauss_order, near_eval, 2.0, 1e-2),
        kernel(k),
        obs_pt(random_pt()),
        obs_n(random_pt()),
        obs_length_scale(get_len_scale(sphere, 0, gauss_order)),
        src_strength(std::vector<double>(3 * sphere.facets.size(), 1.0))
    {}

    double go() {
        Problem<3> p = {sphere, sphere, kernel, src_strength};

        return eval_integral_equation(p, qs, {obs_length_scale, obs_pt, obs_n});
    }
    double go_row() {
        Problem<3> p = {sphere, sphere, kernel, src_strength};

        auto row = integral_equation_vector(p, qs, {obs_length_scale, obs_pt, obs_n});
        double row_sum = 0.0;
        for(std::size_t i = 0; i < row.size(); i++) {
            row_sum += row[i] * src_strength[i];
        }
        return row_sum;
    }

    Mesh3D sphere;
    QuadStrategy<3> qs;
    Kernel<3> kernel;
    Vec3<double> obs_pt;
    Vec3<double> obs_n;
    double obs_length_scale;
    std::vector<double> src_strength;
}; 

TEST(EvalIntegralEquationSphereSurfaceArea) {
    EvalProb ep(5, 3, 2, one<3>);
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
    CHECK_CLOSE(result, -1.0, 1e-3);
    CHECK_CLOSE(result2, -1.0, 1e-3);
}

TEST(ConstantLaplaceBoundary) {
    EvalProb ep(1, 3, 3, laplace_double);
    for (auto f: ep.sphere.facets) {
        for (auto v: f.vertices) {
            ep.obs_pt = v;
            ep.obs_n = -normalized(ep.obs_pt);
            double result = ep.go();
            double result2 = ep.go();
            CHECK_CLOSE(result, -1.0, 1e-2);
            CHECK_CLOSE(result2, -1.0, 1e-2);
        }
    }
}

TEST(MatrixRowVsEval) {
    EvalProb ep(4, 3, 2, laplace_double);
    double result = ep.go();
    double result2 = ep.go_row();
    CHECK_CLOSE(result, -1.0, 1e-3);
    CHECK_CLOSE(result2, -1.0, 1e-3);
}

TEST(MassTerm) {
    auto sphere = sphere_mesh({0,0,0}, 1.0);
    std::vector<double> str(3 * sphere.facets.size(), 1.0);
    for (std::size_t i = 0; i < sphere.facets.size(); i++) {
        for (int d = 0; d < 3; d++) {
            if (sphere.facets[i].vertices[d][0] > 0.5) {
                str[3 * i + d] = 0.0;
            }
        }
    }
    Problem<3> p = {sphere, sphere, one<3>, str};
    QuadStrategy<3> qs(2);
    auto res = mass_term(p, qs);
    CHECK_EQUAL(res.size(), 3 * sphere.facets.size());
    double true_area = 0.0;
    for (auto f: sphere.facets) {
        true_area += tri_area(f.vertices);
    }
    double mass_area = 0.0;
    for (auto r: res) {
        mass_area += r;
    }   
    CHECK_CLOSE(mass_area, (5.0 / 6.0) * true_area, 1e-12);
}


//TODO: Fixture for this and the next one.
TEST(DirectInteractConstantLaplace) {
    auto sphere = sphere_mesh({0,0,0}, 1.0).refine_repeatedly(2);
    int n_dofs = 3 * sphere.facets.size();
    std::vector<double> str(n_dofs, 1.0);

    QuadStrategy<3> qs(2, 2, 3, 4, 3.0, 1e-3);
    Problem<3> p_double = {sphere, sphere, laplace_double, str};
    Problem<3> p_single = {sphere, sphere, laplace_single, str};
    auto res0 = direct_interact(p_double, qs);
    for (std::size_t i = 0; i < res0.size(); i++) {
        res0[i] = -res0[i];
    }
    auto res1 = direct_interact(p_single, qs);

    Problem<3> p_mass = {sphere, sphere, one<3>, str};
    auto res2 = mass_term(p_mass, qs);
    CHECK_ARRAY_CLOSE(res0, res2, n_dofs, 3e-2);
    CHECK_ARRAY_CLOSE(res1, res2, n_dofs, 3e-2);
    CHECK_ARRAY_CLOSE(res0, res1, n_dofs, 3e-2);
}

/* Two dimensional test cases */
double exact_single(double obs_x, double obs_y) {
    return 0.0795775 * (
        -2 * obs_y * atan((-1 + obs_x) / obs_y) +
        2 * obs_y * atan((1 + obs_x) / obs_y) -
        (-1 + obs_x) * (-2 + log(pow((1 - obs_x), 2) + pow(obs_y, 2))) +
        (1 + obs_x) * (-2 + log(pow((1 + obs_x), 2) + pow(obs_y, 2))));
}

double exact_double(double x, double y) {
    return -(atan((1 - x) / y) + atan((1 + x) / y)) / (2 * M_PI);
}

/* Check that integrals over a single element are being properly
 * computed in 2D. Do this for many different observation points and two
 * different kernels (laplace_single, laplace_double)
 */
TEST(OneSegment2D) {
    std::array<double, 2> v0 = {-1.0, 0.0};
    std::array<double, 2> v1 = {1.0, 0.0};

    auto quad = gauss(15);
    std::vector<std::function<double (double, double)>> exact =
        {exact_single, exact_double};
    std::vector<Kernel<2>> kernel = {laplace_single2d, laplace_double2d};
    Facet<2> facet{{v0, v1}};
    FaceInfo<2> face(facet);
    CHECK_EQUAL(face.jacobian, 1.0);
    CHECK_EQUAL(face.area, 2.0);

    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < 20; i++) {
            for (int j = 0; j < 20; j++) {
                double obs_x = -5.0 + 10 * (i / 19.0);
                double obs_y = -5.0 + 10 * (j / 19.0);
                Vec2<double> obs_loc = {obs_x, obs_y};
                Vec2<double> obs_normal = {0.0, 0.0};
                double result = integrate<double,1>(quad, 
                    [&](const Vec<double,1> x_hat) {
                        auto eval = eval_quad_pt<2>(x_hat, kernel[k], face,
                                                    obs_loc, obs_normal);
                        return (eval[0] + eval[1]);
                    });
                        
                double exact_val = exact[k](obs_x, obs_y);
                CHECK_CLOSE(result, exact_val, 1e-4);
            }
        }
    }
}

TEST(ConstantLaplace2D) {
    int refine = 5;
    Vec2<double> center = {20.0, 0.0};
    Mesh<2> src_circle = circle_mesh(center, 19.0).refine_repeatedly(refine);
    QuadStrategy<2> qs(3, 3, 3, 5, 3.0, 1e-3);
    std::vector<double> u(2 * src_circle.facets.size(), 1.0);
    for (double i = 1.0; i < 18.0; i++) {
        Mesh<2> obs_circle = circle_mesh(center, i).refine_repeatedly(refine);
        Problem<2> p{src_circle, obs_circle, laplace_double2d, u};

        // Do it via eval_integral_equation for each vertex.
        for (std::size_t i = 0; i < obs_circle.facets.size(); i++) {
            ObsPt<2> pt = {0.390, obs_circle.facets[i].vertices[0], {0,0}}; 
            double result = eval_integral_equation(p, qs, pt);
            CHECK_CLOSE(-result, 1.0, 1e-4);
        }

        // Now, do all of the observation quadrature points using direct_interact
        // But, we need to scale by the length of the element of the observation
        // mesh, because the new values are for element interactions, not pt
        // interactions.
        double scaling_factor = 0.5 * dist(obs_circle.facets[0].vertices[0],
                                           obs_circle.facets[0].vertices[1]);
        auto results = direct_interact(p, qs);
        std::vector<double> all_ones(results.size(), -scaling_factor);
        CHECK_ARRAY_CLOSE(results, all_ones, results.size(), 1e-4);
    }
}

template <int dim>
void direct_interact_one_test(const Mesh<dim>& mesh,
                              double correct) {
    int n_dofs = dim * mesh.facets.size();
    std::vector<double> str(n_dofs, 1.0);

    QuadStrategy<dim> qs(2, 2, 3, 3, 3.0, 1e-2);
    Problem<dim> p = {mesh, mesh, one<dim>, str};
    auto res = direct_interact(p, qs);
    auto matrix = interact_matrix(p, qs);

    auto res2 = bem_mat_mult(matrix, n_dofs, str);
    double total = 0.0;
    for (auto r: res) {
        total += r;
    }
    double total2 = 0.0;
    for (auto r2: res2) {
        total2 += r2;
    }
    double sa2 = std::pow(correct, 2);
    double error = std::fabs((sa2 - total) / sa2);
    double error2 = std::fabs((sa2 - total2) / sa2);
    CHECK_CLOSE(error, 0.0, 6.4e-2);
    CHECK_CLOSE(error2, 0.0, 6.4e-2);
}

TEST(DirectInteractOne2d) {
    auto circle = circle_mesh({0,0}, 1.0).refine_repeatedly(4);
    direct_interact_one_test(circle, 2 * M_PI);
}

TEST(DirectInteractOne3d) {
    auto sphere = sphere_mesh({0,0,0}, 1.0).refine_repeatedly(3);
    direct_interact_one_test(sphere, 4 * M_PI);
}


int main(int, char const *[])
{
    return UnitTest::RunAllTests();
    // return RunOneTest("ConstantLaplace2D");
}
