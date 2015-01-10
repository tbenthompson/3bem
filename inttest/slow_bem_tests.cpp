#include "UnitTest++.h"
#include "bem.h"
#include "quadrature.h"
#include "mesh_gen.h"
#include "laplace_kernels.h"

using namespace tbem;

TEST(ConstantLaplaceBoundary) {
    Vec3<double> center{0, 0, 0};
    double r = 3.0;
    int refine_level = 2;
    int near_eval = 4;
    int gauss_order = 4;
    auto sphere = sphere_mesh(center, r, refine_level);
    QuadStrategy<3> qs(gauss_order, gauss_order, near_eval, 2.0, 1e-3);
    std::vector<double> src_strength(3 * sphere.facets.size(), 1.0);
    double obs_length_scale = get_len_scale<3>(sphere, 0, gauss_order);
    for (auto f: sphere.facets) {
        for (auto v: f.vertices) {
            auto obs_pt = v;
            auto obs_n = -normalized(obs_pt);
            auto p = make_problem<3>(sphere, sphere, LaplaceDouble<3>(), src_strength);

            double result = eval_integral_equation(p, qs, 
                {obs_length_scale, obs_pt, obs_n, obs_n});
            CHECK_CLOSE(result, -1.0, 1e-2);
        }
    }
}

TEST(DirectInteractConstantLaplace) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 2);
    int n_dofs = 3 * sphere.facets.size();
    std::vector<double> str(n_dofs, 1.0);

    QuadStrategy<3> qs(2, 2, 4, 3.0, 1e-3);
    auto p_double = make_problem<3>(sphere, sphere, LaplaceDouble<3>(), str);
    auto p_single = make_problem<3>(sphere, sphere, LaplaceSingle<3>(), str);
    auto res0 = direct_interact(p_double, qs);
    for (std::size_t i = 0; i < res0.size(); i++) {
        res0[i] = -res0[i];
    }
    auto res1 = direct_interact(p_single, qs);

    auto p_mass = make_problem<3>(sphere, sphere, OneKernel<3>(), str);
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
template <typename KT>
void test_one_segment2d_integration(const KT& k,
                                    std::function<double(double,double)> exact) {
    auto quad = gauss(15);
    std::array<double, 2> v0 = {-1.0, 0.0};
    std::array<double, 2> v1 = {1.0, 0.0};
    Facet<2> facet{{v0, v1}};
    auto face = FacetInfo<2>::build(facet);
    CHECK_EQUAL(face.jacobian, 1.0);
    CHECK_EQUAL(face.area_scale, 4.0);

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            double obs_x = -5.0 + 10 * (i / 19.0);
            double obs_y = -5.0 + 10 * (j / 19.0);
            Vec2<double> obs_loc = {obs_x, obs_y};
            Vec2<double> obs_normal = {0.0, 0.0};
            double result = integrate<double,1>(quad, 
                [&](const Vec<double,1> x_hat) {
                    auto eval = eval_point_influence<2>(x_hat, k, face,
                                                obs_loc, obs_normal);
                    return (eval[0] + eval[1]);
                });
                    
            double exact_val = exact(obs_x, obs_y);
            CHECK_CLOSE(result, exact_val, 1e-4);
        }
    }
}

TEST(OneSegment2D) {
    test_one_segment2d_integration(LaplaceSingle<2>(), exact_single);
    test_one_segment2d_integration(LaplaceDouble<2>(), exact_double);
}

//TODO: slow test
TEST(ConstantLaplace2D) {
    int refine = 6;
    Vec2<double> center = {20.0, 0.0};
    Mesh<2> src_circle = circle_mesh(center, 19.0, refine);
    QuadStrategy<2> qs(3, 3, 5, 3.0, 1e-3);
    std::vector<double> u(2 * src_circle.facets.size(), 7.0);
    for (double i = 1.0; i < 19.0; i++) {
        Mesh<2> obs_circle = circle_mesh(center, i, refine);
        auto p = make_problem<2>(src_circle, obs_circle, LaplaceDouble<2>(), u);

        // Do it via eval_integral_equation for each vertex.
        for (std::size_t i = 0; i < obs_circle.facets.size(); i++) {
            ObsPt<2> pt = {0.390, obs_circle.facets[i].vertices[0], {0,0}, {0,0}}; 
            double result = eval_integral_equation(p, qs, pt);
            CHECK_CLOSE(result, -7.0, 1e-4);
        }

        // Now, do all of the observation quadrature points using direct_interact
        // But, we need to scale by the length of the element of the observation
        // mesh, because the new values are for element interactions, not pt
        // interactions.
        double scaling_factor = 0.5 * dist(obs_circle.facets[0].vertices[0],
                                           obs_circle.facets[0].vertices[1]);
        auto results = direct_interact(p, qs);
        std::vector<double> all_ones(results.size(), -7.0 * scaling_factor);
        CHECK_ARRAY_CLOSE(results, all_ones, results.size(), 1e-3);
    }
}

template <int dim>
void direct_interact_one_test(const Mesh<dim>& mesh,
                              double correct) {
    int n_dofs = dim * mesh.facets.size();
    std::vector<double> str(n_dofs, 1.0);

    QuadStrategy<dim> qs(2, 2, 3, 3.0, 1e-2);
    auto p = make_problem<dim>(mesh, mesh, OneKernel<dim>(), str);
    std::vector<double> res = direct_interact(p, qs);
    auto matrix = interact_matrix(p, qs);

    std::vector<double> res2 = bem_mat_mult(matrix, OneKernel<dim>(), n_dofs, str);
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
    auto circle = circle_mesh({0,0}, 1.0, 4);
    direct_interact_one_test<2>(circle, 2 * M_PI);
}

TEST(DirectInteractOne3d) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 3);
    direct_interact_one_test<3>(sphere, 4 * M_PI);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
    // return RunOneTest("ConstantLaplace2D");
}
