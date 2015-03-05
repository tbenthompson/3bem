#include "UnitTest++.h"
#include "dense_builder.h"
#include "quadrature.h"
#include "mesh_gen.h"
#include "laplace_kernels.h"
#include "test_shared.h"
#include "vectorx.h"

using namespace tbem;

TEST(ConstantLaplaceBoundary) {
    Vec3<double> center{0, 0, 0};
    double r = 3.0;
    int refine_level = 2;
    int near_eval = 4;
    int gauss_order = 4;
    auto sphere = sphere_mesh(center, r, refine_level);
    QuadStrategy<3> qs(gauss_order, gauss_order, near_eval, 2.0, 1e-3);
    std::vector<double> src_strength(sphere.n_dofs(), 1.0);
    double obs_length_scale = std::sqrt(tri_area(sphere.facets[0]));
    for (auto f: sphere.facets) {
        for (auto v: f) {
            auto obs_pt = v;
            auto obs_n = -normalized(obs_pt);
            LaplaceDouble<3> double_kernel;
            auto p = make_boundary_integral<3>(sphere, sphere, double_kernel);
            ObsPt<3> obs{obs_length_scale, obs_pt, obs_n, obs_n};
            auto op = mesh_to_point_operator(p, qs, obs);
            auto result = op.apply({src_strength})[0];
            CHECK_CLOSE(result[0], -1.0, 1e-2);
        }
    }
}

TEST(GalerkinMatrixConstantLaplace) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 2);
    int n_dofs = sphere.n_dofs();
    std::vector<double> str(n_dofs, 1.0);

    QuadStrategy<3> qs(2, 2, 4, 3.0, 1e-3);
    LaplaceDouble<3> double_kernel;
    auto p_double = make_boundary_integral<3>(sphere, sphere, double_kernel);
    LaplaceSingle<3> single_kernel;
    auto p_single = make_boundary_integral<3>(sphere, sphere, single_kernel);
    auto mat0 = mesh_to_mesh_operator(p_double, qs);
    auto res0 = mat0.apply({str})[0];
    for (std::size_t i = 0; i < res0.size(); i++) {
        res0[i] = -res0[i];
    }
    auto mat1 = mesh_to_mesh_operator(p_single, qs);
    auto res1 = mat1.apply({str})[0];

    IdentityScalar<3> identity;
    auto p_mass = make_boundary_integral<3>(sphere, sphere, identity);
    auto mass_op = mass_operator(p_mass, qs);
    auto res2 = mass_op.apply({str})[0];
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
            auto result = integrate<Vec1<Vec1<double>>,1>(quad, 
                [&](const Vec<double,1> x_hat) {
                    auto eval = eval_point_influence<2>(x_hat, k, face,
                                                obs_loc, obs_normal);
                    return (eval[0] + eval[1]);
                });
                    
            double exact_val = exact(obs_x, obs_y);
            CHECK_CLOSE(result[0][0], exact_val, 1e-4);
        }
    }
}

TEST(OneSegment2D) {
    LaplaceSingle<2> single_kernel;
    test_one_segment2d_integration(single_kernel, exact_single);
    LaplaceDouble<2> double_kernel;
    test_one_segment2d_integration(double_kernel, exact_double);
}

TEST(ConstantLaplace2D) {
    int refine = 6;
    Vec2<double> center = {20.0, 0.0};
    Mesh<2> src_circle = circle_mesh(center, 19.0, refine);
    QuadStrategy<2> qs(3, 3, 5, 3.0, 1e-3);
    LaplaceDouble<2> double_kernel;
    std::vector<double> u(src_circle.n_dofs(), 7.0);
    for (double i = 1.0; i < 19.0; i++) {
        Mesh<2> obs_circle = circle_mesh(center, i, refine);
        auto p = make_boundary_integral<2>(obs_circle, src_circle, double_kernel);

        // Do it via eval_integral_equation for each vertex.
        for (std::size_t i = 0; i < obs_circle.n_facets(); i++) {
            ObsPt<2> pt = {0.390, obs_circle.facets[i][0], {0,0}, {0,0}}; 
            auto op = mesh_to_point_operator(p, qs, pt);
            double result = op.apply({u})[0][0];
            CHECK_CLOSE(result, -7.0, 1e-4);
        }

        // Now, do all of the observation quadrature points using mesh_to_mesh_operator
        // But, we need to scale by the length of the element of the observation
        // mesh, because the new values are for element interactions, not pt
        // interactions. For example, Integral(1)_{0 to 0.2} == 0.2
        auto results_op = mesh_to_mesh_operator(p, qs);
        auto results = results_op.apply({u})[0];
        std::vector<double> all_ones(results.size());
        const double integral_of_basis_fnc = 0.5;
        for (size_t j = 0; j < obs_circle.facets.size(); j++) {
            double integral_of_one = integral_of_basis_fnc * 
                dist(obs_circle.facets[j][0], obs_circle.facets[j][1]);
            double val = -7.0 * integral_of_one;
            all_ones[2 * j] = val;
            all_ones[2 * j + 1] = val;
        }
        CHECK_ARRAY_CLOSE(results, all_ones, results.size(), 1e-3);
    }
}

template <size_t dim>
void galerkin_matrix_one_test(const Mesh<dim>& mesh,
                              double correct) {
    std::vector<double> str(mesh.n_dofs(), 1.0);
    IdentityScalar<dim> identity;
    auto p = make_boundary_integral<dim>(mesh, mesh, identity);
    QuadStrategy<dim> qs(2);

    auto matrix = mesh_to_mesh_operator(p, qs);
    auto res = matrix.apply({str})[0]; 
    double total = 0.0;
    for (auto r: res) {
        total += r;
    }

    double sa_squared = std::pow(correct, 2);
    double error = std::fabs((sa_squared - total) / sa_squared);
    CHECK_CLOSE(error, 0.0, 6.4e-2);
}

TEST(GalerkinMatrixOne2d) {
    auto circle = circle_mesh({0,0}, 1.0, 4);
    galerkin_matrix_one_test<2>(circle, 2 * M_PI);
}

TEST(GalerkinMatrixOne3d) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 3);
    galerkin_matrix_one_test<3>(sphere, 4 * M_PI);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
    // return RunOneTest("ConstantLaplace2D");
}
