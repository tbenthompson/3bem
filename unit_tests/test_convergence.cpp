#include "catch.hpp"
#include "integral_term.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "util.h"
#include "gte_wrapper.h"
#include "mesh_gen.h"
#include "integral_operator.h"
#include "continuity_builder.h"
#include "constraint_matrix.h"

/*
 * This file has tests at a slightly higher level that make sure that the 
 * quadrature methods converge sufficiently quickly for a variety of pieces:
 *
 * nearfield observation quadrature (FAILURE!)
 * farfield observation quadrature 
 * far field inner integral quadrature 
 * near field inner integral quadrature -- corner 
 * near field inner integral quadrature -- noncorner 
 * singular inner integral quadrature 
 */
//TODO: These tests are slow.
//TODO: Make these tests faster by refining the mesh less.

using namespace tbem;

template <size_t dim>
Vec<Vec<double,dim>,dim> unit_facet() 
{
    auto facet = zeros<Vec<Vec<double,dim>,dim>>::make();
    for (size_t d = 0; d < dim - 1; d++) {
        facet[d + 1][d] = 1.0;
    }
    return facet;
}

template <typename F2>
std::vector<double> convergence_test(size_t steps, const F2& calc_fnc) 
{
    std::vector<double> out(steps, 0.0);
    std::vector<double> prev_result;
    for (size_t i = 0; i < steps; i++) {
        auto result = calc_fnc(i);
        if (i > 0) {
            std::vector<double> diff(result.size());
            for (size_t j = 0; j < result.size(); j++) {
                if (std::fabs(result[j]) < 1e-10) {
                    diff[j] = 0.0;
                } else {
                    diff[j] = fabs((result[j] - prev_result[j]) / result[j]);
                }
            }
            auto max = *std::max_element(diff.begin(), diff.end());
            out[i] = max;
        }
        prev_result = result;
    }
    return out;
}

template <typename F2>
std::vector<double> integral_convergence_test(
    size_t steps, const Vec<double,2>& loc, const F2& calc_fnc) 
{
    return convergence_test(steps, [&](size_t step)
        {
            auto facet = unit_facet<2>();
            auto facet_info = FacetInfo<2>::build(facet);
            ObsPt<2> obs{loc, facet_info.normal, facet_info.normal * 0.3};
            IntegralTerm<2,1,1> term{obs, facet_info};

            auto result_big = calc_fnc(term, step);
            return std::vector<double>{result_big[0][0][0], result_big[1][0][0]};
        }
    );
}

void convergent_sinh(double x_loc) 
{
    for (size_t i = 0; i < 7; i++) {
        double y_dist = std::pow(0.1, static_cast<double>(i));
        auto error = integral_convergence_test(6, {x_loc, y_dist},
            [] (IntegralTerm<2,1,1>& term, size_t step) 
            {
                LaplaceHypersingular<2> K;
                size_t sinh_order = 2 + 2 * step; 
                SinhIntegrator<2,1,1> integrator(sinh_order);
                auto nearest_pt = closest_pt_facet(term.obs.loc, term.src_face.facet);
                return integrator.compute_nearfield(K, term, nearest_pt);
            }
        );
        REQUIRE(error.back() < std::pow(0.1, 1e-6));
    }
}

TEST_CASE("sinh convergence 2D non-corner", "[convergence]")
{
    convergent_sinh(0.6);
}

TEST_CASE("sinh convergence 2D corner", "[convergence]")
{
    convergent_sinh(0.98);
}

TEST_CASE("farfield gauss convergence", "[convergence]") 
{
    double dist = 3.0;
    auto error = integral_convergence_test(5, {dist, dist},
        [] (IntegralTerm<2,1,1>& term, size_t step) {
            LaplaceHypersingular<2> K;
            auto integrator = make_sinh_integrator(1, 1, 1, step + 1, 1, 0.0, K);
            return integrator.compute_farfield(term);
        }
    );
    REQUIRE(error.back() < 1e-4);
}

TEST_CASE("inner integrals limit steps", "[convergence]")
{
    auto xs = linspace(0.05, 0.95, 6);
    std::vector<double> ys{0.0, 0.001, 0.1, 0.5, 1.0, 2.0, 5.0};
    for (size_t i = 0; i < xs.size(); i++) {
        for (size_t j = 0; j < ys.size(); j++) {
            auto error = integral_convergence_test(11, {xs[i], ys[j]},
                [] (IntegralTerm<2,1,1>& term, size_t step) {
                    LaplaceHypersingular<2> K;
                    auto integrator = make_sinh_integrator(
                        step + 2, 1, 1, step + 3, step + 2, 3.0, K
                    );
                    return integrator.compute_term(term);
                }
            );
            // std::cout << xs[i] << " " << ys[j] << " " << error.back() << std::endl;
            REQUIRE(error.back() < 1e-8);
        }
    }
}

template <typename F>
std::vector<double> operator_convergence_test(size_t steps, Mesh<2>& obs_mesh,
    Mesh<2>& src_mesh, const F& integrator_builder)
{
    //Applying constraints is essential because otherwise
    auto obs_continuity = from_constraints(
        convert_to_constraints(mesh_continuity(obs_mesh.begin()))
    );
    auto src_continuity = from_constraints(
        convert_to_constraints(mesh_continuity(src_mesh.begin()))
    );
    auto error = convergence_test(steps, [&](size_t step)
        {
            auto integrator = integrator_builder(step);
            auto op = dense_boundary_operator(
                obs_mesh, src_mesh, integrator, src_mesh
            );
            auto condensed_op = condense_matrix(obs_continuity, src_continuity, op);
            return condensed_op.data();
        }
    );
    return error;
}

TEST_CASE("observation quadrature convergence", "[convergence]")
{
    //To see convergence in this test, I had to make two modifications to
    //the original quadrature parameters. First, I had to ensure continuity.
    //Otherwise, there are singularities in the integrand of the outer
    //observation point integral. 
    //Second, I had to increase the source quadrature orders. Otherwise, 
    //convergence plateaued at ~1e-5 error. 
    LaplaceHypersingular<2> K;
    auto m = circle_mesh({0, 0}, 1.0, 0); 
    auto error = operator_convergence_test(8, m, m,
        [&] (size_t step) 
        {
            return make_sinh_integrator(12, 4 * step + 3, 2, 10, 8, 3.0, K); 
        }
    );
    for (auto e: error) {
        std::cout << e << std::endl;
    }
    //TODO: Nearfield quadrature convergences sufficiently assuming continuity.
    //But, it does so slowly. How can it be sped up? What about a tanh-sinh 
    //quadrature rule? At the very least, nearfield observation order should
    //be separated from farfield observation order.
    REQUIRE(error.back() < 5e-5);
}

TEST_CASE("observation quadrature convergence farfield", "[convergence]")
{
    LaplaceHypersingular<2> K;
    double sep = 4.0;
    auto m1 = circle_mesh({0, 0}, 1.0, 0); 
    auto m2 = circle_mesh({sep + 2.0, 0}, 1.0, 0); 
    auto error = operator_convergence_test(4, m1, m2,
        [&] (size_t step) 
        {
            return make_sinh_integrator(12, 2, 4 * step + 3, 2, 8, 3.0, K); 
        }
    );
    for (size_t i = 0; i < error.size(); i++) {
        REQUIRE(error[i] < 5e-5);
    }
}
