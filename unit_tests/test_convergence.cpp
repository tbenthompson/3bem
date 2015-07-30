#include "catch.hpp"
#include "integral_term.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "util.h"
#include "gte_wrapper.h"
#include "mesh_gen.h"
#include "integral_operator.h"

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
            auto integrator = make_sinh_integrator(1, 1, step + 1, 1, 0.0, K);
            return integrator.compute_farfield(term);
        }
    );
    REQUIRE(error.back() < 1e-4);
}

TEST_CASE("inner integrals limit steps", "[convergence]")
{
    auto xs = linspace(0.05, 0.95, 10);
    for (size_t i = 0; i < xs.size(); i++) {
        std::vector<double> ys{0.0, 0.0001, 0.001, 0.01, 0.1, 0.5};
        for (size_t j = 0; j < ys.size(); j++) {
            auto error = integral_convergence_test(10, {xs[i], ys[j]},
                [] (IntegralTerm<2,1,1>& term, size_t step) {
                    LaplaceHypersingular<2> K;
                    auto integrator = make_sinh_integrator(12, 1, 1, step + 2, 3.0, K);
                    return integrator.compute_term(term);
                }
            );
            std::cout << xs[i] << " " << ys[j] << " " << error.back() << std::endl;
            // REQUIRE(error.back() < 1e-8);
        }
    }
}

template <typename F>
std::vector<double> operator_convergence_test(size_t steps, Mesh<2>& obs_mesh,
    Mesh<2>& src_mesh, const F& integrator_builder)
{
    auto error = convergence_test(steps, [&](size_t step)
        {
            auto integrator = integrator_builder(step);
            auto op = dense_boundary_operator(
                obs_mesh, src_mesh, integrator, src_mesh
            );
            return op.data();
        }
    );
    return error;
}

TEST_CASE("observation quadrature convergence", "[convergence]")
{
    LaplaceHypersingular<2> K;
    auto m = circle_mesh({0, 0}, 1.0, 3); 
    auto error = operator_convergence_test(4, m, m,
        [&] (size_t step) 
        {
            return make_sinh_integrator(12, 4 * step + 3, 2, 8, 3.0, K); 
        }
    );
    //TODO: FIX THE NEARFIELD OBSERVATION QUADRATURE!
    for (size_t i = 0; i < error.size(); i++) {
        REQUIRE(error[i] < 2e-6);
    }
}

TEST_CASE("observation quadrature convergence farfield", "[convergence]")
{
    LaplaceHypersingular<2> K;
    double sep = 4.0;
    auto m1 = circle_mesh({0, 0}, 1.0, 3); 
    auto m2 = circle_mesh({sep + 2.0, 0}, 1.0, 3); 
    auto error = operator_convergence_test(4, m1, m2,
        [&] (size_t step) 
        {
            return make_sinh_integrator(12, 4 * step + 3, 2, 8, 3.0, K); 
        }
    );
    for (size_t i = 0; i < error.size(); i++) {
        REQUIRE(error[i] < 2e-6);
    }
}
