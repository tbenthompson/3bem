#include "catch.hpp"
#include "util.h"
#include "numerics.h"
#include "quadrature.h"
#include "adaptive_quad.h"

using namespace tbem;

/* A helper function for integrating a given function using a quadrature rule.
 * Via templating, can be used with 1D, 2D, double, Vec3<double> quadrature.
 */
template <typename T, size_t dim>
T integrate(const QuadRule<dim>& qr, 
    const std::function<T(std::array<double,dim>)>& fnc) 
{
    T integral_val = qr[0].w * fnc(qr[0].x_hat);;
    for (size_t i = 1; i < qr.size(); i++) {
        integral_val += qr[i].w * fnc(qr[i].x_hat);
    }
    return integral_val;
}

TEST_CASE("GaussQuadrature", "[quadrature]") 
{
    auto qr = gauss(4);
    REQUIRE(qr.size() == 4);

    SECTION("Check weights") {
        //Check that the weights are in the right range and sum properly
        double wt_sum = 0.0;
        for (auto xw: qr) {
            REQUIRE(std::fabs(xw.x_hat[0]) <= 1);
            wt_sum += xw.w;
        }
        REQUIRE_CLOSE(wt_sum, 2.0, 1e-12);
    }

    SECTION("Quadratic") {
        double result = integrate<double,1>(qr, [](std::array<double,1> x) {
                return 3 * x[0] * x[0];
            });
        REQUIRE_CLOSE(2.0, result, 1e-12);
    }
}

TEST_CASE("High order gauss quadrature", "[quadrature]") 
{
    auto qr_high_order = gauss(50);
    REQUIRE(qr_high_order.size() == 50);
    auto result = integrate<double,1>(qr_high_order, [](std::array<double,1> x) {
            return 101 * std::pow(x[0], 100);
        });
    REQUIRE_CLOSE(2.0, result, 1e-11);
}

TEST_CASE("Odd size gauss quadrature", "[quadrature]") 
{
    auto qr_odd = gauss(5);
    REQUIRE(qr_odd.size() == 5);
}

TEST_CASE("GaussExactness", "[quadrature]") {
    for (size_t i = 0; i < 10; i++) {
        size_t n = random<size_t>(0, 250);
        if (n == 0) {
            continue;
        }
        int g = 2 * n - 1;
        auto q = gauss(n);
        double result = integrate<double,1>(q, [&](std::array<double,1> x) {
                return (g + 1) * pow(x[0], g);
            });
        double exact = 2.0 * ((g + 1) % 2);
        REQUIRE(std::fabs(exact - result) < 1e-13);
    }
}

void test_tri_integrate(QuadRule<2> q2d_tri) {
    double result = integrate<double,2>(q2d_tri, [](std::array<double,2> x) {
        return std::exp(x[0] / (x[1] - 1.1));
    });
    REQUIRE_CLOSE(result, 0.337429, 1e-6);
    result = integrate<double,2>(q2d_tri, [](std::array<double,2> x) {
            return std::exp(x[0] / (x[1] + 1.1));
        });
    REQUIRE_CLOSE(result, 0.656602, 1e-6);
}

TEST_CASE("TriangleIntegrate", "[quadrature]") {
    auto q2d_tri = tri_gauss(8);
    test_tri_integrate(q2d_tri);
}

TEST_CASE("AdaptiveQuad", "[quadrature]") {
    double result = adaptive_integrate<double>([] (double x) {
            return std::exp(x);
        }, 0, 2, 1e-3);
    REQUIRE_CLOSE(result, std::exp(2) - 1, 1e-6);
}

TEST_CASE("AdaptiveQuadCloseToSingular", "[quadrature]") {
    double a = 0.5;
    double result = adaptive_integrate<double>([a] (double x) {
            return 1.0 / std::sqrt(x*x + a);
        }, -0.5, 0.5, 1e-3);
    double exact = -std::log(-1 + std::sqrt(1 + 4 * a)) +
                    std::log(1 + std::sqrt(1 + 4 * a));
    REQUIRE_CLOSE(result, exact, 1e-6);
}

TEST_CASE("AdaptiveQuadVec", "[quadrature]") {
    Vec3<double> result = adaptive_integrate<Vec3<double>>([] (double x) {
            return Vec3<double>{
                std::exp(x), std::exp(2*x), std::exp(3*x)
            };
        }, 0, 2, 1e-6);
    Vec3<double> exact = {
        std::exp(2) - 1, (std::exp(4) - 1) / 2.0, (std::exp(6) - 1) / 3.0
    };
    REQUIRE_ARRAY_CLOSE(result, exact, 3, 1e-4);
}

TEST_CASE("AdaptiveQuadTensor", "[quadrature]") {
    auto result = adaptive_integrate<Vec2<Vec2<double>>>([] (double x) {
            return Vec2<Vec2<double>>{{
                {{x, 1}}, {{0, 2 * x}} 
            }};
        }, 0, 2, 1e-6);
    Vec2<Vec2<double>> correct{{
        {{2, 2}}, {{0, 4}}
    }};
    REQUIRE_ARRAY_CLOSE((double*)(&result[0]), (double*)(&correct[0]), 4, 1e-12);
}

