#include "catch.hpp"
#include "adaptive_quad.h"

using namespace tbem;

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
