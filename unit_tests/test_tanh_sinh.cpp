#include "catch.hpp"
#include "tanh_sinh.h"
#include "integrate.h"
#include "numbers.h"

using namespace tbem;

TEST_CASE("Tanh-sinh simple", "[tanhsinh]")
{
    size_t nh = 10;
    for (size_t nq = 5; nq < 30; nq += 5) {
        auto q = tanh_sinh(nq);
        auto result = integrate<double,1>(q, [] (std::array<double,1> pt) {
            return std::pow(std::sin(pt[0]), 2);
        });
        auto exact = 1.0 - (std::sin(2.0) / 2.0);
        std::cout << nq << " " << (result - exact) << std::endl;
        // REQUIRE_CLOSE(result, exact, 1e-10);
    }
}

TEST_CASE("Tanh-sinh endpoint singular", "[tanhsinh]")
{
    auto q = tanh_sinh(20);
    auto result = integrate<double,1>(q, [] (std::array<double,1> pt) {
        return 1.0 / std::sqrt(1.0 - pt[0] * pt[0]);
    });
    std::cout << (result - M_PI) << std::endl;
    REQUIRE_CLOSE(result, M_PI, 1e-10);
}
