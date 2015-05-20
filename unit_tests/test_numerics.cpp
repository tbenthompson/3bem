#include "catch.hpp"
#include "numerics.h"
#include "numbers.h"
#include <iostream>
#include "util.h"

using namespace tbem;


TEST_CASE("LinearMapping", "[numerics]") {
    for (size_t i = 0; i < 100; i++) {
        double real = random<double>(0, 100);
        double v0 = random<double>(0, 100);
        double v1 = random<double>(0, 100);
        double ref = real_to_ref(real, v0, v1);
        std::array<Vec2<double>,2> locs = {{{v0,0}, {v1,0}}};
        double real2 = ref_to_real({ref}, locs)[0];
        REQUIRE(std::fabs(real - real2) < 1e-8);
    }
}

TEST_CASE("AreaTri", "[numerics]") {
    double result = tri_area({{{0,0,0},{1,0,0},{0,1,0}}});
    REQUIRE_CLOSE(result, 0.5, 1e-12);
    result = tri_area({{{1,1,1},{3,1,1},{3,3,1}}});
    REQUIRE_CLOSE(result, 2.0, 1e-12);
}


TEST_CASE("LinearInterp", "[numerics]") {
    REQUIRE_CLOSE(linear_interp<3>({0,0},{1,0,0}), 1.0, 1e-12);
    REQUIRE_CLOSE(linear_interp<3>({1,0},{0,1,0}), 1.0, 1e-12);
    REQUIRE_CLOSE(linear_interp<3>({0,1},{0,0,1}), 1.0, 1e-12);
    REQUIRE_CLOSE(linear_interp<3>({0.5,0.5},{0,0,1}), 0.5, 1e-12);
    REQUIRE_CLOSE(linear_interp<3>({0.0,0.5},{0,0,1}), 0.5, 1e-12);
}

TEST_CASE("LinearInterpOnes", "[numerics]") {
    for (size_t i = 0; i < 30; i++) {
        double x_hat = random<double>(0, 1);
        double y_hat = random<double>(0, 1);
        double result = linear_interp<3>({x_hat, y_hat}, {1,1,1});
        double exact = 1.0;
        REQUIRE(std::fabs(result - exact) < 1e-12);
    }
}
