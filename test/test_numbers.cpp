#include "catch.hpp"
#include "numbers.h"

using namespace tbem;

TEST_CASE("Range", "[numbers]") {
    auto nats5 = range(-1, 8);
    double correct[] = {-1, 0, 1, 2, 3, 4, 5, 6, 7};
    REQUIRE_ARRAY_EQUAL(nats5, correct, 9);
}

TEST_CASE("Linspace", "[numbers]") {
    auto vals = linspace(0.0, 1.0, 10);
    REQUIRE(vals.size() == 10);
    REQUIRE(vals[0] == 0.0);
    REQUIRE(vals[vals.size() - 1] == 1.0);
}
