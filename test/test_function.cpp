#include "catch.hpp"
#include "vectorx.h"

using namespace tbem;

TEST_CASE("VectorX", "[function]") 
{
    BlockVectorX a{{1,2}, {3,4}};
    BlockVectorX b{{-1,-2}, {-3,-4}};
    BlockVectorX c{{0,0}, {0,0}};

    SECTION("Size") {
        REQUIRE(a.size() == 2);
        REQUIRE(a[0].size() == 2);
    }

    SECTION("Index") {
        REQUIRE(a[0][0] == 1.0);
    }

    SECTION("Add") {
        a += b;
        REQUIRE(a == c);
    }

    SECTION("Subtract and negate") {
        a -= -b;
        REQUIRE(a == c);
    }

    SECTION("Mul") {
        a *= 2;
        a += b * 2;
        REQUIRE(a == c);
    }

    SECTION("More complex math") {
        auto d = a * b + 2 - (-c);
        BlockVectorX correct{{1, -2}, {-7, -14}};
        REQUIRE(d == correct);
    }
}
