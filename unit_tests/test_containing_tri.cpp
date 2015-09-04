#include "catch.hpp"
#include "containing_tri.h"

using namespace tbem;

TEST_CASE("find containing tri", "[containing_tri]")
{
    std::vector<Vec<size_t,3>> tris{
        {0, 1, 2}, {2, 1, 3},
        {0, 4, 5}, {0, 5, 1},
        {1, 5, 6}, {1, 6, 3},
        {3, 6, 7}, {2, 7, 3},
        {2, 7, 4}, {2, 4, 0}
    };
    std::vector<Vec<double,2>> pts{
        {-1, -1}, {1, -1}, {-1, 1}, {1, 1}, {-2, -2}, {2, -2}, {2, 2}, {-2, 2}
    };

    SECTION("inside") {
        auto result = find_containing_tri_idx({-0.5, 0}, tris, pts);
        REQUIRE(result.first);
        REQUIRE(result.second == 0);
    }

    SECTION("outside") {
        auto result = find_containing_tri_idx({-2.5, 0}, tris, pts);
        REQUIRE(!result.first);
    }
}
