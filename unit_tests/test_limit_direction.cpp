#include "catch.hpp"
#include "limit_direction.h"
#include "mesh.h"
#include "util.h"
#include "nearest_neighbors.h"

using namespace tbem;


TEST_CASE("limit direction isolated edge", "[limit_direction]") 
{
    Facet<2> f{{{0, 0}, {0, 1}}};
    for (size_t i = 0; i < 100; i++) {
        auto y_val = random<double>(0, 1);
        Vec<double,2> p{0, y_val};
        auto dir = decide_limit_dir(p, {{}, f, p, 0.0});
        REQUIRE(dir == (Vec<double,2>{-1, -y_val + 0.5}));
    }
}


TEST_CASE("limit direction intersection", "[limit_direction]") 
{
    Facet<2> f{{{0, 0}, {1, 0}}};
    Facet<2> f2{{{0, 1}, {0, 0}}};
    Vec<double,2> p{0, 0};
    
    auto dir = decide_limit_dir(p, {{f2}, f, p, 0.0});
    REQUIRE(dir == (Vec<double,2>{0.5, 1.0}));

    auto dir2 = decide_limit_dir(p, {{f}, f2, p, 0.0});
    REQUIRE(dir2 == (Vec<double,2>{1.0, 0.5}));
}

TEST_CASE("limit direction reflex angle", "[limit_direction]")
{
    Facet<2> f{{{1, 0}, {0, 0}}};
    Facet<2> f2{{{0, 0}, {-1, 1}}};

    Vec<double,2> p{-0.001, 0.001};
    auto dir = decide_limit_dir(p, {{f}, f2, p, 0});
    REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{-1.499, -0.501}, 2, 1e-12);
}

TEST_CASE("limit direction distance cutoff", "[limit_direction]")
{
    Facet<2> f{{{1, 0}, {0, 0}}};

    SECTION("zero distance") {
        Vec<double,2> p{0.5, 0};
        auto dir = decide_limit_dir(p, {{}, f, p, 0});
        REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, -1}, 2, 1e-12);
    }

    SECTION("very small distance") {
        Vec<double,2> p{0.5, -1e-10};
        auto dir = decide_limit_dir(p, {{}, f, p, 1e-10});
        //TODO: Maybe implement a cutoff distance at some point for small
        //distances with respect to the element length -- like 1e-10
        REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, 0}, 2, 1e-12);
    }

    SECTION("large distance") {
        Vec<double,2> p{0.5, 1};
        auto dir = decide_limit_dir(p, {{}, f, p, 1});
        REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, 0}, 2, 1e-12);
    }
}


TEST_CASE("limit direction acute angle", "[limit_direction]")
{
    Facet<2> f{{{0, 0}, {1, 0}}};
    Facet<2> f2{{{1, 0.1}, {0, 0}}};

    Vec<double,2> p{0.5, 0};
    auto dir = decide_limit_dir(p, {{f2}, f, p, 0});
    REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, 0.025}, 2, 1e-12);
}
