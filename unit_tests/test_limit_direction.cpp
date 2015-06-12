#include "catch.hpp"
#include "limit_direction.h"
#include "mesh.h"
#include "util.h"
#include "nearest_neighbors.h"

using namespace tbem;

TEST_CASE("nearfield facets on line", "[limit_direction]") 
{
    Facet<2> f{{{1,1},{2,1}}};
    Mesh<2> m{{f}};
    auto result = NearfieldFacetFinder<2>(m.facets, 1.0).find({0, 0});
    REQUIRE(result.facet_indices.size() == 0);
    REQUIRE(result.nearest_facet_idx == 0);
    REQUIRE(result.pt == m.facets[0][0]);
    REQUIRE(result.distance == std::sqrt(2));
}

TEST_CASE("nearfield facets out of line", "[limit_direction]") 
{
    Facet<2> f{{{-1,-1},{1,-1}}};
    Mesh<2> m{{f}};
    auto result = NearfieldFacetFinder<2>(m.facets, 1.0).find({0, 0});
    REQUIRE(result.nearest_facet_idx == 0);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0,-1}, 2, 1e-14);
    REQUIRE(result.distance == 1.0);
}

TEST_CASE("nearfield facets two edges", "[limit_direction]") 
{
    Facet<2> f{{{0.5,-1},{1,-1}}};
    Facet<2> f2{{{-1,-1},{0.5,-1}}};
    Mesh<2> m{{f,f2}};
    auto result = NearfieldFacetFinder<2>(m.facets, 1.0).find({0, 0});
    REQUIRE(result.nearest_facet_idx == 1);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0,-1}, 2, 1e-14);
    REQUIRE(result.distance == 1.0);
}

TEST_CASE("nearfield facets at intersection", "[limit_direction]") 
{
    Facet<2> f{{{0,0},{1,0}}};
    Facet<2> f2{{{0,1},{0,0}}};
    Mesh<2> m{{f,f2}};
    auto res = NearfieldFacetFinder<2>(m.facets, 1.0).find({0, 0});
    REQUIRE(res.facet_indices.size() == 2);
    bool correct_facets = (res.facet_indices[0] == 0 || res.facet_indices[1] == 0) &&
        (res.facet_indices[0] == 1 || res.facet_indices[1] == 1);
    REQUIRE(correct_facets);
    REQUIRE_ARRAY_CLOSE(res.pt, Vec<double,2>{0,0}, 2, 1e-14);
    REQUIRE(res.distance == 0.0);
}

TEST_CASE("nearfield facets non-one far threshold", "[limit_direction]")
{
    Mesh<2> m{{{{{0, -1}, {1, -1}}}, {{{0, 1}, {0, 0}}}, {{{0, 5}, {0, 3}}}}};
    auto result = NearfieldFacetFinder<2>(m.facets, 3.0).find({0, 3});
    REQUIRE(result.facet_indices.size() == 2);
}

TEST_CASE("limit direction isolated edge", "[limit_direction]") 
{
    Facet<2> f{{{0, 0}, {0, 1}}};
    for (size_t i = 0; i < 100; i++) {
        auto y_val = random<double>(0, 1);
        Vec<double,2> p{0, y_val};
        auto dir = decide_limit_dir(p, {{0}, 0, p, 0.0}, {f}, 0.5);
        REQUIRE(dir == (Vec<double,2>{-0.5, (-y_val / 2) + 0.25}));
    }
}


TEST_CASE("limit direction intersection", "[limit_direction]") 
{
    Facet<2> f{{{0, 0}, {1, 0}}};
    Facet<2> f2{{{0, 1}, {0, 0}}};
    Vec<double,2> p{0, 0};
    
    auto dir = decide_limit_dir(p, {{0, 1}, 0, p, 0.0}, {f, f2}, 0.5);
    REQUIRE(dir == (Vec<double,2>{0.25, 0.5}));

    auto dir2 = decide_limit_dir(p, {{0, 1}, 1, p, 0.0}, {f, f2}, 0.5);
    REQUIRE(dir2 == (Vec<double,2>{0.5, 0.25}));
}

TEST_CASE("limit direction reflex angle", "[limit_direction]")
{
    Facet<2> f{{{1, 0}, {0, 0}}};
    Facet<2> f2{{{0, 0}, {-1, 1}}};

    Vec<double,2> p{-0.001, 0.001};
    auto dir = decide_limit_dir(p, {{0, 1}, 1, p, 0}, {f, f2}, 0.5);
    REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{-0.7495, -0.2505}, 2, 1e-12);
}

TEST_CASE("limit direction distance cutoff", "[limit_direction]")
{
    Facet<2> f{{{1, 0}, {0, 0}}};

    SECTION("zero distance") {
        Vec<double,2> p{0.5, 0};
        auto dir = decide_limit_dir(p, {{0}, 0, p, 0}, {f}, 0.5);
        REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, -0.5}, 2, 1e-12);
    }

    SECTION("very small distance") {
        Vec<double,2> p{0.5, -1e-12};
        auto dir = decide_limit_dir(p, {{0}, 0, p, 1e-12}, {f}, 0.5, 1e-11);
        REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, -0.5}, 2, 1e-10);
    }

    SECTION("large distance") {
        Vec<double,2> p{0.5, 1};
        auto dir = decide_limit_dir(p, {{0}, 0, p, 1}, {f}, 0.5);
        REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, 0}, 2, 1e-12);
    }
}


TEST_CASE("limit direction acute angle", "[limit_direction]")
{
    Facet<2> f{{{0, 0}, {1, 0}}};
    Facet<2> f2{{{1, 0.1}, {0, 0}}};

    Vec<double,2> p{0.5, 0};
    auto dir = decide_limit_dir(p, {{0, 1}, 0, p, 0}, {f, f2}, 0.5);
    REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, 0.0125}, 2, 1e-12);
}

TEST_CASE("limit direction radius is equal to separation", "[limit_direction]")
{

    Facet<2> f{{{0, 0}, {1, 0}}};
    Facet<2> f2{{{0, 0.5}, {1, 0.5}}};

    Vec<double,2> p{0.5, 0};
    auto dir = decide_limit_dir(p, {{0, 1}, 0, p, 0}, {f, f2}, 0.5);
    REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{0, 0.125}, 2, 1e-12);
}

TEST_CASE("limit direction no information", "[limit_direction]")
{
    Vec<double,2> p{0.5, 0};
    auto dir = decide_limit_dir(p, {{}, 0, p, 0}, {}, 0.5);
    REQUIRE_ARRAY_CLOSE(dir, zeros<Vec<double,2>>::make(), 3, 1e-12);
}

TEST_CASE("limit direction at element junction", "[limit_direction]")
{
    Facet<2> f{{{0, 0}, {1, 0}}};
    Facet<2> f2{{{0.5, 0.5}, {0, 0}}};
    Vec<double,2> p{0, 0};
    
    auto dir = decide_limit_dir(p, {{0, 1}, 0, p, 0.0}, {f, f2}, 0.5);
    std::cout << dir << std::endl;
    REQUIRE(dir == (Vec<double,2>{0.25, 0.5}));

    auto dir2 = decide_limit_dir(p, {{0, 1}, 1, p, 0.0}, {f, f2}, 0.5);
    std::cout << dir2 << std::endl;
    REQUIRE(dir2 == (Vec<double,2>{0.5, 0.25}));
}
