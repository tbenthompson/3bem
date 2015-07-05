#include "catch.hpp"
#include "mesh_preprocess.h"

using namespace tbem;

TEST_CASE("two intersecting segments", "[mesh_preprocess]")
{
    Vec<Vec<double,2>,2> f0{{{0, 0}, {1, 0}}};
    Vec<Vec<double,2>,2> f1{{{0.5, -1}, {0.5, 1}}};
    auto intersections = MeshPreprocessor<2>().find_intersections({f0}, {f1});
    REQUIRE(intersections.size() == 1);
    REQUIRE(intersections[0].facet_idx_A == 0);
    REQUIRE(intersections[0].facet_idx_B == 0);
    REQUIRE(intersections[0].pts.size() == 1);
    REQUIRE(intersections[0].pts[0] == (Vec<double,2>{0.5, 0}));
}

TEST_CASE("three intersections", "[mesh_preprocess]")
{
    Vec<Vec<double,2>,2> f0{{{0, 0}, {1, 0}}};
    Vec<Vec<double,2>,2> f1{{{0.6, 0.5}, {1, 0.5}}};
    Vec<Vec<double,2>,2> f2{{{0.5, -1}, {0.5, 1}}};
    Vec<Vec<double,2>,2> f3{{{0.6, -1}, {0.9, 1}}};
    auto intersections = MeshPreprocessor<2>().find_intersections({f0, f1}, {f2, f3});
    REQUIRE(intersections.size() == 3);
}

TEST_CASE("touching intersection", "[mesh_preprocess]")
{
    Vec<Vec<double,2>,2> f0{{{0, 0}, {1, 0}}};
    Vec<Vec<double,2>,2> f1{{{1, 0}, {2, 0}}};
    auto intersections = MeshPreprocessor<2>().find_intersections({f0}, {f1});
    REQUIRE(intersections.size() == 1);
}

TEST_CASE("two intersecting triangles", "[mesh_preprocess]")
{
    Vec<Vec<double,3>,3> f0{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> f1{{{0.5, 0.5, 0.0}, {0, 0, 0}, {0, 0, 1}}};
    auto intersections = MeshPreprocessor<3>().find_intersections({f0}, {f1});
    REQUIRE(intersections.size() == 1);
    REQUIRE(intersections[0].facet_idx_A == 0);
    REQUIRE(intersections[0].facet_idx_B == 0);
    REQUIRE(intersections[0].pts.size() == 2);
    REQUIRE(intersections[0].pts[0] == (Vec<double,3>{0.0, 0.0, 0.0}));
    REQUIRE(intersections[0].pts[1] == (Vec<double,3>{0.5, 0.5, 0.0}));
}

TEST_CASE("split at intersection 2D", "[mesh_preprocess]")
{
    Vec<Vec<double,2>,2> f0{{{0, 0}, {1, 0}}};
    Vec<Vec<double,2>,2> f1{{{0.5, -1}, {0.5, 1}}};
    MeshPreprocessor<2> mp;
    auto intersections = mp.find_intersections({f0}, {f1});
    auto split_facets = mp.split_facets_at_intersections({f0}, intersections);
    REQUIRE(split_facets.size() == 2);
}

TEST_CASE("don't split at end points 2D", "[mesh_preprocess]")
{
    Vec<Vec<double,2>,2> f0{{{0, 0}, {1, 0}}};
    Vec<Vec<double,2>,2> f1{{{1.0, -1}, {1.0, 1}}};
    MeshPreprocessor<2> mp;
    auto intersections = mp.find_intersections({f0}, {f1});
    auto split_facets = mp.split_facets_at_intersections({f0}, intersections);
    REQUIRE(split_facets.size() == 1);

    auto intersections_swapped = mp.find_intersections({f1}, {f0});
    auto split_facets_swapped = mp.split_facets_at_intersections({f1}, intersections);
    REQUIRE(split_facets_swapped.size() == 2);
}
