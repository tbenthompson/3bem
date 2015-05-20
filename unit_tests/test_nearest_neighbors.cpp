#include "catch.hpp"
#include "nearest_neighbors.h"
#include "util.h"
#include "mesh.h"

using namespace tbem;

TEST_CASE("SimpleIdenticalPoints", "[nearest_neighbors]") 
{
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    };
    auto oct = build_octree(es, 1);
    auto pts = nearby_points({1.0, 2.0, 0.0}, es, oct, 0.0);
    REQUIRE(pts.size() == 1);
    REQUIRE(pts[0] == 0);
}

TEST_CASE("TwoIdenticalPoints", "[nearest_neighbors]") 
{
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    };
    auto oct = build_octree(es, 1);
    auto pts = nearby_points({1.0, 2.0, 0.0}, es, oct, 0.0);
    REQUIRE(pts.size() == 2);
}

TEST_CASE("ClosePoints", "[nearest_neighbors]")
{
    double eps = 1e-13;
    std::vector<Vec<double,3>> es{
        {0.0, -2.0, 3.0}, {1.0, 2.0, 0.0}, {1.0, 2.0 - (eps / 10.), 0.0}
    };
    auto oct = build_octree(es, 1);
    auto pts = nearby_points({1.0, 2.0, 0.0}, es, oct, eps);
    REQUIRE(pts.size() == 2);
}

std::vector<Vec<double,2>> line_pts(size_t n) 
{
    std::vector<Vec<double,2>> pts;     
    for (size_t i = 0; i < n; i++) {
        pts.push_back({static_cast<double>(i), 0});
        pts.push_back({static_cast<double>(i + 1), 0});
    }
    return pts;
}

TEST_CASE("AllPairs", "[nearest_neighbors]") {
    size_t n = 500;
    auto pts = line_pts(n);
    auto oct = build_octree(pts, 50);
    auto result = nearby_points_all_pairs(pts, pts, oct, oct, 0.0);
    // 1000 pairs like (k, k). 998 pairs like (k, k+1) and (k+1, k)
    REQUIRE(result.size() == 1998);
}

TEST_CASE("AllPairsPerformance", "[nearest_neighbors]") 
{
    size_t n = 50000;
    auto pts = line_pts(n);
    auto oct = build_octree(pts, 50);
    //TODO: Capacity test
    auto result = nearby_points_all_pairs(pts, pts, oct, oct, 0.0);
}

TEST_CASE("NearestFacetsOneFacetEndpoint", "[nearest_neighbors]") {
    Facet<2> f{{{1,1},{2,1}}};
    Mesh<2> m{{f}};
    auto result = nearest_facets({0, 0}, m.facets);
    REQUIRE(result.facets[0] == m.facets[0]);
    REQUIRE(result.pt == m.facets[0][0]);
    REQUIRE(result.distance == std::sqrt(2));
}

TEST_CASE("NearestFacetsOneFacetNotEndpoint", "[nearest_neighbors]") {
    Facet<2> f{{{-1,-1},{1,-1}}};
    Mesh<2> m{{f}};
    auto result = nearest_facets({0, 0}, m.facets);
    REQUIRE(result.facets[0] == m.facets[0]);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0,-1}, 2, 1e-14);
    REQUIRE(result.distance == 1.0);
}

TEST_CASE("NearestFacetsTwoFacetsNotEndpoint", "[nearest_neighbors]") {
    Facet<2> f{{{0.5,-1},{1,-1}}};
    Facet<2> f2{{{-1,-1},{0.5,-1}}};
    Mesh<2> m{{f,f2}};
    auto result = nearest_facets({0, 0}, m.facets);
    REQUIRE(result.facets[0] == m.facets[1]);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0,-1}, 2, 1e-14);
    REQUIRE(result.distance == 1.0);
}

TEST_CASE("NearestFacetsIntersection", "[nearest_neighbors]") {
    Facet<2> f{{{0,0},{1,0}}};
    Facet<2> f2{{{0,1},{0,0}}};
    Mesh<2> m{{f,f2}};
    auto result = nearest_facets({0, 0}, m.facets);
    REQUIRE(result.facets.size() == 2);
    bool correct_facets = result.facets[0] == f || result.facets[1] == f &&
        result.facets[0] == f2 || result.facets[1] == f2; 
    REQUIRE(correct_facets);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0,0}, 2, 1e-14);
    REQUIRE(result.distance == 0.0);
}

