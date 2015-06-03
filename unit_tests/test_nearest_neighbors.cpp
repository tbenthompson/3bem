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
    auto oct = make_octree(es, 1);
    auto pts = nearby_points({1.0, 2.0, 0.0}, es, oct, 0.0);
    REQUIRE(pts.size() == 1);
    REQUIRE(pts[0] == 0);
}

TEST_CASE("TwoIdenticalPoints", "[nearest_neighbors]") 
{
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    };
    auto oct = make_octree(es, 1);
    auto pts = nearby_points({1.0, 2.0, 0.0}, es, oct, 0.0);
    REQUIRE(pts.size() == 2);
}

TEST_CASE("ClosePoints", "[nearest_neighbors]")
{
    double eps = 1e-13;
    std::vector<Vec<double,3>> es{
        {0.0, -2.0, 3.0}, {1.0, 2.0, 0.0}, {1.0, 2.0 - (eps / 10.), 0.0}
    };
    auto oct = make_octree(es, 1);
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

TEST_CASE("AllPairs", "[nearest_neighbors]") 
{
    size_t n = 500;
    auto pts = line_pts(n);
    auto oct = make_octree(pts, 50);
    auto result = nearby_points_all_pairs(pts, pts, oct, oct, 0.0);
    // 1000 pairs like (k, k). 998 pairs like (k, k+1) and (k+1, k)
    REQUIRE(result.size() == 1998);
}

TEST_CASE("AllPairsPerformance", "[nearest_neighbors]") 
{
    size_t n = 50000;
    auto pts = line_pts(n);
    auto oct = make_octree(pts, 50);
    //TODO: Capacity test
    auto result = nearby_points_all_pairs(pts, pts, oct, oct, 0.0);
}

TEST_CASE("test nearest facet", "[nearest_neighbors]") 
{
    //Test the fast version against the brute force version through many
    //random examples
    size_t n_facets = 10;
    for (size_t i = 0; i < 100; i++) {
        std::vector<Facet<2>> fs;
        for (size_t j = 0; j < n_facets; j++) {
            auto vs = random_pts<2>(2);
            fs.push_back({vs[0], vs[1]});
        }
        // auto brute_force = nearest_facet_brute_force(
    }
}

TEST_CASE("nearest facets on line", "[nearest_neighbors]") 
{
    Facet<2> f{{{1,1},{2,1}}};
    Mesh<2> m{{f}};
    auto result = nearest_facets({0, 0}, m.facets);
    REQUIRE(result.facets.size() == 0);
    REQUIRE(result.nearest_facet == f);
    REQUIRE(result.pt == m.facets[0][0]);
    REQUIRE(result.distance == std::sqrt(2));
}

TEST_CASE("nearest facets out of line", "[nearest_neighbors]") 
{
    Facet<2> f{{{-1,-1},{1,-1}}};
    Mesh<2> m{{f}};
    auto result = nearest_facets({0, 0}, m.facets);
    REQUIRE(result.nearest_facet == m.facets[0]);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0,-1}, 2, 1e-14);
    REQUIRE(result.distance == 1.0);
}

TEST_CASE("nearest facets two edges", "[nearest_neighbors]") 
{
    Facet<2> f{{{0.5,-1},{1,-1}}};
    Facet<2> f2{{{-1,-1},{0.5,-1}}};
    Mesh<2> m{{f,f2}};
    auto result = nearest_facets({0, 0}, m.facets);
    REQUIRE(result.nearest_facet == f2);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0,-1}, 2, 1e-14);
    REQUIRE(result.distance == 1.0);
}

TEST_CASE("nearest facets at intersection", "[nearest_neighbors]") 
{
    Facet<2> f{{{0,0},{1,0}}};
    Facet<2> f2{{{0,1},{0,0}}};
    Mesh<2> m{{f,f2}};
    auto result = nearest_facets({0, 0}, m.facets);
    REQUIRE(result.facets.size() == 1);
    bool correct_facets = (result.facets[0] == f || result.nearest_facet == f) &&
        (result.facets[0] == f2 || result.nearest_facet == f2);
    REQUIRE(correct_facets);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0,0}, 2, 1e-14);
    REQUIRE(result.distance == 0.0);
}

