#include "catch.hpp"
#include "intersect_balls.h"

using namespace tbem;

TEST_CASE("identical points", "[intersect_balls]") 
{
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    };
    auto oct = make_octree(es, 1);
    auto pts = intersect_balls({1.0, 2.0, 0.0}, es, oct, 0.0);
    REQUIRE(pts.size() == 1);
    REQUIRE(pts[0] == 0);
}

TEST_CASE("two identical points", "[intersect_balls]") 
{
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    };
    auto oct = make_octree(es, 1);
    auto pts = intersect_balls({1.0, 2.0, 0.0}, es, oct, 0.0);
    REQUIRE(pts.size() == 2);
}

TEST_CASE("very close points", "[intersect_balls]")
{
    double eps = 1e-13;
    std::vector<Vec<double,3>> es{
        {0.0, -2.0, 3.0}, {1.0, 2.0, 0.0}, {1.0, 2.0 - (eps / 10.), 0.0}
    };
    auto oct = make_octree(es, 1);
    auto pts = intersect_balls({1.0, 2.0, 0.0}, es, oct, eps);
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

TEST_CASE("all pairs", "[intersect_balls]") 
{
    size_t n = 500;
    auto pts = line_pts(n);
    auto oct = make_octree(pts, 50);
    auto result = intersect_balls_all_pairs(pts, pts, oct, oct, 0.0);
    // 1000 pairs like (k, k). 998 pairs like (k, k+1) and (k+1, k)
    REQUIRE(result.size() == 1998);
}

TEST_CASE("all pairs performance", "[intersect_balls]") 
{
    size_t n = 50000;
    auto pts = line_pts(n);
    auto oct = make_octree(pts, 50);
    //TODO: Capacity test
    auto result = intersect_balls_all_pairs(pts, pts, oct, oct, 0.0);
}

