#include "catch.hpp"
#include "intersect_balls.h"
#include "geometry.h"
#include "octree.h"

using namespace tbem;

TEST_CASE("identical points", "[intersect_balls]") 
{
    auto bs = balls_from_centers_radii<3>({
        {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    }, {0.0, 0.0, 0.0});
    auto oct = make_octree(bs, 1);
    auto pts = intersect_balls(Ball<3>{{1.0, 2.0, 0.0}, 0.0}, bs, oct);
    REQUIRE(pts.size() == 1);
    REQUIRE(pts[0] == 0);
}

TEST_CASE("two identical points", "[intersect_balls]") 
{
    auto bs = balls_from_centers_radii<3>({
        {1.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    }, {0, 0, 0, 0});
    auto oct = make_octree(bs, 1);
    auto pts = intersect_balls(Ball<3>{{1.0, 2.0, 0.0}, 0.0}, bs, oct);
    REQUIRE(pts.size() == 2);
}

TEST_CASE("very close points", "[intersect_balls]")
{
    double eps = 1e-13;
    auto bs = balls_from_centers_radii<3>({
        {0.0, -2.0, 3.0}, {1.0, 2.0, 0.0}, {1.0, 2.0 - (eps / 10.), 0.0}
    }, {0, 0, 0});
    auto oct = make_octree(bs, 1);
    auto pts = intersect_balls(Ball<3>{{1.0, 2.0, 0.0}, eps}, bs, oct);
    REQUIRE(pts.size() == 2);
}

std::vector<Ball<2>> line_pts(size_t n) 
{
    std::vector<Ball<2>> pts;     
    for (size_t i = 0; i < n; i++) {
        pts.push_back({{static_cast<double>(i), 0}, 0});
        pts.push_back({{static_cast<double>(i + 1), 0}, 0});
    }
    return pts;
}

TEST_CASE("all pairs", "[intersect_balls]") 
{
    size_t n = 500;
    auto pts = line_pts(n);
    auto result = intersect_balls_all_pairs(pts, pts);
    // 1000 pairs like (k, k). 998 pairs like (k, k+1) and (k+1, k)
    REQUIRE(result.size() == 1998);
}
