#include "catch.hpp"
#include "gte_wrapper.h"

using namespace tbem;

TEST_CASE("seg seg no intersection", "[gte_wrapper]") 
{
    auto pts = seg_seg_intersection({{{0, 0}, {1, 0}}}, {{{-1, -1}, {-1, 1}}});
    REQUIRE(pts.size() == 0);
}

TEST_CASE("seg seg one intersection", "[gte_wrapper]") 
{
    auto pts = seg_seg_intersection({{{0, 0}, {1, 0}}}, {{{0.5, -1}, {0.5, 1}}});
    REQUIRE(pts.size() == 1);
    REQUIRE_ARRAY_CLOSE(pts[0], (Vec<double,2>{0.5, 0}), 2, 1e-12);
}

TEST_CASE("seg seg same segment", "[gte_wrapper]") 
{
    auto pts = seg_seg_intersection({{{0, 0}, {0, 1}}}, {{{0, 0}, {0, 1}}});
    REQUIRE(pts.size() == 2);
    REQUIRE_ARRAY_CLOSE(pts[0], (Vec<double,2>{0, 0}), 2, 1e-12);
    REQUIRE_ARRAY_CLOSE(pts[1], (Vec<double,2>{0, 1}), 2, 1e-12);
}

TEST_CASE("seg tri no intersection", "[gte_wrapper]") 
{
    auto pts = seg_tri_intersection(
        {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
        {{{-0.1, -0.1, -1}, {-0.1, -0.1, 1}}}
    );
    REQUIRE(pts.size() == 0);
}

TEST_CASE("seg tri intersection", "[gte_wrapper]") 
{
    auto pts = seg_tri_intersection(
        {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
        {{{0.1, 0.1, -1}, {0.1, 0.1, 1}}}
    );
    REQUIRE(pts.size() == 1);
    REQUIRE_ARRAY_CLOSE(pts[0], (Vec<double,3>{0.1, 0.1, 0}), 3, 1e-12);
}

TEST_CASE("in polygon", "[gte_wrapper]")
{
    std::vector<Vec<double,2>> poly{
        {0, 0}, {2, 0}, {2, 2}, {0, 2}
    };

    SECTION("yes") {
        in_polygon(poly, {1, 1});
    }

    SECTION("no") {
        in_polygon(poly, {-1, 0});
    }
}
