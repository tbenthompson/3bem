#include "catch.hpp"
#include "gte_wrapper.h"
#include "numerics.h"
#include "geometry.h"

using namespace tbem;

TEST_CASE("nearest point endpoint", "[gte_wrapper]") 
{
    Vec<Vec<double,2>,2> f{{{1,1},{2,1}}};
    auto result = closest_pt_facet<2>({0, 0}, f);
    REQUIRE(result.pt == f[0]);
}

TEST_CASE("nearest point not endpoint", "[gte_wrapper]") 
{
    Vec<Vec<double,2>,2> f{{{-1,-1},{1,-1}}};
    auto result = closest_pt_facet<2>({0, 0}, f);
    REQUIRE(result.pt == (Vec<double,2>{0,-1}));
}

TEST_CASE("nearest point triangle above vertex", "[gte_wrapper]") 
{
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = closest_pt_facet<3>({0, 0, 1}, f);
    REQUIRE(result.pt == (Vec<double,3>{0,0,0}));
}

TEST_CASE("nearest point triangle vertex", "[gte_wrapper]") 
{
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = closest_pt_facet<3>({0, 0, 1}, f);
    REQUIRE(result.pt == (Vec<double,3>{0,0,0}));
}

TEST_CASE("nearest point triangle interior", "[gte_wrapper]") 
{
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = closest_pt_facet<3>({0.5, 0.1, 1}, f);
    REQUIRE(result.pt == (Vec<double,3>{0.5, 0.1, 0}));
}

TEST_CASE("nearest point triangle outside", "[gte_wrapper]") 
{
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = closest_pt_facet<3>({1.5, 1.1, 1}, f);
    REQUIRE_CLOSE(result.pt, Vec<double,3>{0.7, 0.3, 0}, 1e-12);
}

TEST_CASE("seg seg no intersection", "[gte_wrapper]") 
{
    auto pts = seg_facet_intersection<2>({{{0, 0}, {1, 0}}}, {{{-1, -1}, {-1, 1}}});
    REQUIRE(pts.size() == 0);
}

TEST_CASE("seg seg one intersection", "[gte_wrapper]") 
{
    auto pts = seg_facet_intersection<2>({{{0, 0}, {1, 0}}}, {{{0.5, -1}, {0.5, 1}}});
    REQUIRE(pts.size() == 1);
    REQUIRE_ARRAY_CLOSE(pts[0], (Vec<double,2>{0.5, 0}), 2, 1e-12);
}

TEST_CASE("seg seg same segment", "[gte_wrapper]") 
{
    auto pts = seg_facet_intersection<2>({{{0, 0}, {0, 1}}}, {{{0, 0}, {0, 1}}});
    REQUIRE(pts.size() == 2);
    REQUIRE_ARRAY_CLOSE(pts[0], (Vec<double,2>{0, 0}), 2, 1e-12);
    REQUIRE_ARRAY_CLOSE(pts[1], (Vec<double,2>{0, 1}), 2, 1e-12);
}

TEST_CASE("seg tri no intersection", "[gte_wrapper]") 
{
    auto pts = seg_facet_intersection<3>(
        {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
        {{{-0.1, -0.1, -1}, {-0.1, -0.1, 1}}}
    );
    REQUIRE(pts.size() == 0);
}

TEST_CASE("seg tri intersection", "[gte_wrapper]") 
{
    auto pts = seg_facet_intersection<3>(
        {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
        {{{0.1, 0.1, -1}, {0.1, 0.1, 1}}}
    );
    REQUIRE(pts.size() == 1);
    REQUIRE_ARRAY_CLOSE(pts[0], (Vec<double,3>{0.1, 0.1, 0}), 3, 1e-12);
}

TEST_CASE("seg tri edge overlap", "[gte_wrapper]")
{
    auto pts = seg_facet_intersection<3>(
        {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
        {{{0, 0, 0}, {1, 0, 0}}}
    );
    REQUIRE(pts.size() == 2);
    REQUIRE(pts[0] == (Vec<double,3>{0, 0, 0}));
    REQUIRE(pts[1] == (Vec<double,3>{1, 0, 0}));
}

TEST_CASE("seg tri edge partial overlap", "[gte_wrapper]")
{
    auto pts = seg_facet_intersection<3>(
        {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
        {{{0, 0, 0}, {0.5, 0, 0}}}
    );
    REQUIRE(pts.size() == 2);
    REQUIRE(pts[0] == (Vec<double,3>{0, 0, 0}));
    REQUIRE(pts[1] == (Vec<double,3>{0.5, 0, 0}));
}

TEST_CASE("seg tri edge partial overlap 2", "[gte_wrapper]")
{
    auto pts = seg_facet_intersection<3>(
        {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}},
        {{{0.4, 0.2, 0}, {0.5, 0.1, 0}}}
    );
    REQUIRE(pts.size() == 2);
    REQUIRE(pts[0] == (Vec<double,3>{0.4, 0.2, 0}));
    REQUIRE(pts[1] == (Vec<double,3>{0.5, 0.1, 0}));
}

TEST_CASE("is point in polygon", "[gte_wrapper]")
{
    std::vector<Vec<double,2>> poly{
        {0, 0}, {2, 0}, {2, 2}, {0, 2}
    };

    SECTION("yes") {
        is_point_in_polygon({1, 1}, poly);
    }

    SECTION("no") {
        is_point_in_polygon({-1, 0}, poly);
    }
}

TEST_CASE("intersect ball box 3d", "[gte_wrapper]") 
{
    Box<3> b{{0, 0, 0}, {1, 1, 1}};

    SECTION("fully contained is intersection") {
        REQUIRE(is_intersection_box_ball<3>(b, {{0.5, 0.5, 0.5}, 0.1}));
    }

    SECTION("partially contained is intersection") {
        REQUIRE(is_intersection_box_ball<3>(b, {{1.1, 0.5, 0.5}, 0.9}));
    }

    SECTION("outside") {
        REQUIRE(!is_intersection_box_ball<3>(b, {{1.5, 0.5, 0.5}, 0.1}));
    }
}

TEST_CASE("intersect ball box 2d", "[gte_wrapper]") 
{
    Box<3> b{{3, 3}, {1, 1}};

    SECTION("fully contained is intersection") {
        REQUIRE(is_intersection_box_ball<3>(b, {{3.5, 3.5}, 0.1}));
        REQUIRE(is_intersection_box_ball<3>(b, {{3.0, 3.0}, 0.0}));
    }

    SECTION("partially contained is intersection") {
        REQUIRE(is_intersection_box_ball<3>(b, {{4.1, 3.5}, 0.9}));
        REQUIRE(is_intersection_box_ball<3>(b, {{3.9, 3.9}, 0.2}));
    }

    SECTION("outside") {
        REQUIRE(!is_intersection_box_ball<3>(b, {{4.5, 3.5}, 0.1}));
    }
}

TEST_CASE("is point in ball", "[gte_wrapper]")
{
    REQUIRE(is_point_in_ball<2>({2.01, 1.05}, {{2, 1}, 0.1}));
}

TEST_CASE("is point in triangle", "[gte_wrapper]")
{
    REQUIRE(is_point_in_triangle({0.1, 1.0}, {{{0.05, 0.0}, {1.0, 1.0}, {0.0, 1.5}}}));
}

TEST_CASE("tri tri parallel", "[gte_wrapper]")
{
    Vec<Vec<double,3>,3> fA{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> fB{{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}}};
    auto result = facet_facet_intersection(fA, fB);
    REQUIRE(result.size() == 0);
}

TEST_CASE("tri tri no intersection", "[gte_wrapper]")
{
    Vec<Vec<double,3>,3> fA{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> fB{{{0, 0, 1}, {1, 0, 1}, {0, 1, 0.1}}};
    REQUIRE(facet_facet_intersection(fA, fB).size() == 0);
}

TEST_CASE("tri tri point intersection", "[gte_wrapper]")
{
    Vec<Vec<double,3>,3> fA{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> fB{{{0, 0, 1}, {1, 0, 1}, {0, 1, 0}}};
    REQUIRE(facet_facet_intersection(fA, fB).size() == 1);
}

TEST_CASE("tri tri point plane intersection not inside triangle", "[gte_wrapper]")
{
    Vec<Vec<double,3>,3> fA{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> fB{{{0, 0, 1}, {0, 0, 2}, {0, 1, 1}}};
    REQUIRE(facet_facet_intersection(fA, fB).size() == 0);
}

TEST_CASE("tri tri line intersection", "[gte_wrapper]")
{
    Vec<Vec<double,3>,3> fA{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> fB{{{0, 0, 1}, {1, 0, 1}, {0, 1, -0.5}}};
    auto result = facet_facet_intersection(fA, fB);
    REQUIRE(result.size() == 2);
}

TEST_CASE("tri tri same tri", "[gte_wrapper]")
{
    Vec<Vec<double,3>,3> fA{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> fB{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    REQUIRE(facet_facet_intersection(fA, fB).size() == 3);
}

TEST_CASE("tri tri coplanar and overlapping", "[gte_wrapper]")
{
    Vec<Vec<double,3>,3> fA{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> fB{{{0.5, 0, 0}, {1.5, 0, 0}, {0.5, 1, 0}}};
    REQUIRE(facet_facet_intersection(fA, fB).size() == 3);
}

TEST_CASE("tri tri only planes intersect", "[gte_wrapper]")
{
    Vec<Vec<double,3>,3> fA{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    Vec<Vec<double,3>,3> fB{{{1.5, 0, 0}, {2.5, 0, 0}, {1.5, 1, 0}}};
    REQUIRE(facet_facet_intersection(fA, fB).size() == 0);
}
