#include "catch.hpp"
#include "gte_wrapper.h"
#include "numerics.h"

using namespace tbem;

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

TEST_CASE("nearest point endpoint", "[gte_wrapper]") 
{
    Vec<Vec<double,2>,2> f{{{1,1},{2,1}}};
    auto result = ref_to_real(closest_pt_facet<2>({0, 0}, f), f);
    REQUIRE(result == f[0]);
}

TEST_CASE("nearest point not endpoint", "[gte_wrapper]") 
{
    Vec<Vec<double,2>,2> f{{{-1,-1},{1,-1}}};
    auto result = ref_to_real(closest_pt_facet<2>({0, 0}, f), f);
    REQUIRE(result == (Vec<double,2>{0,-1}));
}

TEST_CASE("nearest point triangle above vertex", "[gte_wrapper]") 
{
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_facet<3>({0, 0, 1}, f), f);
    REQUIRE(result == (Vec<double,3>{0,0,0}));
}

TEST_CASE("nearest point triangle vertex", "[gte_wrapper]") 
{
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_facet<3>({0, 0, 1}, f), f);
    REQUIRE(result == (Vec<double,3>{0,0,0}));
}

TEST_CASE("nearest point triangle interior", "[gte_wrapper]") 
{
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_facet<3>({0.5, 0.1, 1}, f), f);
    REQUIRE(result == (Vec<double,3>{0.5, 0.1, 0}));
}

TEST_CASE("nearest point triangle outside", "[gte_wrapper]") 
{
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_facet<3>({1.5, 1.1, 1}, f), f);
    REQUIRE_CLOSE(result, Vec<double,3>{0.7, 0.3, 0}, 1e-12);
}
