#include "catch.hpp"
#include "closest_pt.h"
#include "vec_ops.h"
#include "numerics.h"

using namespace tbem;

TEST_CASE("NearestPointOneSegEndpoint", "[block_dof_map]") {
    Vec<Vec<double,2>,2> f{{{1,1},{2,1}}};
    auto result = ref_to_real(closest_pt_seg({0, 0}, f), f);
    REQUIRE(result == f[0]);
}

TEST_CASE("NearestPointOneSegNotEndpoint", "[block_dof_map]") {
    Vec<Vec<double,2>,2> f{{{-1,-1},{1,-1}}};
    auto result = ref_to_real(closest_pt_seg({0, 0}, f), f);
    REQUIRE(result == (Vec<double,2>{0,-1}));
}

TEST_CASE("NearestPointTriangle", "[block_dof_map]") {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_tri({0, 0, 1}, f), f);
    REQUIRE(result == (Vec<double,3>{0,0,0}));
}

TEST_CASE("NearestPointTriangleEasyVertex", "[block_dof_map]") {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_tri({0, 0, 1}, f), f);
    REQUIRE(result == (Vec<double,3>{0,0,0}));
}

TEST_CASE("NearestPointTriangleInterior", "[block_dof_map]") {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_tri({0.5, 0.1, 1}, f), f);
    REQUIRE(result == (Vec<double,3>{0.5, 0.1, 0}));
}

TEST_CASE("NearestPointTriangleOutside", "[block_dof_map]") {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_tri({1.5, 1.1, 1}, f), f);
    REQUIRE_CLOSE(result, Vec<double,3>{0.7, 0.3, 0}, 1e-12);
}
