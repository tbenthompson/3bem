#include "UnitTest++.h"
#include "closest_pt.h"
#include "vec_ops.h"
#include "numerics.h"

using namespace tbem;

TEST(NearestPointOneSegEndpoint) {
    Vec<Vec<double,2>,2> f{{{1,1},{2,1}}};
    auto result = ref_to_real(closest_pt_seg({0, 0}, f), f);
    CHECK_EQUAL(result, f[0]);
}

TEST(NearestPointOneSegNotEndpoint) {
    Vec<Vec<double,2>,2> f{{{-1,-1},{1,-1}}};
    auto result = ref_to_real(closest_pt_seg({0, 0}, f), f);
    CHECK_EQUAL(result, (Vec<double,2>{0,-1}));
}

TEST(NearestPointTriangle) {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_tri({0, 0, 1}, f), f);
    CHECK_EQUAL(result, (Vec<double,3>{0,0,0}));
}

TEST(NearestPointTriangleEasyVertex) {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_tri({0, 0, 1}, f), f);
    CHECK_EQUAL(result, (Vec<double,3>{0,0,0}));
}

TEST(NearestPointTriangleInterior) {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_tri({0.5, 0.1, 1}, f), f);
    CHECK_EQUAL(result, (Vec<double,3>{0.5, 0.1, 0}));
}

TEST(NearestPointTriangleOutside) {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real(closest_pt_tri({1.5, 1.1, 1}, f), f);
    CHECK_CLOSE(result, (Vec<double,3>{0.7, 0.3, 0}), 1e-12);
}

int main() {
    return UnitTest::RunAllTests();
}
