#include "UnitTest++.h"
#include "vec.h"
#include "numerics.h"
#include "util.h"

using namespace tbem;

struct Data {
    Vec3<double> a = {{1.0, 1.0, 2.0}};
    Vec3<double> b = {{2.0, 0.5, -1.0}};
};

TEST_FIXTURE(Data, VecAdd) {
    auto c = a + b; 
    b += a;
    double exact[3] = {3.0, 1.5, 1.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-12);
    CHECK_ARRAY_CLOSE(b, exact, 3, 1e-12);
}

TEST_FIXTURE(Data, VecSub) {
    auto c = a - b; 
    double exact[3] = {-1.0, 0.5, 3.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-12);
}

TEST_FIXTURE(Data, VecMul) {
    auto c = a * b; 
    double exact[3] = {2.0, 0.5, -2.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-12);
}

TEST_FIXTURE(Data, VecDiv) {
    auto c = a / b; 
    double exact[3] = {0.5, 2.0, -2.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-12);
}

TEST_FIXTURE(Data, VecNorm) {
    auto c = normalized(b);
    normalize(b);
    double m = std::sqrt(5.25);
    double exact[3] = {2 / m, 0.5 / m, -1 / m};
    CHECK_ARRAY_CLOSE(b, exact, 3, 1e-6);
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-6);
}

TEST_FIXTURE(Data, VecNegate) {
    auto c = -b;
    double exact[3] = {-2.0, -0.5, 1.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-6);
}

TEST_FIXTURE(Data, VecCross) {
    auto c = cross(a, b);
    double exact[3] = {-2, 5, -1.5};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-6);
}

TEST(VecPrint) {
    Vec3<double> a = {1.0, 2.0, 3.0};
    std::stringstream output_buf;
    output_buf << a;
    CHECK_EQUAL(output_buf.str(), "(1, 2, 3)");
}

TEST(WhichSidePT3D) {
    auto val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,-1});
    CHECK_EQUAL(val, BEHIND);
    val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,1});
    CHECK_EQUAL(val, FRONT);
    val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,0});
    CHECK_EQUAL(val, INTERSECT);
}

TEST(WhichSidePT2D) {
    auto val = which_side_point<2>({{{0,0}, {1,0}}}, {0,-1});
    CHECK_EQUAL(val, BEHIND);
    val = which_side_point<2>({{{0,0}, {1,0}}}, {0,1});
    CHECK_EQUAL(val, FRONT);
    val = which_side_point<2>({{{0,0}, {1,0}}}, {0,0});
    CHECK_EQUAL(val, INTERSECT);
}

TEST(SegmentSide) {
    CHECK_EQUAL(facet_side<2>({FRONT, BEHIND}), INTERSECT);
    CHECK_EQUAL(facet_side<2>({FRONT, FRONT}), FRONT);
    CHECK_EQUAL(facet_side<2>({FRONT, INTERSECT}), FRONT);
}

TEST(TriSide) {
    CHECK_EQUAL(facet_side<3>({FRONT, FRONT, FRONT}), FRONT);
    CHECK_EQUAL(facet_side<3>({INTERSECT, FRONT, FRONT}), FRONT);
    CHECK_EQUAL(facet_side<3>({INTERSECT, INTERSECT, FRONT}), FRONT);
    CHECK_EQUAL(facet_side<3>({BEHIND, INTERSECT, BEHIND}), BEHIND);
}

TEST(OuterProductVectorVal) {
    CHECK_EQUAL((outer_product<double>(Vec<double,2>{1.0, 1.0}, 0.5)), 
                (Vec2<double>{0.5, 0.5}));
}

TEST(OuterProductVectorVal3D) {
    CHECK_EQUAL((outer_product<double>(Vec<double,3>{1.0, 1.0, -2.0}, 0.5)), 
                (Vec3<double>{0.5, 0.5, -1.0}));
}

TEST(OuterProductVectorVector) {
    Vec2<double> K = {1.0, 1.0};
    Vec2<double> x = {3.0, 4.0};
    auto result = outer_product(K, x);
    Vec2<Vec2<double>> correct{{{3.0, 4.0}, {3.0, 4.0}}};
    CHECK_EQUAL(result, correct);
}

TEST(OuterProductTensorVector) {
    Vec2<Vec2<double>> right{{{3.0, 0.0}, {0.0,4.0}}};
    Vec2<double> left = {1.0, 1.0};
    Vec2<Vec2<Vec2<double>>> correct{{
        {{{3.0, 0.0}, {0.0, 4.0}}},
        {{{3.0, 0.0}, {0.0, 4.0}}}
    }};
    auto result = outer_product(left, right);
    CHECK_EQUAL(result, correct);
}

TEST(OuterProductTensorVector3d) {
    Vec3<Vec3<double>> right{{
        {3.0, 0.0, 1.0}, {0.0, 4.0, -2.0}, {0.5, 7.0, -3.0}
    }};
    Vec3<double> left = {1.0, -1.0, 1.0};
    Vec3<Vec3<Vec3<double>>> correct{{
        {{{3.0, 0.0, 1.0}, {0.0, 4.0, -2.0}, {0.5, 7.0, -3.0}}},
        {{{-3.0, 0.0, -1.0}, {0.0, -4.0, 2.0}, {-0.5, -7.0, 3.0}}},
        {{{3.0, 0.0, 1.0}, {0.0, 4.0, -2.0}, {0.5, 7.0, -3.0}}}
    }};
    auto result = outer_product(left, right);
    CHECK_EQUAL(result, correct);
}

TEST(InnerProductVecVec) {
    Vec2<double> right{{3.0, 4.0}};
    Vec2<double> left = {1.0, 1.0};
    double correct = 7.0;
    auto result = dot_product(left, right);
    CHECK_EQUAL(result, correct);
}

TEST(ZerosTensor) {
    auto z = zeros<Vec2<Vec2<double>>>::make();
    Vec2<Vec2<double>> c{{{0.0, 0.0}, {0.0, 0.0}}};
    CHECK_EQUAL(z, c);
}

TEST(OnesTensor) {
    auto z = ones<Vec2<Vec2<double>>>::make();
    Vec2<Vec2<double>> c{{{1.0, 1.0}, {1.0, 1.0}}};
    CHECK_EQUAL(z, c);
}

TEST(ConstantTensor) {
    auto z = constant<Vec2<Vec2<double>>>::make(2.2);
    Vec2<Vec2<double>> c{{{2.2, 2.2}, {2.2, 2.2}}};
    CHECK_EQUAL(z, c);
}

int main() {
    return UnitTest::RunAllTests();
}
