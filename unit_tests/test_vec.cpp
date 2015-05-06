#include "catch.hpp"
#include "vec.h"
#include "numerics.h"
#include "util.h"

using namespace tbem;

struct Data {
    Vec3<double> a = {{1.0, 1.0, 2.0}};
    Vec3<double> b = {{2.0, 0.5, -1.0}};
};

TEST_CASE("VecOperations", "[vec]") 
{
    Vec3<double> a = {{1.0, 1.0, 2.0}};
    Vec3<double> b = {{2.0, 0.5, -1.0}};

    SECTION("Add") {
        auto c = a + b; 
        b += a;
        Vec<double,3> exact{3.0, 1.5, 1.0};
        REQUIRE_CLOSE(b, exact, 1e-12);
        REQUIRE_CLOSE(c, exact, 1e-12);
    }

    SECTION("Sub") {
        auto c = a - b; 
        Vec<double,3> exact{-1.0, 0.5, 3.0};
        REQUIRE_CLOSE(c, exact, 1e-12);
    }
    
    SECTION("Mul") {
        auto c = a * b; 
        Vec<double,3> exact{2.0, 0.5, -2.0};
        REQUIRE_CLOSE(c, exact, 1e-12);
    }

    SECTION("Div") {
        auto c = a / b; 
        Vec<double,3> exact{0.5, 2.0, -2.0};
        REQUIRE_CLOSE(c, exact, 1e-12);
    }

    SECTION("Norm") {
        auto c = normalized(b);
        normalize(b);
        double m = std::sqrt(5.25);
        Vec<double,3> exact{2 / m, 0.5 / m, -1 / m};
        REQUIRE_CLOSE(b, exact, 1e-6);
        REQUIRE_CLOSE(c, exact, 1e-6);
    }

    SECTION("Negate") {
        auto c = -b;
        Vec<double,3> exact{-2.0, -0.5, 1.0};
        REQUIRE_CLOSE(c, exact, 1e-6);
    }

    SECTION("Cross") {
        auto c = cross(a, b);
        Vec<double,3> exact{-2, 5, -1.5};
        REQUIRE_CLOSE(c, exact, 1e-6);
    }
}




TEST_CASE("VecPrint", "[vec]") {
    Vec3<double> a = {1.0, 2.0, 3.0};
    std::stringstream output_buf;
    output_buf << a;
    REQUIRE(output_buf.str() == "(1, 2, 3)");
}

TEST_CASE("WhichSidePT3D", "[vec]") {
    auto val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,-1});
    REQUIRE(val == BEHIND);
    val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,1});
    REQUIRE(val == FRONT);
    val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,0});
    REQUIRE(val == INTERSECT);
}

TEST_CASE("WhichSidePT2D", "[vec]") {
    auto val = which_side_point<2>({{{0,0}, {1,0}}}, {0,-1});
    REQUIRE(val == BEHIND);
    val = which_side_point<2>({{{0,0}, {1,0}}}, {0,1});
    REQUIRE(val == FRONT);
    val = which_side_point<2>({{{0,0}, {1,0}}}, {0,0});
    REQUIRE(val == INTERSECT);
}

TEST_CASE("SegmentSide", "[vec]") {
    REQUIRE(facet_side<2>({FRONT, BEHIND}) == INTERSECT);
    REQUIRE(facet_side<2>({FRONT, FRONT}) == FRONT);
    REQUIRE(facet_side<2>({FRONT, INTERSECT}) == FRONT);
}

TEST_CASE("TriSide", "[vec]") {
    REQUIRE(facet_side<3>({FRONT, FRONT, FRONT}) == FRONT);
    REQUIRE(facet_side<3>({INTERSECT, FRONT, FRONT}) == FRONT);
    REQUIRE(facet_side<3>({INTERSECT, INTERSECT, FRONT}) == FRONT);
    REQUIRE(facet_side<3>({BEHIND, INTERSECT, BEHIND}) == BEHIND);
}

TEST_CASE("OuterProductVectorVal", "[vec]") {
    auto outer = outer_product<double>(Vec<double,2>{1.0, 1.0}, 0.5);
    Vec2<double> exact{0.5, 0.5};
    REQUIRE(outer == exact);
}

TEST_CASE("OuterProductVectorVal3D", "[vec]") {
    auto outer = outer_product<double>(Vec<double,3>{1.0, 1.0, -2.0}, 0.5);
    Vec3<double> exact{0.5, 0.5, -1.0};
    REQUIRE(outer == exact);
}

TEST_CASE("OuterProductVectorVector", "[vec]") {
    Vec2<double> K = {1.0, 1.0};
    Vec2<double> x = {3.0, 4.0};
    auto result = outer_product(K, x);
    Vec2<Vec2<double>> correct{{{3.0, 4.0}, {3.0, 4.0}}};
    REQUIRE(result == correct);
}

TEST_CASE("OuterProductTensorVector", "[vec]") {
    Vec2<Vec2<double>> right{{{3.0, 0.0}, {0.0,4.0}}};
    Vec2<double> left = {1.0, 1.0};
    Vec2<Vec2<Vec2<double>>> correct{{
        {{{3.0, 0.0}, {0.0, 4.0}}},
        {{{3.0, 0.0}, {0.0, 4.0}}}
    }};
    auto result = outer_product(left, right);
    REQUIRE(result == correct);
}

TEST_CASE("OuterProductTensorVector3d", "[vec]") {
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
    REQUIRE(result == correct);
}

TEST_CASE("InnerProductVecVec", "[vec]") {
    Vec2<double> right{{3.0, 4.0}};
    Vec2<double> left = {1.0, 1.0};
    double correct = 7.0;
    auto result = dot_product(left, right);
    REQUIRE(result == correct);
}

TEST_CASE("ZerosTensor", "[vec]") {
    auto z = zeros<Vec2<Vec2<double>>>::make();
    Vec2<Vec2<double>> c{{{0.0, 0.0}, {0.0, 0.0}}};
    REQUIRE(z == c);
}

TEST_CASE("OnesTensor", "[vec]") {
    auto z = ones<Vec2<Vec2<double>>>::make();
    Vec2<Vec2<double>> c{{{1.0, 1.0}, {1.0, 1.0}}};
    REQUIRE(z == c);
}

TEST_CASE("ConstantTensor", "[vec]") 
{
    auto z = constant<Vec2<Vec2<double>>>::make(2.2);
    Vec2<Vec2<double>> c{{{2.2, 2.2}, {2.2, 2.2}}};
    REQUIRE(z == c);
}

TEST_CASE("Declaration ordering shouldn't matter", "[vec]")
{
    Vec2<Vec1<Vec1<double>>> x;
    Vec2<Vec1<Vec1<double>>> y;
    x[0][0][0] = 5.0;
    x[1][0][0] = 10.0;
    y[0][0][0] = 15.0;
    y[1][0][0] = 10.0;
    x += y;
    REQUIRE(x[0][0][0] == 20.0);
    REQUIRE(x[1][0][0] == 20.0);
}

TEST_CASE("Inequality", "[vec]") {
    Vec2<double> a{1, 0};
    auto res = a < 0.5;
    REQUIRE(res == (Vec2<bool>{false, true}));
}
