#include "catch.hpp"
#include "vec_ops.h"

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

    SECTION("Negate") {
        auto c = -b;
        Vec<double,3> exact{-2.0, -0.5, 1.0};
        REQUIRE_CLOSE(c, exact, 1e-6);
    }
}

TEST_CASE("VecPrint", "[vec]") {
    Vec3<double> a = {1.0, 2.0, 3.0};
    std::stringstream output_buf;
    output_buf << a;
    REQUIRE(output_buf.str() == "(1, 2, 3)");
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
