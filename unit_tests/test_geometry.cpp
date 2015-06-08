#include "catch.hpp"
#include "geometry.h"
#include "util.h"

using namespace tbem;

TEST_CASE("bounding box", "[octree]") 
{
    std::vector<Ball<3>> es{
        {{1.0, 2.0, 0.0}, 0.0}, {{-1.0, 0.0, -3.0}, 0.0}, {{0.0, -2.0, 3.0}, 0.0}
    };
    auto bb = Box<3>::bounding_box(es);
    double hw[] = {1.0, 2.0, 3.0};
    double center[] = {0.0, 0.0, 0.0};
    REQUIRE_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    REQUIRE_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
}

TEST_CASE("degenerate bounding box", "[octree]") 
{
    std::vector<Ball<3>> es = {{{0.0, 0.0, 0.0}, 0.0}};
    auto bb = Box<3>::bounding_box(es);
    double hw[] = {0.0, 0.0, 0.0};
    double center[] = {0.0, 0.0, 0.0};
    REQUIRE_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    REQUIRE_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
}

TEST_CASE("bounding box no pts", "[octree]")
{
    auto bb = Box<3>::bounding_box({});
    REQUIRE(bb.center == (zeros<Vec<double,3>>::make()));
    REQUIRE(bb.half_width == (zeros<Vec<double,3>>::make()));
}

TEST_CASE("bounding box nonzero radius one pt", "[octree]")
{
    Ball<2> b{{1.0, 0.0}, 1.0};
    auto bb = Box<2>::bounding_box({b});
    REQUIRE(bb.center == (Vec<double,2>{1.0, 0.0}));
    REQUIRE(bb.half_width == (Vec<double,2>{1.0, 1.0}));
    REQUIRE(bb.in_box(b));
}

TEST_CASE("bounding box nonzero radius contains its pts", "[octree]")
{
    for (size_t i = 0; i < 100; i++) {
        auto balls = random_balls<2>(10, 0.05); 
        auto box = Box<2>::bounding_box(balls).expand_by_max_axis_multiple(1e-3);
        for (auto b: balls) {
            if (!box.in_box(b)) {
                for (auto b2: balls) {
                    std::cout << b2.center << " " << b2.radius << std::endl;
                }
                std::cout << box.center << " " << box.half_width << std::endl;
                std::cout << b.center << " " << b.radius << std::endl;
            }
            REQUIRE(box.in_box(b));
        }
    }
}

TEST_CASE("get subcell", "[octree]") 
{
    Box<3> box{{0, 1, 0}, {2, 2, 2}};
    auto child = box.get_subcell({1,0,1});
    REQUIRE_ARRAY_CLOSE(child.center, (Vec<double,3>{1,0,1}), 3, 1e-14);
    REQUIRE_ARRAY_CLOSE(child.half_width, (Vec<double,3>{1,1,1}), 3, 1e-14);
}

TEST_CASE("in box circle", "[octree]")
{
    Box<2> box{{0, 0}, {1, 1}};
    SECTION("inside") {
        REQUIRE(box.in_box(Ball<2>{{0, 0}, 0.5}));
    }
    SECTION("outside") {
        REQUIRE(!box.in_box(Ball<2>{{3, 0}, 0.5}));
    }
    SECTION("intersection") {
        REQUIRE(!box.in_box(Ball<2>{{1, 0}, 0.5}));
    }
}

TEST_CASE("in box", "[octree]") 
{
    Box<3> box{{0,1,0}, {2,2,2}};
    REQUIRE(box.in_box({0,2,0}, {false, false, false}));
    REQUIRE(!box.in_box({0,4,0}, {false, false, false}));
}

TEST_CASE("in box inclusive", "[octree]") 
{
    Box<3> box{{0,1,0}, {2,2,2}};
    REQUIRE(box.in_box({0,3,0}, {false, true, false}));
    REQUIRE(!box.in_box({0,3,0}, {false, false, false}));
}

TEST_CASE("containing subcell box", "[geometry]")
{
    Box<3> c{{0, 0, 0}, {1, 1, 1}};
    REQUIRE(c.find_containing_subcell({0.1, 0.1, 0.1}) == 7);
    REQUIRE(c.find_containing_subcell({0.1, -0.1, -0.1}) == 4);
    REQUIRE(c.find_containing_subcell({-0.1, -0.1, -0.1}) == 0);
}

TEST_CASE("Simple geometry", "[geometry]")
{
    Vec3<double> a = {{1.0, 1.0, 2.0}};
    Vec3<double> b = {{2.0, 0.5, -1.0}};

    SECTION("Norm") {
        auto c = normalized(b);
        normalize(b);
        double m = std::sqrt(5.25);
        Vec<double,3> exact{2 / m, 0.5 / m, -1 / m};
        REQUIRE_CLOSE(b, exact, 1e-6);
        REQUIRE_CLOSE(c, exact, 1e-6);
    }

    SECTION("Cross") {
        auto c = cross(a, b);
        Vec<double,3> exact{-2, 5, -1.5};
        REQUIRE_CLOSE(c, exact, 1e-6);
    }
}

TEST_CASE("WhichSidePT3D", "[geometry]") 
{
    auto val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,-1});
    REQUIRE(val == BEHIND);
    val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,1});
    REQUIRE(val == FRONT);
    val = which_side_point<3>({{{0,0,0}, {1,0,0}, {0,1,0}}}, {0,0,0});
    REQUIRE(val == INTERSECT);
}

TEST_CASE("WhichSidePT2D", "[geometry]") 
{
    auto val = which_side_point<2>({{{0,0}, {1,0}}}, {0,-1});
    REQUIRE(val == BEHIND);
    val = which_side_point<2>({{{0,0}, {1,0}}}, {0,1});
    REQUIRE(val == FRONT);
    val = which_side_point<2>({{{0,0}, {1,0}}}, {0,0});
    REQUIRE(val == INTERSECT);
}

TEST_CASE("SegmentSide", "[geometry]") 
{
    REQUIRE(facet_side<2>({FRONT, BEHIND}) == INTERSECT);
    REQUIRE(facet_side<2>({FRONT, FRONT}) == FRONT);
    REQUIRE(facet_side<2>({FRONT, INTERSECT}) == FRONT);
}

TEST_CASE("TriSide", "[geometry]") 
{
    REQUIRE(facet_side<3>({FRONT, FRONT, FRONT}) == FRONT);
    REQUIRE(facet_side<3>({INTERSECT, FRONT, FRONT}) == FRONT);
    REQUIRE(facet_side<3>({INTERSECT, INTERSECT, FRONT}) == FRONT);
    REQUIRE(facet_side<3>({BEHIND, INTERSECT, BEHIND}) == BEHIND);
}

TEST_CASE("OuterProductVectorVal", "[geometry]") 
{
    auto outer = outer_product<double>(Vec<double,2>{1.0, 1.0}, 0.5);
    Vec2<double> exact{0.5, 0.5};
    REQUIRE(outer == exact);
}

TEST_CASE("OuterProductVectorVal3D", "[geometry]") 
{
    auto outer = outer_product<double>(Vec<double,3>{1.0, 1.0, -2.0}, 0.5);
    Vec3<double> exact{0.5, 0.5, -1.0};
    REQUIRE(outer == exact);
}

TEST_CASE("OuterProductVectorVector", "[geometry]") 
{
    Vec2<double> K = {1.0, 1.0};
    Vec2<double> x = {3.0, 4.0};
    auto result = outer_product(K, x);
    Vec2<Vec2<double>> correct{{{3.0, 4.0}, {3.0, 4.0}}};
    REQUIRE(result == correct);
}

TEST_CASE("OuterProductTensorVector", "[geometry]") 
{
    Vec2<Vec2<double>> right{{{3.0, 0.0}, {0.0,4.0}}};
    Vec2<double> left = {1.0, 1.0};
    Vec2<Vec2<Vec2<double>>> correct{{
        {{{3.0, 0.0}, {0.0, 4.0}}},
        {{{3.0, 0.0}, {0.0, 4.0}}}
    }};
    auto result = outer_product(left, right);
    REQUIRE(result == correct);
}

TEST_CASE("OuterProductTensorVector3d", "[geometry]") 
{
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

TEST_CASE("InnerProductVecVec", "[geometry]") 
{
    Vec2<double> right{{3.0, 4.0}};
    Vec2<double> left = {1.0, 1.0};
    double correct = 7.0;
    auto result = dot_product(left, right);
    REQUIRE(result == correct);
}

TEST_CASE("Point in ball", "[geometry]") 
{
    Ball<3> b{{0,1,2}, 1};

    SECTION("in ball") {
        REQUIRE(b.inside({0, 0.5, 2.1}));
    }

    SECTION("not in ball") {
        REQUIRE(!b.inside({0, 0.0, 2.1}));
    }
}

TEST_CASE("balls from centers radii", "[geometry]")
{
    auto bs = balls_from_centers_radii<2>({{0, 2}, {1, 0}}, {3, 4});
    REQUIRE(bs.size() == 2);
    REQUIRE(bs[0].center == (Vec<double,2>{0, 2}));
    REQUIRE(bs[0].radius == 3);
}

TEST_CASE("Centroid 3D", "[geometry]") 
{
    Vec<Vec<double,3>,3> f{{
        {-1, 1, 0}, {4, -1, 0}, {0, 0, 9}
    }};
    REQUIRE_ARRAY_CLOSE(centroid(f), Vec<double,3>{1, 0, 3}, 3, 1e-15);
}

TEST_CASE("Bounding circle 2D", "[geometry]")
{
    Vec<Vec<double,2>,2> f{{
        {-1, 1}, {3, -1}
    }};
    REQUIRE_CLOSE(facet_ball(f).radius, std::sqrt(4 + 1), 1e-15);
}

TEST_CASE("Bounding circle 3D", "[geometry]")
{
    Vec<Vec<double,3>,3> f{{
        {-1, 1, 0}, {4, -1, 0}, {0, 0, 9}
    }};
    REQUIRE_CLOSE(facet_ball(f).radius, std::sqrt(37), 1e-15);
}

TEST_CASE("intersecting balls", "[geometry]") 
{
    auto bs = balls_from_centers_radii<2>({{0, 0}, {2, 0}}, {1, 1.1});
    REQUIRE(balls_intersect(bs[0], bs[1]));
}

TEST_CASE("not intersecting balls", "[geometry]") 
{
    auto bs = balls_from_centers_radii<2>({{0, 0}, {2, 0}}, {1, 0.5});
    REQUIRE(!balls_intersect(bs[0], bs[1]));
}
