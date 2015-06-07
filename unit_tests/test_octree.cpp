#include "catch.hpp"
#include "octree.h"
#include "util.h"
#include "vec_ops.h"
#include "geometry.h"
#include <set>
using namespace tbem;

template <size_t dim>
std::vector<Ball<dim>> random_balls(size_t n, double r_max)
{
    auto pts = random_pts<dim>(n);
    auto r = random_list(n, 0, r_max);
    std::vector<Ball<dim>> balls(n);
    for (size_t i = 0; i < n; i++) {
        balls[i] = {pts[i], r[i]};
    }
    return balls;
}

std::vector<Ball<3>> three_pts() {
    return {
        {{1.0, 2.0, 0.0}, 0.0}, {{-1.0, 0.0, -3.0}, 0.0}, {{0.0, -2.0, 3.0}, 0.0}
    };
}
TEST_CASE("bounding box", "[octree]") 
{
    auto es = three_pts();
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

TEST_CASE("make child idx 3d", "[octree]") 
{
    REQUIRE(make_child_idx<3>(0) == (Vec3<size_t>{0, 0, 0}));
    REQUIRE(make_child_idx<3>(2) == (Vec3<size_t>{0, 1, 0}));
    REQUIRE(make_child_idx<3>(7) == (Vec3<size_t>{1, 1, 1}));
}

TEST_CASE("make child idx 2d", "[octree]") 
{
    REQUIRE(make_child_idx<2>(1) == (Vec2<size_t>{0, 1}));
    REQUIRE(make_child_idx<2>(3) == (Vec2<size_t>{1, 1}));
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

TEST_CASE("check law of large numbers", "[octree]") 
{
    int n = 1000000;
    auto pts = random_pts<3>(n);
    //TODO: Make a octree capacity test
    auto tree = make_octree(pts, 100);
    for (size_t i = 0; i < 8; i++) {
        int n_pts = tree.children[i]->data.indices.size();
        int diff = abs(n_pts - (n / 8));
        CHECK(diff < (n / 32));
    }
}

TEST_CASE("one level octree", "[octree]") 
{
    auto es = three_pts();
    auto oct = make_octree(es, 4);
    REQUIRE(oct.data.level == 0);
    for (size_t i = 0; i < 8; i++) {
        REQUIRE(oct.children[i] == nullptr);
    }
    REQUIRE(oct.data.indices.size() == 3);
}

TEST_CASE("many level octree", "[octree]") 
{
    auto pts = random_pts<3>(1000);
    auto oct = make_octree(pts, 4);
    REQUIRE(oct.data.indices.size() == 1000);
    REQUIRE(oct.children[0]->data.level == 1);
}

TEST_CASE("degenerate line octree in 2d", "[octree]") 
{
    std::vector<Vec<double,2>> pts;
    for (size_t i = 0; i < 100; i++) {
        pts.push_back({static_cast<double>(i), 0.0});
    }
    auto oct = make_octree(pts, 1);
    REQUIRE(!oct.is_leaf());
}

template <size_t dim>
size_t n_pts(const Octree<dim>& cell) 
{
    if (cell.is_leaf()) {
        return cell.data.indices.size();
    }
    size_t n_child_pts = 0;
    for (size_t c = 0; c < Octree<dim>::split; c++) {
        if (cell.children[c] == nullptr) {
            continue;
        }
        n_child_pts += n_pts(*cell.children[c]);
    }
    return n_child_pts;
}

template <size_t dim>
void check_n_pts(const Octree<dim>& cell) 
{
    REQUIRE(n_pts(cell) == cell.data.indices.size());
    for (size_t c = 0; c < Octree<dim>::split; c++) {
        if (cell.children[c] == nullptr) {
            continue;
        }
        if (cell.children[c]->is_leaf()) {
            continue;
        }
        check_n_pts(*cell.children[c]);
    }
}

TEST_CASE("check octree cell counts", "[octree]") 
{
    auto pts = random_pts<3>(1000);
    auto oct = make_octree(pts, 4);
    check_n_pts(oct);
}

TEST_CASE("check octree cell counts for degenerate line", "[octree]") 
{
    size_t n = 100;
    std::vector<Vec<double,2>> pts;
    for (size_t i = 0; i < n; i++) {
        pts.push_back({static_cast<double>(i), 0.0});
    }
    auto oct = make_octree(pts, 26);
    CHECK(n_pts(oct) == n);
}

TEST_CASE("make octree with two identical points", "[octree]") 
{
    std::vector<Vec<double,3>> es{{1.0, 2.0, 0.0}, {1.0, 2.0, 0.0}};
    auto oct = make_octree(es, 1);
}

TEST_CASE("make octree with two very similar points", "[octree]") 
{
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {1.0, 2.0 - 1e-20, 0.0}, {0.0, 0.0, 0.0}
    };
    auto oct = make_octree(es, 1);
}

template <size_t dim>
size_t count_children(const Octree<dim>& cell) 
{
    if (cell.is_leaf()) {
        return 1;
    }

    size_t n_c = 0;
    for (auto& c: cell.children) {
        if (c == nullptr) {
            continue;
        }
        n_c += 1 + count_children(*c);
    }
    return n_c;
}

TEST_CASE("count children for degenerate line octree", "[octree]") 
{
    size_t n = 10;
    std::vector<Vec<double,2>> pts;     
    for (size_t i = 0; i < n; i++) {
        pts.push_back({static_cast<double>(i), 0});
        pts.push_back({static_cast<double>(i + 1), 0});
    }

    auto oct = make_octree(pts, 1);
    REQUIRE(count_children(oct) == 31); 
}

template <size_t dim>
std::set<size_t> check_indices_unique(const Octree<dim>& oct,
    const std::set<size_t>& indices)
{
    auto new_set = indices;
    auto ret = new_set.emplace(oct.data.index);
    CHECK(ret.second);
    for (auto& c: oct.children) {
        if (c == nullptr) {
            continue;
        }
        new_set = check_indices_unique(*c, new_set);
    }
    return new_set;
}

TEST_CASE("check cells have unique indices", "[octree]") 
{
    auto pts = random_pts<3>(1000);
    auto oct = make_octree(pts, 4);
    check_indices_unique(oct, std::set<size_t>{});
}

TEST_CASE("non zero ball radius", "[octree]")
{
    auto balls = random_balls<3>(10, 0.1);
    auto oct = make_octree(balls, 1);
}

template <size_t dim>
void check_true_bounds_contain_balls(const Octree<dim>& oct,
    const std::vector<Ball<dim>>& balls)
{
    for (size_t i = 0; i < oct.data.indices.size(); i++) {
        auto ball_idx = oct.data.indices[i];
        REQUIRE(oct.data.true_bounds.in_box(balls[ball_idx]));
    }
    for (auto& c: oct.children) {
        if (c == nullptr) {
            continue;
        }
        check_true_bounds_contain_balls(*c, balls);
    }
}

TEST_CASE("test true bounds non-zero radii", "[octree]")
{
    auto balls = random_balls<3>(100, 0.05);
    auto oct = make_octree(balls, 1);
    check_true_bounds_contain_balls(oct, balls);
}

TEST_CASE("impossible subdivision", "[octree]")
{
    std::vector<Ball<2>> bs{
        {{0, 0}, 0.5},
        {{0, 0}, 1.0}
    };
    /* auto oct = make_octree(bs, 1); */
}
