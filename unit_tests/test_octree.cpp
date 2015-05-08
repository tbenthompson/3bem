#include "catch.hpp"
#include "octree.h"
#include "util.h"
#include "vec_ops.h"
#include <set>
using namespace tbem;

std::vector<Vec<double,3>> three_pts() {
    return {{1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}};
}
TEST_CASE("BoundingBox", "[octree]") {
    auto es = three_pts();
    auto bb = Box<3>::bounding_box(es);
    double hw[] = {1.0, 2.0, 3.0};
    double center[] = {0.0, 0.0, 0.0};
    REQUIRE_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    REQUIRE_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
}

TEST_CASE("DegenerateBoundingBox", "[octree]") {
    std::vector<Vec<double,3>> es = {{0.0, 0.0, 0.0}};
    auto bb = Box<3>::bounding_box(es);
    double hw[] = {0.0, 0.0, 0.0};
    double center[] = {0.0, 0.0, 0.0};
    REQUIRE_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    REQUIRE_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
}

TEST_CASE("MakeChildIdx", "[octree]") {
    REQUIRE(make_child_idx<3>(0) == (Vec3<size_t>{0,0,0}));
    REQUIRE(make_child_idx<3>(2) == (Vec3<size_t>{0,1,0}));
    REQUIRE(make_child_idx<3>(7) == (Vec3<size_t>{1,1,1}));
    REQUIRE(make_child_idx<2>(1) == (Vec2<size_t>{0,1}));
    REQUIRE(make_child_idx<2>(3) == (Vec2<size_t>{1,1}));
}

TEST_CASE("MakeChildBox", "[octree]") {
    Box<3> box{{0,1,0}, {2,2,2}};
    auto child = box.get_subcell({1,0,1});
    REQUIRE_ARRAY_CLOSE(child.center, (Vec<double,3>{1,0,1}), 3, 1e-14);
    REQUIRE_ARRAY_CLOSE(child.half_width, (Vec<double,3>{1,1,1}), 3, 1e-14);
}

TEST_CASE("InBox", "[octree]") {
    Box<3> box{{0,1,0}, {2,2,2}};
    REQUIRE(box.in_box({0,2,0}, {false, false, false}));
    REQUIRE(!box.in_box({0,4,0}, {false, false, false}));
}

TEST_CASE("InBoxInclusive", "[octree]") {
    Box<3> box{{0,1,0}, {2,2,2}};
    REQUIRE(box.in_box({0,3,0}, {false, true, false}));
    REQUIRE(!box.in_box({0,3,0}, {false, false, false}));
}

TEST_CASE("CheckLawOfLargeNumbers", "[octree]") {
    int n = 1000000;
    auto pts = random_pts<3>(n);
    //TODO: Make a octree capacity test
    auto tree = build_octree(pts, 100);
    for (size_t i = 0; i < 8; i++) {
        int n_pts = tree.children[i]->data.indices.size();
        int diff = abs(n_pts - (n / 8));
        CHECK(diff < (n / 32));
    }
}

TEST_CASE("SmallOctree", "[octree]") {
    auto es = three_pts();
    auto oct = build_octree(es, 4);
    REQUIRE(oct.data.level == 0);
    for (size_t i = 0; i < 8; i++) {
        REQUIRE(oct.children[i] == nullptr);
    }
    REQUIRE(oct.data.indices.size() == 3);
}

TEST_CASE("NotOneLevel", "[octree]") {
    auto pts = random_pts<3>(1000);
    auto oct = build_octree(pts, 4);
    REQUIRE(oct.data.indices.size() == 1000);
    REQUIRE(oct.children[0]->data.level == 1);
}

TEST_CASE("OctreePlane", "[octree]") {
    std::vector<Vec<double,2>> pts;
    for (size_t i = 0; i < 100; i++) {
        pts.push_back({static_cast<double>(i), 0.0});
    }
    auto oct = build_octree(pts, 1);
    REQUIRE(!oct.is_leaf());
}

template <size_t dim>
size_t n_pts(const Octree<dim>& cell) {
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
void check_n_pts(const Octree<dim>& cell) {
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

TEST_CASE("SumPointsRandom", "[octree]") {
    auto pts = random_pts<3>(1000);
    auto oct = build_octree(pts, 4);
    check_n_pts(oct);
}

TEST_CASE("SumPointsPlane", "[octree]") {
    size_t n = 100;
    std::vector<Vec<double,2>> pts;
    for (size_t i = 0; i < n; i++) {
        pts.push_back({static_cast<double>(i), 0.0});
    }
    auto oct = build_octree(pts, 26);
    CHECK(n_pts(oct) == n);
}

TEST_CASE("OctreeIdenticalPoints", "[octree]") {
    std::vector<Vec<double,3>> es{{1.0, 2.0, 0.0}, {1.0, 2.0, 0.0}};
    auto oct = build_octree(es, 1);
}

TEST_CASE("OctreeVerySimilarPoints", "[octree]") {
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {1.0, 2.0 - 1e-20, 0.0}, {0.0, 0.0, 0.0}
    };
    auto oct = build_octree(es, 1);
}

template <size_t dim>
size_t count_children(const Octree<dim>& cell) {
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

TEST_CASE("OctreeLine", "[octree]") {
    size_t n = 10;
    std::vector<Vec<double,2>> pts;     
    for (size_t i = 0; i < n; i++) {
        pts.push_back({static_cast<double>(i), 0});
        pts.push_back({static_cast<double>(i + 1), 0});
    }

    auto oct = build_octree(pts, 1);
    REQUIRE(count_children(oct) == 31); 
}

template <size_t dim>
std::set<size_t> check_nonoverlapping_indices(const Octree<dim>& oct,
    const std::set<size_t>& indices)
{
    auto new_set = indices;
    auto ret = new_set.emplace(oct.data.index);
    CHECK(ret.second);
    for (auto& c: oct.children) {
        if (c == nullptr) {
            continue;
        }
        new_set = check_nonoverlapping_indices(*c, new_set);
    }
    return new_set;
}

TEST_CASE("NonOverlappingIndices", "[octree]") {
    auto pts = random_pts<3>(1000);
    auto oct = build_octree(pts, 4);
    check_nonoverlapping_indices(oct, std::set<size_t>{});
}
