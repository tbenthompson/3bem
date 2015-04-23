#include "UnitTest++.h"
#include "octree.h"
#include "util.h"
#include "vec_ops.h"
#include "test_shared.h"
using namespace tbem;

std::vector<Vec<double,3>> three_pts() {
    return {{1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}};
}
TEST(BoundingBox) {
    auto es = three_pts();
    auto bb = Box<3>::bounding_box(es);
    double hw[] = {1.0, 2.0, 3.0};
    double center[] = {0.0, 0.0, 0.0};
    CHECK_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    CHECK_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
}

TEST(DegenerateBoundingBox) {
    std::vector<Vec<double,3>> es = {{0.0, 0.0, 0.0}};
    auto bb = Box<3>::bounding_box(es);
    double hw[] = {0.0, 0.0, 0.0};
    double center[] = {0.0, 0.0, 0.0};
    CHECK_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    CHECK_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
}

TEST(MakeChildIdx) {
    CHECK_EQUAL(make_child_idx<3>(0), (Vec3<size_t>{0,0,0}));
    CHECK_EQUAL(make_child_idx<3>(2), (Vec3<size_t>{0,1,0}));
    CHECK_EQUAL(make_child_idx<3>(7), (Vec3<size_t>{1,1,1}));
    CHECK_EQUAL(make_child_idx<2>(1), (Vec2<size_t>{0,1}));
    CHECK_EQUAL(make_child_idx<2>(3), (Vec2<size_t>{1,1}));
}

TEST(MakeChildBox) {
    Box<3> box{{0,1,0}, {2,2,2}};
    auto child = box.get_subcell({1,0,1});
    CHECK_ARRAY_CLOSE(child.center, (Vec<double,3>{1,0,1}), 3, 1e-14);
    CHECK_ARRAY_CLOSE(child.half_width, (Vec<double,3>{1,1,1}), 3, 1e-14);
}

TEST(InBox) {
    Box<3> box{{0,1,0}, {2,2,2}};
    CHECK(box.in_box({0,2,0}, {false, false, false}));
    CHECK(!box.in_box({0,4,0}, {false, false, false}));
}

TEST(InBoxInclusive) {
    Box<3> box{{0,1,0}, {2,2,2}};
    CHECK(box.in_box({0,3,0}, {false, true, false}));
    CHECK(!box.in_box({0,3,0}, {false, false, false}));
}

TEST(CheckLawOfLargeNumbers) {
    int n = 1000000;
    auto pts = random_pts<3>(n);
    TIC
    auto tree = build_octree(pts, 100);
    TOC("Build " + std::to_string(n));
    for (size_t i = 0; i < 8; i++) {
        int n_pts = tree.children[i]->data.indices.size();
        int diff = abs(n_pts - (n / 8));
        CHECK(diff < (n / 32));
    }
}

TEST(SmallOctree) {
    auto es = three_pts();
    auto oct = build_octree(es, 4);
    CHECK_EQUAL(oct.data.level, 0);
    for (size_t i = 0; i < 8; i++) {
        CHECK(oct.children[i] == nullptr);
    }
    CHECK_EQUAL(oct.data.indices.size(), 3);
}

TEST(NotOneLevel) {
    auto pts = random_pts<3>(1000);
    auto oct = build_octree(pts, 4);
    CHECK_EQUAL(oct.data.indices.size(), 1000);
    CHECK_EQUAL(oct.children[0]->data.level, 1);
}

TEST(OctreePlane) {
    std::vector<Vec<double,2>> pts;
    for (size_t i = 0; i < 100; i++) {
        pts.push_back({static_cast<double>(i), 0.0});
    }
    auto oct = build_octree(pts, 1);
    CHECK(!oct.is_leaf());
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
    CHECK_EQUAL(n_pts(cell), cell.data.indices.size());
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

TEST(SumPointsRandom) {
    auto pts = random_pts<3>(1000);
    auto oct = build_octree(pts, 4);
    check_n_pts(oct);
}

TEST(SumPointsPlane) {
    size_t n = 100;
    std::vector<Vec<double,2>> pts;
    for (size_t i = 0; i < n; i++) {
        pts.push_back({static_cast<double>(i), 0.0});
    }
    auto oct = build_octree(pts, 26);
    CHECK(n_pts(oct) == n);
}

TEST(OctreeIdenticalPoints) {
    std::vector<Vec<double,3>> es{{1.0, 2.0, 0.0}, {1.0, 2.0, 0.0}};
    auto oct = build_octree(es, 1);
}

TEST(OctreeVerySimilarPoints) {
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

TEST(OctreeLine) {
    size_t n = 10;
    std::vector<Vec<double,2>> pts;     
    for (size_t i = 0; i < n; i++) {
        pts.push_back({static_cast<double>(i), 0});
        pts.push_back({static_cast<double>(i + 1), 0});
    }

    auto oct = build_octree(pts, 1);
    CHECK_EQUAL(count_children(oct), 31); 
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
