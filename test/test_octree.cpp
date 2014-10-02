#include "UnitTest++.h"
#include "octree.h"
#include "test_shared.h"

TEST(BoundingBox) {
    auto es = three_pts();
    auto bb = bounding_box(es, 0, 3);
    double hw[] = {1.0, 2.0, 3.0};
    double center[] = {0.0, 0.0, 0.0};
    CHECK_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    CHECK_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
    for (int d = 0; d < 3; d++) {
        CHECK_CLOSE(bb.center[d] - bb.half_width[d], bb.min_corner[d], 1e-12);
        CHECK_CLOSE(bb.center[d] + bb.half_width[d], bb.max_corner[d], 1e-12);
    }
}

TEST(DegenerateBoundingBox) {
    std::vector<std::array<double,3>> pts = {{0.0, 0.0, 0.0}};
    auto bb = bounding_box(pts, 0, 1);
    double hw[] = {0.0, 0.0, 0.0};
    double center[] = {0.0, 0.0, 0.0};
    CHECK_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    CHECK_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
}

TEST(ToOctreeSpace) {
    const double n = 10;
    const double sub_n = 10;
    //TODO: This could be a property based test.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < sub_n; j++) {
            double loc = (i / n) + j / (sub_n * n);
            if (j == 0) {
                // Behavior is undefined on the boundary. This is okay.
                loc += 0.001;
            }
            int val = to_octree_space(loc, 0.5, 0.5, n);
            CHECK(val == i);
        }
    }

    // The function fails the edge case at the boundary. Read the
    // note [interior] in the octree header
    CHECK(to_octree_space(1.0, 0.5, 0.5, 16) == 16);
}

TEST(ToOctreeSpace1Cell) {
    int val = to_octree_space(0.3, 0.5, 0.5, 1);
    CHECK(val == 0);
    val = to_octree_space(0.8, 0.5, 0.5, 1);
    CHECK(val == 0);
}

void check_nonoverlapping(OctreeCell& cell, Octree tree) {

    bool first = true;
    OctreeCell* prev_cell = NULL;
    OctreeCell* cur_cell = NULL;
    for (int i = 0; i < 8; i++) {
        if (cell.children[i] == -1) {
            continue;
        }
        cur_cell = &tree.cells[cell.children[i]];
        if (first) {
            CHECK(cell.begin == cur_cell->begin);
            first = false;
        } else {
            CHECK(prev_cell->end == cur_cell->begin);
        }
        prev_cell = cur_cell;
    }
    if(&cur_cell != NULL) {
        CHECK(cur_cell->end == cell.end);
    }
}

TEST(NonOverlapping) {
    int n = 100;
    for (int i = 0; i < 10; i++) {
        auto pts = random_pts(n);
        Octree tree(pts, 5);
        check_nonoverlapping(tree.get_root(), tree);
    }
}

void check_children_contained(OctreeCell& cell, Octree tree) {
    for(int i = 0; i < 8; i++) {
        int cur = cell.children[i];
        if (cur == -1) {
            continue;
        }
        // std::cout << cur << std::endl;
        for(int d = 0; d < 3; d++) {
            // std::cout << cell.bounds.min_corner[d] << "  " << tree.cells[cur].bounds.min_corner[d] << std::endl;
            // std::cout << cell.bounds.max_corner[d] << "  " << tree.cells[cur].bounds.max_corner[d] << std::endl;
            CHECK(cell.bounds.min_corner[d] <= tree.cells[cur].bounds.min_corner[d]);
            CHECK(cell.bounds.max_corner[d] >= tree.cells[cur].bounds.max_corner[d]);
        }
        check_children_contained(tree.cells[cur], tree); 
    }
}

TEST(ChildrenContainedInParents) {
    int n = 100;
    auto pts = random_pts(n);
    Octree tree(pts, 5);
    check_children_contained(tree.get_root(), tree);
}

TEST(CheckLawOfLargeNumbers) {
    int n = 100000;
    auto pts = random_pts(n);
    Octree tree(pts, n / 4);
    for (int i = 0; i < 8; i++) {
        int n_pts = tree.cells[i].end - tree.cells[i].begin;
        CHECK(abs(n_pts - (n / 8)) < 1000);
    }
}

TEST(SmallOctree) {
    auto es = three_pts();
    Octree oct(es, 4);
    CHECK(oct.cells.size() == 1);
    CHECK(oct.cells[0].begin == 0);
    CHECK(oct.cells[0].end == 3);
}
