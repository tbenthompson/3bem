#include "UnitTest++.h"
#include "octree.h"
#include "test_shared.h"

TEST(BoundingBox) {
    auto es = three_pts();
    auto bb = bounding_box(es);
    double hw[] = {1.001, 2.002, 3.003};
    double center[] = {0.0, 0.0, 0.0};
    CHECK_ARRAY_CLOSE(bb.half_width, hw, 3, 1e-14);
    CHECK_ARRAY_CLOSE(bb.center, center, 3, 1e-14);
    for (int d = 0; d < 3; d++) {
        CHECK_CLOSE(bb.center[d] - bb.half_width[d], bb.min_corner[d], 1e-12);
        CHECK_CLOSE(bb.center[d] + bb.half_width[d], bb.max_corner[d], 1e-12);
    }
}

TEST(GetChildBox) {
    Box parent = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 
                  {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}};
    Box child = get_child_box({1, 0, 1}, parent);
    double correct[2][3] = {
        {0.5, -0.5, 0.5}, //center
        {0.5, 0.5, 0.5} //half_width
    };
    CHECK_ARRAY_CLOSE(child.center, correct[0], 3, 1e-6);
    CHECK_ARRAY_CLOSE(child.half_width, correct[1], 3, 1e-6);
    for (int d = 0; d < 3; d++) {
        CHECK_CLOSE(child.center[d] - child.half_width[d], child.min_corner[d], 1e-12);
        CHECK_CLOSE(child.center[d] + child.half_width[d], child.max_corner[d], 1e-12);
    }
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

TEST(BigBig2) {
    const int n = 1e4;
    auto pts = random_pts(n);
    TIC
    Octree octree(pts, 500000); 
    TOC("Octree assembly");
    // unsigned int total_indices = check_invariant(octree.root);
    // CHECK(total_indices == octree.elements->size());
}

