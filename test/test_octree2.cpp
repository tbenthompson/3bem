#include "UnitTest++.h"
#include "octree2.h"
#include "geom.h"
#include "test_shared.h"

TEST(BigBig2) {
    auto es = random_pts<3>(10000000);
    TIC
    Octree2<Vec<3>, 3> octree(es, 7);
    TOC("Octree2 assembly");
    // unsigned int total_indices = check_invariant(octree.root);
    // CHECK(total_indices == octree.elements->size());
}

TEST(ToOctreeSpace) {
    const double n = 10;
    const double sub_n = 10;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < sub_n; j++) {
            double loc = (i / n) + j / (sub_n * n);
            int val = to_octree_space(loc, 0.5, 0.5, n);
            CHECK(val == i);
        }
    }

    // The function fails the edge case at the boundary. Read the
    // note [interior] in the octree header
    CHECK(to_octree_space(1.0, 0.5, 0.5, 16) == 16);
}
