#include "UnitTest++.h"
#include "octree.h"
#include "test_shared.h"

TEST(BoundingBox) {
    auto es = three_pts();
    auto bb = bounding_box(es);
    double min[] = {-1.0, -2.0, -3.0};
    double max[] = {1.0, 2.0, 3.0};
    CHECK_ARRAY_CLOSE((bb.center + bb.half_width).loc, max, 3, 1e-14);
    CHECK_ARRAY_CLOSE((bb.center - bb.half_width).loc, min, 3, 1e-14);
}

TEST(Naturals) {
    auto nats5 = naturals(5);
    double correct[] = {0, 1, 2, 3, 4};
    CHECK_ARRAY_EQUAL(nats5, correct, 5);
}

TEST(BigBig2) {
    const int n = (int)1e5;
    std::array<std::vector<double>,3> es =
        {random_list(n), random_list(n), random_list(n)};
    TIC
    Octree octree(es, 7);
    TOC("Octree assembly");
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
