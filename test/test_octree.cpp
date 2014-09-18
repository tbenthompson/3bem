#include "UnitTest++.h"
#include "octree.h"
#include "mpi.h"
#include "test_shared.h"

TEST(BoundingBox) {
    auto es = three_pts();
    auto bb = bounding_box(es);
    double min[] = {-1.0, -2.0, -3.0};
    double max[] = {1.0, 2.0, 3.0};
    CHECK_ARRAY_CLOSE(((bb.center + bb.half_width) / 1.001).loc, max, 3, 1e-14);
    CHECK_ARRAY_CLOSE(((bb.center - bb.half_width) / 1.001).loc, min, 3, 1e-14);
}

TEST(Naturals) {
    auto nats5 = naturals(3, 8);
    double correct[] = {3, 4, 5, 6, 7};
    CHECK_ARRAY_EQUAL(nats5, correct, 5);
}

TEST(BigBig2) {
    //TODO: Replace with an mpi_setup procedure.
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    const int n = (int)(1e4 / p);
    std::array<std::vector<double>,3> es =
        {random_list(n), random_list(n), random_list(n)};
    TIC
    Octree octree(es, 5); 
    TOC("Octree assembly");
    // unsigned int total_indices = check_invariant(octree.root);
    // CHECK(total_indices == octree.elements->size());
}

TEST(MortonEncode) {
    int split = 2;
    //TODO: This could be a property based test.
    uint64_t center = morton_encode(0, 0, split);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                uint64_t val = morton_encode(k, j, i); 
                if (i < split) {
                    CHECK(val < center);
                } else if (i == split && j == 0 && k == 0) {
                    CHECK(val == center); 
                } else {
                    CHECK(val > center);
                }
            }
        }
    }
}

TEST(MortonEncodeCrossLevels) {
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
