#include "UnitTest++.h"
#include "octree.h"
#include "test_shared.h"

std::unique_ptr<std::vector<Vec<3> > > three_pts() {
    std::unique_ptr<std::vector<Vec<3> > > es(new std::vector<Vec<3> >);
    es->push_back(Vec<3>({1.0, 0.0, 3.0}));
    es->push_back(Vec<3>({0.0, 2.0, -3.0}));
    es->push_back(Vec<3>({-1.0, -2.0, 0.0}));
    return es;
}

unsigned int check_invariant(OctreeNode<Vec<3>, 3> node) {
    if (node.children.size() == 0) {
        return node.indices.size();
    }
    unsigned int contained_indices = 0;
    for (auto n: node.children) {
        contained_indices += check_invariant(n);
    }
    CHECK(contained_indices > 0);
    return contained_indices;
}

TEST(BigBig) {
    auto es = random_pts<1>(100000);
    TIC
    Octree<Vec<1>, 1> octree(es, 10);
    TOC("Octree assembly");
    // unsigned int total_indices = check_invariant(octree.root);
    // CHECK(total_indices == octree.elements->size());
}

TEST(CreateOnlyLeaf) {
    auto es = three_pts();
    Octree<Vec<3>, 3> octree(es, 4);
    if (es) {
        CHECK(false);
    }
    CHECK_ARRAY_EQUAL(octree.root.indices, naturals(3), 3);
}

TEST(OneLevel) {
    auto es = three_pts();
    Octree<Vec<3>, 3> octree(es, 1);
    CHECK(octree.root.indices.size() == 0);
    CHECK(octree.root.children[0].indices[0] == 0);
    CHECK(octree.root.children[4].indices[0] == 1);
    CHECK(octree.root.children[3].indices[0] == 2);
}

double centers[8][3] = {
    {0.5, 0.5, 0.5},
    {-0.5, 0.5, 0.5},
    {0.5, -0.5, 0.5},
    {-0.5, -0.5, 0.5},
    {0.5, 0.5, -0.5},
    {-0.5, 0.5, -0.5},
    {0.5, -0.5, -0.5},
    {-0.5, -0.5, -0.5}
};

TEST(FindOctant) {
    auto zero = Vec<3>({0.0, 0.0, 0.0});
    for(int i = 0; i < 8; ++i) {
        Vec<3> vec = {{centers[i][0], centers[i][1], centers[i][2]}};
        CHECK(find_octant(vec, zero) == i);
    }
}

TEST(BoundingBox) {
    std::vector<int> idxs = {0, 1, 2};
    auto bb = bounding_box_subset(*three_pts(), idxs);
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
