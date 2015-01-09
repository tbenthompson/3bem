#include "UnitTest++.h"
#include "continuity_builder.h"
#include "mesh_gen.h"
#include "shared.h"

using namespace tbem;

Mesh<2> disjoint_mesh() {
    std::vector<Facet<2>> facets;
    for (int i = 0; i < 10; i++) {
        double val = 2 * i;
        facets.push_back(Facet<2>{{{
            {val, -val}, {val + 1, -val - 1}
        }}});
    }
    return {facets, false, nullptr};
}

Mesh<2> connected_mesh() {
    return line_mesh({0, 0}, {0, 1}).refine_repeatedly(2);
}

TEST(FindOverlappingVerticesDifferentMeshes) {
    auto m0 = disjoint_mesh();
    auto m1 = connected_mesh();
    auto overlaps = find_overlapping_vertices(m0.begin(), m1.begin());
    CHECK_EQUAL(overlaps.size(), 1);
    CHECK(overlaps.find(m0.begin())->second == m1.begin());
}

TEST(FindOverlappingVertices) {
    auto m = disjoint_mesh();
    auto overlaps = find_overlapping_vertices(m.begin(), m.begin());
    int n_verts = 2 * m.facets.size();
    CHECK_EQUAL(overlaps.size(), n_verts);
    for (int i = 0; i < n_verts; i++) {
        auto v_it = m.begin() + i;
        CHECK(overlaps.find(v_it)->second == v_it);
    }
}

TEST(FindOverlappingVerticesSameMeshDisjoint) {
    auto m = disjoint_mesh();
    auto overlaps = find_overlapping_vertices_same_mesh(m.begin());
    CHECK_EQUAL(overlaps.size(), 0);
}

TEST(FindOverlappingVerticesSameMeshConnected) {
    auto m = connected_mesh();
    auto overlaps = find_overlapping_vertices_same_mesh(m.begin());
    CHECK_EQUAL(overlaps.size(), m.facets.size() - 1);

    int n_verts = 2 * m.facets.size();
    for (int i = 1; i < n_verts - 2; i += 2) {
        auto v_it = m.begin() + i;
        CHECK(overlaps.find(v_it)->second == v_it + 1);
    }
}

TEST(ContinuityBuilderBuild) {
    auto m = connected_mesh();    
    ContinuityBuilder<2> cb(m.begin());
    auto constraints = cb.build();
    CHECK_EQUAL(constraints.size(), m.facets.size() - 1);
}

TEST(ConstraintMesh) {
    auto sphere = sphere_mesh({0, 0, 0}, 1).refine_repeatedly(2);
    ContinuityBuilder<3> cb(sphere.begin());
    auto constraints = cb.build();
    auto matrix = ConstraintMatrix::from_constraints(constraints);
    CHECK_EQUAL(3 * sphere.facets.size(), 384);
    CHECK_EQUAL(matrix.map.size(), 318);
    auto my_c = matrix.map.begin()->second.terms;
    CHECK_EQUAL(my_c[0].weight, 1);
}

TEST(RectMeshContinuity) {
    auto zplane = rect_mesh({-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0})
        .refine_repeatedly(1);
    ContinuityBuilder<3> cb(zplane.begin());
    auto constraints = cb.build();
    CHECK_EQUAL(constraints.size(), 29);
}

TEST(CutWithDiscontinuity) {
    auto zplane = rect_mesh({-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0})
        .refine_repeatedly(1);
    auto xplane = rect_mesh({0, -1, -1}, {0, 1, -1}, {0, 1, 1}, {0, -1, 1})
        .refine_repeatedly(1);
    ContinuityBuilder<3> cb(zplane.begin());
    cb.apply_discontinuities(xplane.begin());
    auto constraints = cb.build();
    CHECK_EQUAL(constraints.size(), 16);
}

TEST(GetReducedToCountTheNumberOfVerticesOnASphereApproximation) {
    auto sphere = sphere_mesh({0, 0, 0}, 1).refine_repeatedly(0);
    ContinuityBuilder<3> cb(sphere.begin());
    auto constraints = cb.build();
    auto matrix = ConstraintMatrix::from_constraints(constraints);
    std::vector<double> all(3 * sphere.facets.size(), 1.0);
    auto reduced = matrix.get_reduced(all);
    CHECK_EQUAL(reduced.size(), 6);
    CHECK_ARRAY_EQUAL(reduced, (std::vector<double>(6, 4.0)), 6);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
    // return RunOneTest("FindOverlappingVerticesDifferentMeshes");
}
