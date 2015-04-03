#include "UnitTest++.h"
#include "constraint_builder.h"
#include "vectorx.h"
#include "mesh_gen.h"
#include "test_shared.h"

using namespace tbem;

Mesh<2> disjoint_mesh() {
    std::vector<Facet<2>> facets;
    for (int i = 0; i < 10; i++) {
        double val = 2 * i;
        facets.push_back({{
            {{val, -val}}, {{val + 1, -val - 1}}
        }});
    }
    return {facets};
}

Mesh<2> connected_mesh() {
    return line_mesh({{0, 0}}, {{0, 1}}).refine_repeatedly(2);
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
    int n_verts = m.n_dofs();
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
    CHECK_EQUAL(overlaps.size(), m.n_facets() - 1);

    int n_verts = m.n_dofs();
    for (int i = 1; i < n_verts - 2; i += 2) {
        auto v_it = m.begin() + i;
        CHECK(overlaps.find(v_it)->second == v_it + 1);
    }
}

TEST(ConstraintMesh) {
    auto sphere = sphere_mesh({{0, 0, 0}}, 1, 2);
    auto continuity = mesh_continuity(sphere.begin());
    auto constraints = convert_to_constraints(continuity);
    auto matrix = from_constraints(constraints);
    CHECK_EQUAL(3 * sphere.facets.size(), 384);
    CHECK_EQUAL(matrix.size(), 318);
    auto my_c = matrix.begin()->second.terms;
    CHECK_EQUAL(my_c[0].weight, 1);
}

TEST(RectMeshContinuity) {
    auto zplane = rect_mesh({{-1, -1, 0}}, {{1, -1, 0}}, {{1, 1, 0}}, {{-1, 1, 0}})
        .refine_repeatedly(1);
    auto continuity = mesh_continuity(zplane.begin());
    auto constraints = convert_to_constraints(continuity);
    CHECK_EQUAL(constraints.size(), 29);
}

TEST(CutWithDiscontinuity) {
    auto zplane = rect_mesh({{-1, -1, 0}}, {{1, -1, 0}}, {{1, 1, 0}}, {{-1, 1, 0}})
        .refine_repeatedly(1);
    auto xplane = rect_mesh({{0, -1, -1}}, {{0, 1, -1}}, {{0, 1, 1}}, {{0, -1, 1}})
        .refine_repeatedly(1);
    auto continuity = mesh_continuity(zplane.begin());
    auto cut_cont = cut_at_intersection(continuity, zplane.begin(), xplane.begin());
    auto constraints = convert_to_constraints(cut_cont);
    CHECK_EQUAL(constraints.size(), 16);
}

TEST(GetReducedToCountTheNumberOfVerticesOnASphereApproximation) {
    auto sphere = sphere_mesh({{0, 0, 0}}, 1, 0);
    auto continuity = mesh_continuity(sphere.begin());
    auto constraints = convert_to_constraints(continuity);
    auto matrix = from_constraints(constraints);
    std::vector<double> all(3 * sphere.facets.size(), 1.0);
    auto reduced = condense_vector(matrix, all);
    CHECK_EQUAL(reduced.size(), 6);
    CHECK_ARRAY_EQUAL(reduced, (std::vector<double>(6, 4.0)), 6);
}

TEST(ImposeNeighborBCs) {
    auto m0 = disjoint_mesh();
    auto m1 = connected_mesh();
    std::vector<double> bcs(m1.n_dofs(), 2.33);
    auto constraints = form_neighbor_bcs(m1.begin(), m0.begin(), bcs); 
    CHECK_EQUAL(constraints.size(), 1);
    CHECK_EQUAL(constraints[0].terms.size(), 1);
    CHECK_EQUAL(constraints[0].terms[0], (LinearTerm{0, 1}));
    CHECK_EQUAL(constraints[0].rhs, 2.33);
}

TEST(InterpolateBCConstraints) {
    auto m0 = line_mesh({-1, 0}, {1, 0});
    auto cs = interpolate_bc_constraints<2>(m0, range(m0.n_dofs()),
        [](const Vec<double,2>& x) {
            return x[0] + 1.0;
        });
    CHECK_EQUAL(cs.size(), 2);
    CHECK_EQUAL(cs[0].terms[0], (LinearTerm{0, 1}));
    CHECK_EQUAL(cs[0].rhs, 0.0);
    CHECK_EQUAL(cs[1].terms[0], (LinearTerm{1, 1}));
    CHECK_EQUAL(cs[1].rhs, 2.0);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
    // return RunOneTest("FindOverlappingVerticesDifferentMeshes");
}
