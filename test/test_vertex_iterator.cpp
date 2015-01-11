#include "UnitTest++.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "vertex_iterator.h"

using namespace tbem;

struct PremadeMesh3D {
    const Mesh<3> m;
    PremadeMesh3D():
        m(sphere_mesh({0, 0, 0}, 1.0, 0))
    { }
};

struct PremadeMesh2D {
    const Mesh<2> m;
    PremadeMesh2D():
        m(circle_mesh({0, 0}, 1.0, 0))
    { }
};

TEST_FIXTURE(PremadeMesh3D, VertexIteratorEnd) {
    auto iter = m.end();
    CHECK_EQUAL(iter.facet_idx, m.n_facets());
    CHECK_EQUAL(iter.vertex_idx, 0);
}

TEST_FIXTURE(PremadeMesh3D, VertexIteratorIsNotEnd) {
    auto iter = m.begin();
    CHECK(!iter.is_end());
    CHECK((iter + m.n_dofs()).is_end());
}

TEST_FIXTURE(PremadeMesh3D, NextVertNextFacet) {
    auto iter = m.begin();
    ++iter;
    CHECK_EQUAL(iter.facet_idx, 0);
    CHECK_EQUAL(iter.vertex_idx, 1);
}

TEST_FIXTURE(PremadeMesh3D, EqualityInequality) {
    auto iter = m.begin();
    auto iter2 = m.begin();
    CHECK(iter == iter2);
    ++iter2;
    CHECK(iter != iter2);
}

TEST_FIXTURE(PremadeMesh3D, InequalityDifferentMeshes) {
    auto m2 = sphere_mesh({0, 0, 0}, 1.0, 0);
    auto iter = m.begin();
    auto iter2 = m2.begin();
    CHECK(iter != iter2);
}

TEST_FIXTURE(PremadeMesh3D, CountVerticesInMesh) {
    int n_verts = 0;
    for (auto iter = m.begin(); iter != m.end(); ++iter) {
        n_verts++;
    }
    CHECK_EQUAL(n_verts, m.n_dofs());
}

TEST_FIXTURE(PremadeMesh2D, DereferenceNonInitial) {
    auto iter = m.begin(); 
    iter.facet_idx = 2;
    iter.vertex_idx = 1;
    auto vert = *iter;
    CHECK_EQUAL(vert, m.facets[2].vertices[1]);
}

TEST_FIXTURE(PremadeMesh3D, AbsoluteIndex3D) {
    VertexIterator<3> iter(m, 1, 1);
    CHECK_EQUAL(iter.absolute_index(), 4);
}

TEST_FIXTURE(PremadeMesh3D, HashEquality) {
    VertexIterator<3> iter(m, 1, 1);
    VertexIterator<3> iter2(m, 1, 1);
    HashVertexIterator<3> hash;
    CHECK(hash(iter) == hash(iter2));
}

TEST_FIXTURE(PremadeMesh2D, HashInequality) {
    VertexIterator<2> iter(m, 1, 1);
    VertexIterator<2> iter2(m, 0, 1);
    HashVertexIterator<2> hash;
    CHECK(hash(iter) != hash(iter2));
}

TEST_FIXTURE(PremadeMesh3D, StepForward) {
    auto iter = m.begin();
    iter += 4;
    CHECK_EQUAL(iter.facet_idx, 1);
    CHECK_EQUAL(iter.vertex_idx, 1);
    iter += 1;
    CHECK_EQUAL(iter.facet_idx, 1);
    CHECK_EQUAL(iter.vertex_idx, 2);
}

TEST_FIXTURE(PremadeMesh3D, PlusOperator) {
    auto iter = m.begin() + 4;
    CHECK_EQUAL(iter.facet_idx, 1);
    CHECK_EQUAL(iter.vertex_idx, 1);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
