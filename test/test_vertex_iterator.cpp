#include "UnitTest++.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "vertex_iterator.h"

using namespace tbem;

struct PremadeMesh {
    const Mesh<3> m;
    PremadeMesh():
        m(sphere_mesh({0, 0, 0}, 1.0))
    { }
};

TEST_FIXTURE(PremadeMesh, VertexIteratorBegin) {
    auto iter = m.begin();
    CHECK_EQUAL(iter.facet_idx, 0);
    CHECK_EQUAL(iter.vertex_idx, 0);
}

TEST_FIXTURE(PremadeMesh, VertexIteratorEnd) {
    auto iter = m.end();
    CHECK_EQUAL(iter.facet_idx, m.facets.size());
    CHECK_EQUAL(iter.vertex_idx, 0);
}

TEST_FIXTURE(PremadeMesh, NextVert) {
    auto iter = m.begin();
    ++iter;
    CHECK_EQUAL(iter.facet_idx, 0);
    CHECK_EQUAL(iter.vertex_idx, 1);
}

TEST_FIXTURE(PremadeMesh, NextVertNextFacet) {
    auto iter = m.begin();
    ++iter; ++iter; ++iter;
    CHECK_EQUAL(iter.facet_idx, 1);
    CHECK_EQUAL(iter.vertex_idx, 0);
}

TEST_FIXTURE(PremadeMesh, Equality) {
    CHECK(m.begin() == m.begin());
}

TEST_FIXTURE(PremadeMesh, Inequality) {
    auto iter = m.begin();
    auto iter2 = m.begin();
    ++iter2;
    CHECK(iter != iter2);
}

TEST_FIXTURE(PremadeMesh, InequalityEnd) {
    CHECK(m.begin() != m.end());
}

TEST_FIXTURE(PremadeMesh, CountVerticesInMesh) {
    int n_verts = 0;
    for (auto iter = m.begin(); iter != m.end(); ++iter) {
        n_verts++;
    }
    CHECK_EQUAL(n_verts, m.facets.size() * 3);
}

TEST_FIXTURE(PremadeMesh, Dereference) {
    auto iter = m.begin(); 
    auto vert = *iter;
    CHECK_EQUAL(vert, m.facets[0].vertices[0]);
}

TEST_FIXTURE(PremadeMesh, DereferenceNonInitial) {
    auto iter = m.begin(); 
    iter.facet_idx = 2;
    iter.vertex_idx = 2;
    auto vert = *iter;
    CHECK_EQUAL(vert, m.facets[2].vertices[2]);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
