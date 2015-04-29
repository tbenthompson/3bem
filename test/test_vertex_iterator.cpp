#include "catch.hpp"
#include "mesh.h"
#include "mesh_gen.h"
#include "vertex_iterator.h"
#include "vec_ops.h"

using namespace tbem;

TEST_CASE("Vertex iterator 3D", "[vertex_iterator]") 
{
    const Mesh<3> m = sphere_mesh({0, 0, 0}, 1.0, 0);

    SECTION("VertexIteratorEnd") {
        auto iter = m.end();
        REQUIRE(iter.facet_idx == m.n_facets());
        REQUIRE(iter.vertex_idx == 0);
    }

    SECTION("VertexIteratorIsNotEnd") {
        auto iter = m.begin();
        REQUIRE(!iter.is_end());
        REQUIRE((iter + m.n_dofs()).is_end());
    }

    SECTION("NextVertNextFacet") {
        auto iter = m.begin();
        ++iter;
        REQUIRE(iter.facet_idx == 0);
        REQUIRE(iter.vertex_idx == 1);
    }

    SECTION("EqualityInequality") {
        auto iter = m.begin();
        auto iter2 = m.begin();
        REQUIRE(iter == iter2);
        ++iter2;
        REQUIRE(iter != iter2);
    }

    SECTION("InequalityDifferentMeshes") {
        auto m2 = sphere_mesh({0, 0, 0}, 1.0, 0);
        auto iter = m.begin();
        auto iter2 = m2.begin();
        REQUIRE(iter != iter2);
    }

    SECTION("CountVerticesInMesh") {
        int n_verts = 0;
        for (auto iter = m.begin(); iter != m.end(); ++iter) {
            n_verts++;
        }
        REQUIRE(n_verts == m.n_dofs());
    }

    SECTION("AbsoluteIndex3D") {
        VertexIterator<3> iter(m, 1, 1);
        REQUIRE(iter.absolute_index() == 4);
    }

    SECTION("HashEquality") {
        VertexIterator<3> iter(m, 1, 1);
        VertexIterator<3> iter2(m, 1, 1);
        HashVertexIterator<3> hash;
        REQUIRE(hash(iter) == hash(iter2));
    }

    SECTION("StepForward") {
        auto iter = m.begin();
        iter += 4;
        REQUIRE(iter.facet_idx == 1);
        REQUIRE(iter.vertex_idx == 1);
        iter += 1;
        REQUIRE(iter.facet_idx == 1);
        REQUIRE(iter.vertex_idx == 2);
    }

    SECTION("PlusOperator") {
        auto iter = m.begin() + 4;
        REQUIRE(iter.facet_idx == 1);
        REQUIRE(iter.vertex_idx == 1);
    }
}

TEST_CASE("Vertex iterator 2D", "[vertex_iterator]") 
{
    const auto m = circle_mesh({0, 0}, 1.0, 0);

    SECTION("DereferenceNonInitial") {
        auto iter = m.begin(); 
        iter.facet_idx = 2;
        iter.vertex_idx = 1;
        auto vert = *iter;
        REQUIRE(vert == m.facets[2][1]);
    }

    SECTION("HashInequality") {
        VertexIterator<2> iter(m, 1, 1);
        VertexIterator<2> iter2(m, 0, 1);
        HashVertexIterator<2> hash;
        REQUIRE(hash(iter) != hash(iter2));
    }
}
