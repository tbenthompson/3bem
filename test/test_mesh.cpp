#include "UnitTest++.h"
#include "util.h"
#include "numerics.h"
#include "mesh.h"
#include "mesh_gen.h"

using namespace tbem;

double perimeter(Mesh<2> m) {
    double p = 0;
    for (auto f: m.facets) {
        p += dist(f.vertices[0], f.vertices[1]);
    }
    return p;
}

TEST(Mesh2D) {
    auto m = circle_mesh({0,0}, 1.0);
    CHECK_CLOSE(m.facets[0].vertices[0][0], 1.0, 1e-15);
    auto m2 = m.refine_repeatedly(3);
    CHECK_EQUAL(m2.facets.size(), 32);
    double length = perimeter(m2);
    CHECK_CLOSE(length, 2 * M_PI, 1e-1);
}

TEST(CircleMesh) {
    std::array<double, 2> center = {20.0, 0.0};
    auto src_circle = circle_mesh(center, 19.0).refine_repeatedly(4);
    CHECK_CLOSE(src_circle.facets[16].vertices[0][0], 20.0, 1e-4);
    CHECK_CLOSE(src_circle.facets[32].vertices[0][0], 1.0, 1e-4);
    CHECK_CLOSE(src_circle.facets[48].vertices[0][0], 20.0, 1e-4);
    CHECK_CLOSE(src_circle.facets[0].vertices[0][0], 39.0, 1e-4);
    CHECK_CLOSE(src_circle.facets[16].vertices[0][1], 19.0, 1e-4);
    CHECK_CLOSE(src_circle.facets[32].vertices[0][1], 0.0, 1e-4);
    CHECK_CLOSE(src_circle.facets[48].vertices[0][1], -19.0, 1e-4);
    CHECK_CLOSE(src_circle.facets[0].vertices[0][1], 0.0, 1e-4);
}

double surface_area(Mesh<3> m) {
    double sa = 0;
    for (auto f: m.facets) {
        sa += tri_area(f.vertices);
    }
    return sa;
}

TEST(Mesh3D) {
    auto m = sphere_mesh({0,0,0}, 1.0); 
    auto m2 = m.refine_repeatedly(4);
    double sa = surface_area(m2);
    CHECK_CLOSE(sa, 4 * M_PI, 1e-1);
}

TEST(FacetField) {
    FacetField<double,3> a = {1,2,3};
    CHECK_EQUAL(a.vertices, (Vec<double,3>{1,2,3}));
}

TEST(MeshField) {
    MeshField<double,3> mf{
        { {0.0, 1.0, 2.0}, {2.0, 3.0, 4.0}, {4.0, 5.0, 6.0} },
        false, nullptr
    };
    auto mf_refined = mf.refine();
    CHECK_EQUAL(mf_refined.facets[0].vertices[1], 0.5);
}

TEST(MeshFieldUnion) {
    MeshField<double,3> mf1{ { {0.0, 1.0, 2.0} }, false, nullptr };
    MeshField<double,3> mf2{ { {2.0, 3.0, 4.0}, {4.0, 5.0, 6.0} }, false, nullptr };
    auto mf_combined = MeshField<double,3>::form_union({mf1, mf2});
    CHECK_EQUAL(mf_combined.facets.size(), 3);
    CHECK_EQUAL(mf_combined.facets[1].vertices[2], 4.0);
}

TEST(CallRefineWithNoFacetsToRefine) {
    MeshField<double,3> mf1{ { {0.0, 1.0, 2.0} }, false, nullptr };
    mf1.refine({});
}


int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
