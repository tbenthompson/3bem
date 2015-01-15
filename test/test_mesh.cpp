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

TEST(Refine2DMesh) {
    auto m2 = circle_mesh({0,0}, 1.0, 3);
    double length = perimeter(m2);
    CHECK_EQUAL(m2.n_facets(), 32);
    CHECK_CLOSE(length, 2 * M_PI, 1e-1);
}

TEST(RefineSubset) {
    auto m2 = line_mesh({0,0}, {1,0}).refine({0}).refine({0});
    CHECK_EQUAL(m2.n_facets(), 3);
    CHECK_EQUAL(m2.facets[1].vertices[0], (Vec2<double>{0.25,0.0}));
}


TEST(CircleMesh) {
    std::array<double, 2> center = {20.0, 0.0};
    auto src_circle = circle_mesh(center, 19.0, 4);
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
    auto m = sphere_mesh({0,0,0}, 1.0, 4); 
    double sa = surface_area(m);
    CHECK_CLOSE(sa, 4 * M_PI, 1e-1);
}


int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
