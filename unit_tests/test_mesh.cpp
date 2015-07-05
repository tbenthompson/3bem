#include "catch.hpp"
#include "util.h"
#include "numerics.h"
#include "mesh.h"
#include "mesh_gen.h"

using namespace tbem;

double perimeter(Mesh<2> m) 
{
    double p = 0;
    for (auto f: m.facets) {
        p += dist(f[0], f[1]);
    }
    return p;
}

TEST_CASE("CallRefineWithNoFacetsToRefine", "[mesh]") 
{
    Mesh<3> mf1{{Facet<3>{{0.0, 1.0, 2.0}}}};
    mf1.refine({});
}

TEST_CASE("Refine2DMesh", "[mesh]") 
{
    auto m2 = circle_mesh({0,0}, 1.0, 3);
    double length = perimeter(m2);
    REQUIRE(m2.n_facets() == 32);
    REQUIRE_CLOSE(length, 2 * M_PI, 1e-1);
}

TEST_CASE("RefineSubset", "[mesh]") 
{
    auto m2 = line_mesh({0,0}, {1,0}).refine({0}).refine({0});
    REQUIRE(m2.n_facets() == 4);
    REQUIRE(m2.facets[1][0] == (Vec2<double>{0.25,0.0}));
}


TEST_CASE("CircleMesh", "[mesh]") 
{
    std::array<double, 2> center = {20.0, 0.0};
    auto src_circle = circle_mesh(center, 19.0, 4);
    REQUIRE_CLOSE(src_circle.facets[16][0][0], 20.0, 1e-4);
    REQUIRE_CLOSE(src_circle.facets[32][0][0], 1.0, 1e-4);
    REQUIRE_CLOSE(src_circle.facets[48][0][0], 20.0, 1e-4);
    REQUIRE_CLOSE(src_circle.facets[0][0][0], 39.0, 1e-4);
    REQUIRE_CLOSE(src_circle.facets[16][0][1], 19.0, 1e-4);
    REQUIRE_CLOSE(src_circle.facets[32][0][1], 0.0, 1e-4);
    REQUIRE_CLOSE(src_circle.facets[48][0][1], -19.0, 1e-4);
    REQUIRE_CLOSE(src_circle.facets[0][0][1], 0.0, 1e-4);
}

double surface_area(Mesh<3> m) 
{
    double sa = 0;
    for (auto f: m.facets) {
        sa += tri_area(f);
    }
    return sa;
}

TEST_CASE("Mesh3D", "[mesh]") 
{
    auto m = sphere_mesh({0,0,0}, 1.0, 4); 
    double sa = surface_area(m);
    REQUIRE_CLOSE(sa, 4 * M_PI, 1e-1);
}

TEST_CASE("refine out of order elements", "[mesh]")
{
    auto m2 = line_mesh({0,0}, {1,0}).refine_repeatedly(4).refine({0, 10, 3, 4, 2});
    REQUIRE(m2.n_facets() == 21);
}
