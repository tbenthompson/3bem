#include "UnitTest++.h"
#include "util.h"
#include "numerics.h"
#include "mesh.h"
#include "mesh_gen.h"

Mesh3D test_mesh() {
    Mesh3D mesh;
    mesh.vertices.push_back({0.0, 1.0, 1.0});
    mesh.vertices.push_back({0.0, 0.0, 1.0});
    mesh.vertices.push_back({0.0, 0.0, 1.0 - 1e-8});
    mesh.vertices.push_back({0.0, 1.0 + 1e-8, 1.0 - 1e-8});
    mesh.vertices.push_back({1.0, 0.0, 0.0});
    mesh.faces.push_back({2, 3, 4});
    return mesh;
}

TEST(DuplicateMap) {
    auto mesh = test_mesh();
    auto map = find_duplicate_map(mesh, 1e-5);
    CHECK_CLOSE(map[0], 1, 1e-5);
    CHECK_CLOSE(map[1], 2, 1e-5);
    CHECK_CLOSE(map[2], 2, 1e-5);
    CHECK_CLOSE(map[3], 1, 1e-5);
    CHECK_CLOSE(map[4], 0, 1e-5);
}

TEST(CountUnique) {
    auto mesh = test_mesh();
    auto map = find_duplicate_map(mesh, 1e-5);
    int unique = count_unique_vertices(map);
    CHECK_CLOSE(unique, 3, 1e-5);
}

TEST(UniqueVertices) {
    auto mesh = test_mesh();
    auto map = find_duplicate_map(mesh, 1e-5);
    auto verts = unique_vertices(mesh, map);
    double correct[3][3] = {{1.0, 0.0, 0.0},
                            {0.0, 1.0 + 1e-8, 1.0 - 1e-8},
                            {0.0, 0.0, 1.0 - 1e-8}};
    CHECK_ARRAY_CLOSE(verts[0], correct[0], 3, 1e-5);
    CHECK_ARRAY_CLOSE(verts[1], correct[1], 3, 1e-5);
    CHECK_ARRAY_CLOSE(verts[2], correct[2], 3, 1e-5);
}

TEST(CleanMesh) {
    auto mesh = test_mesh();
    auto cleaned = clean_mesh(mesh);
    int correct[3] = {2, 1, 0};
    CHECK_ARRAY_CLOSE(cleaned.faces[0], correct, 3, 1e-5);
}

TEST(RefineSphereMesh) {
    auto unrefined = sphere_mesh({0,0,0}, 1.0);
    //TODO: This could be tested easily for general meshes using some kind of
    //property based test
    auto refined = refine_mesh(unrefined, naturals(unrefined.faces.size()));
    int correct = unrefined.faces.size() * 3 + unrefined.vertices.size();
    int correct_faces = unrefined.faces.size() * 4;
    CHECK_EQUAL(refined.vertices.size(), correct);
    CHECK_EQUAL(refined.faces.size(), correct_faces);

    auto cleaned = clean_mesh(refined);
    int correct_clean = 18;
    CHECK_EQUAL(cleaned.vertices.size(), correct_clean);
    CHECK_EQUAL(cleaned.faces.size(), correct_faces);
    for (unsigned int i = 0; i < refined.faces.size(); i++) {
        for (int d = 0; d < 3; d++) {
            auto v = refined.vertices[refined.faces[i][d]];
            double dist = hypot2(v);
            CHECK_CLOSE(dist, 1.0, 1e-8);
        }
    }
}

double perimeter(NewMesh<2> m) {
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

double surface_area(NewMesh<3> m) {
    double sa = 0;
    for (auto f: m.facets) {
        sa += tri_area(f.vertices);
    }
    return sa;
}

TEST(Mesh3D) {
    auto m = sphere_mesh_new({0,0,0}, 1.0); 
    auto m2 = m.refine_repeatedly(4);
    double sa = surface_area(m2);
    CHECK_CLOSE(sa, 4 * M_PI, 1e-1);
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
