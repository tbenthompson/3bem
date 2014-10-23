#include "UnitTest++.h"
#include "util.h"
#include "numerics.h"
#include "mesh.h"

Mesh test_mesh() {
    Mesh mesh;
    mesh.vertices.push_back({0.0, 1.0, 1.0});
    mesh.vertices.push_back({0.0, 0.0, 1.0});
    mesh.vertices.push_back({0.0, 0.0, 1.0 - 1e-8});
    mesh.vertices.push_back({0.0, 1.0 + 1e-8, 1.0 - 1e-8});
    mesh.vertices.push_back({1.0, 0.0, 0.0});
    mesh.faces.push_back({2, 3, 4});
    return mesh;
}

TEST(DuplicateMap) {
    Mesh mesh = test_mesh();
    auto map = find_duplicate_map(mesh, 1e-5);
    CHECK_CLOSE(map[0], 1, 1e-5);
    CHECK_CLOSE(map[1], 2, 1e-5);
    CHECK_CLOSE(map[2], 2, 1e-5);
    CHECK_CLOSE(map[3], 1, 1e-5);
    CHECK_CLOSE(map[4], 0, 1e-5);
}

TEST(CountUnique) {
    Mesh mesh = test_mesh();
    auto map = find_duplicate_map(mesh, 1e-5);
    int unique = count_unique_vertices(map);
    CHECK_CLOSE(unique, 3, 1e-5);
}

TEST(UniqueVertices) {
    Mesh mesh = test_mesh();
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
    Mesh mesh = test_mesh();
    Mesh cleaned = clean_mesh(mesh);
    int correct[3] = {2, 1, 0};
    CHECK_ARRAY_CLOSE(cleaned.faces[0], correct, 3, 1e-5);
}

TEST(RefineSphereMesh) {
    Mesh unrefined = sphere_mesh({0,0,0}, 1.0);
    //TODO: This could be tested easily for general meshes using some kind of
    //property based test
    Mesh refined = refine_mesh(unrefined, naturals(unrefined.faces.size()));
    int correct = unrefined.faces.size() * 3 + unrefined.vertices.size();
    int correct_faces = unrefined.faces.size() * 4;
    CHECK_EQUAL(refined.vertices.size(), correct);
    CHECK_EQUAL(refined.faces.size(), correct_faces);

    Mesh cleaned = clean_mesh(refined);
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

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
