#include "UnitTest++.h"
#include "test_shared.h"
#include "mesh_3d.h"

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
    Mesh3D mesh = test_mesh();
    auto map = find_duplicate_map(mesh, 1e-5);
    CHECK_CLOSE(map[0], 1, 1e-5);
    CHECK_CLOSE(map[1], 2, 1e-5);
    CHECK_CLOSE(map[2], 2, 1e-5);
    CHECK_CLOSE(map[3], 1, 1e-5);
    CHECK_CLOSE(map[4], 0, 1e-5);
}

TEST(CountUnique) {
    Mesh3D mesh = test_mesh();
    auto map = find_duplicate_map(mesh, 1e-5);
    int unique = count_unique_vertices(map);
    CHECK_CLOSE(unique, 3, 1e-5);
}

TEST(UniqueVertices) {
    Mesh3D mesh = test_mesh();
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
    Mesh3D mesh = test_mesh();
    Mesh3D cleaned = clean_mesh(mesh);
    int correct[3] = {2, 1, 0};
    CHECK_ARRAY_CLOSE(cleaned.faces[0], correct, 3, 1e-5);
}

TEST(RefineMesh) {
    Mesh3D mesh;
    mesh.vertices.push_back({0.0, 0.0, 0.0});
    mesh.vertices.push_back({1.0, 0.0, 0.0});
    mesh.vertices.push_back({0.0, 1.0, 0.0});
    mesh.faces.push_back({{0,1,2}});
    Mesh3D refined = refine_mesh(mesh, {0});
    for (int i = 0; i < refined.faces.size(); i++) {
        for (int d = 0; d < 3; d++) {
            auto v = refined.vertices[refined.faces[i][d]];
            std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
        }
        std::cout << " " << std::endl;
    }
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
