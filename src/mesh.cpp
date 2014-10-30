#include "mesh.h"
#include "util.h"
#include "vec.h"

std::unordered_map<int, int> find_duplicate_map(const Mesh& mesh, double eps) {
    std::unordered_map<int, int> old_to_new;
    for (unsigned int i = 0; i < mesh.vertices.size(); i++) {
        //TODO: FIX THE O(N^2) problem
        for (unsigned int j = i + 1; j < mesh.vertices.size(); j++) {
            if (!same_vertex(mesh.vertices[i], mesh.vertices[j], eps)) {
                continue;
            } 
            // Map the deleted vertex to the label of the retained vertex.
            // The retained vertex will always be the last in the set of 
            // identical vertices. So, if indices 1,3,9,11 are all the same
            // point, only vertex 11 will be kept.
            old_to_new[i] = j;
        }
        auto found = old_to_new.find(i);
        if (found == old_to_new.end()) {
            old_to_new[i] = i;
        }
    }

    int next_new_vertex = 0;
    // Traverse the vertices in reverse, so that when we get to a deleted
    // vertex it has already been assigned its new location.
    for (int i = (int)mesh.vertices.size() - 1; i >= 0; i--) {
        if (old_to_new[i] != (int)i) {
            old_to_new[i] = old_to_new[old_to_new[i]];
        } else {
            old_to_new[i] = next_new_vertex;
            next_new_vertex++;
        }
    }

    return old_to_new;
}

int count_unique_vertices(std::unordered_map<int,int>& old_to_new) {
    int max_index = 0;
    for (auto it = old_to_new.begin(); it != old_to_new.end(); it++) {
        max_index = std::max(max_index, it->second);
    }
    return max_index + 1;
}

std::vector<std::array<double,3>> unique_vertices(
        const Mesh& mesh, std::unordered_map<int,int>& old_to_new) {
    int n_unique = count_unique_vertices(old_to_new);
    std::vector<std::array<double,3>> vertices(n_unique);
    for (unsigned int i = 0; i < mesh.vertices.size(); i++) {
        vertices[old_to_new[i]] = mesh.vertices[i]; 
    }
    return vertices;
}
            
Mesh clean_mesh(const Mesh& mesh, double vertex_smear) {
    //Find duplicate vertices
    auto old_to_new = find_duplicate_map(mesh, vertex_smear);
    auto new_vertices = unique_vertices(mesh, old_to_new);
    // std::cout << mesh.vertices.size() << " " << new_vertices.size() << " " << 
        // count_unique_vertices(old_to_new) << std::endl;

    std::vector<std::array<int,3>> new_faces;
    for (unsigned int i = 0; i < mesh.faces.size(); i++) {
        std::array<int,3> face;
        for (int v = 0; v < 3; v++) {
            face[v] = old_to_new[mesh.faces[i][v]];
        }
        new_faces.push_back(face);
    }

    Mesh retval{new_vertices, new_faces, mesh.has_refine_mod, mesh.refine_mod};
    return retval;
}

void refine_face(Mesh& new_mesh, std::array<int, 3> face) {
    // Find the new vertex and faces.
    const auto v0 = new_mesh.vertices[face[0]];
    const auto v1 = new_mesh.vertices[face[1]];
    const auto v2 = new_mesh.vertices[face[2]];

    // Calculate the midpoints of each edge of the triangle. These
    // are used to refine the triangle
    auto midpt01 = midpt(v0, v1);
    auto midpt12 = midpt(v1, v2);
    auto midpt20 = midpt(v2, v0);
    if (new_mesh.has_refine_mod) {
        midpt01 = new_mesh.refine_mod(midpt01);
        midpt12 = new_mesh.refine_mod(midpt12);
        midpt20 = new_mesh.refine_mod(midpt20);
    }

    // Add the vertices while grabbing their indices.
    int midpt01_idx = new_mesh.vertices.size();
    new_mesh.vertices.push_back(midpt01);
    int midpt12_idx = new_mesh.vertices.size();
    new_mesh.vertices.push_back(midpt12);
    int midpt20_idx = new_mesh.vertices.size();
    new_mesh.vertices.push_back(midpt20);

    //Maintain the orientation. Since vertex 1 is "next" after vertex 0 
    //in the original triangle, midpt01 should be "next" after vertex 0
    //in the refined triangle. Following this principle for all the triangles
    //gives this set of new faces
    new_mesh.faces.push_back({face[0], midpt01_idx, midpt20_idx});
    new_mesh.faces.push_back({face[1], midpt12_idx, midpt01_idx});
    new_mesh.faces.push_back({face[2], midpt20_idx, midpt12_idx});
    new_mesh.faces.push_back({midpt01_idx, midpt12_idx, midpt20_idx});
}

Mesh refine_mesh(const Mesh& m, std::vector<int> refine_these) {
    Mesh new_mesh {m.vertices, {}, m.has_refine_mod, m.refine_mod};
    new_mesh.vertices = m.vertices;

    // Sort the refined edges so that we only have to check the
    // next one at any point in the loop.
    std::sort(refine_these.begin(), refine_these.end());

    // The next refined edge.
    int current = 0;

    for (int i = 0; i < (int)m.faces.size(); i++) {
        if (i == refine_these[current]) {
            current += 1;
            refine_face(new_mesh, m.faces[i]);
        } else {
            new_mesh.faces.push_back(m.faces[i]);
        }
    }

    return new_mesh;
}

Mesh refine_mesh(const Mesh& m) {
    return refine_mesh(m, naturals(m.faces.size()));
}

Mesh cube_mesh() {
    std::vector<std::array<double, 3>> vertices = {
        {0.0,0.0,0.0}, {0.0,0.0,1.0}, {0.0,1.0,1.0},
        {1.0,1.0,0.0}, {0.0,0.0,0.0}, {0.0,1.0,0.0},
        {1.0,0.0,1.0}, {0.0,0.0,0.0}, {1.0,0.0,0.0},
        {1.0,1.0,0.0}, {1.0,0.0,0.0}, {0.0,0.0,0.0},
        {0.0,0.0,0.0}, {0.0,1.0,1.0}, {0.0,1.0,0.0},
        {1.0,0.0,1.0}, {0.0,0.0,1.0}, {0.0,0.0,0.0},
        {0.0,1.0,1.0}, {0.0,0.0,1.0}, {1.0,0.0,1.0},
        {1.0,1.0,1.0}, {1.0,0.0,0.0}, {1.0,1.0,0.0},
        {1.0,0.0,0.0}, {1.0,1.0,1.0}, {1.0,0.0,1.0},
        {1.0,1.0,1.0}, {1.0,1.0,0.0}, {0.0,1.0,0.0},
        {1.0,1.0,1.0}, {0.0,1.0,0.0}, {0.0,1.0,1.0},
        {1.0,1.0,1.0}, {0.0,1.0,1.0}, {1.0,0.0,1.0}
    };

    std::vector<std::array<int, 3>> faces;
    for (int i = 0; i < 12; i++) {
        faces.push_back({3 * i, 3 * i + 1, 3 * i + 2});
    }

    Mesh cube{vertices, faces, false, nullptr};
    cube = clean_mesh(cube);

    return cube;
}

Mesh sphere_mesh(std::array<double,3> center, double r, bool interior) {
    std::vector<std::array<double,3>> vertices =
    {
        {0.0, -r, 0.0}, {r, 0.0, 0.0}, {0.0, 0.0, r},
        {-r, 0.0, 0.0}, {0.0, 0.0, -r}, {0.0, r, 0.0}
    };

    for (unsigned int i = 0; i < vertices.size(); i++) {
        for (int d = 0; d < 3; d++) {
            vertices[i][d] += center[d];
        }
    }

    std::vector<std::array<int,3>> faces = 
    {
        {1, 0, 2}, {2, 0, 3}, {3, 0, 4}, {4, 0, 1},
        {5, 1, 2}, {5, 2, 3}, {5, 3, 4}, {5, 4, 1}
    };
    if (!interior) {
        for (unsigned int f = 0; f < faces.size(); f++) {
            std::swap(faces[f][0], faces[f][1]);
        }
    }

    Mesh octahedron{vertices, faces, true, 
        [=](std::array<double,3> x) {
            double dist = std::sqrt(dist2(x, center));
            x[0] = (r / dist) * (x[0] - center[0]) + center[0];
            x[1] = (r / dist) * (x[1] - center[1]) + center[1];
            x[2] = (r / dist) * (x[2] - center[2]) + center[2];
            return x;
        }
    };

    return octahedron;
}
