#ifndef __MESH_3D_H
#define __MESH_3D_H
#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <GL/gl.h>

//TODO: Much of this code is identical for a 2D mesh (in bem.h), see if
//integrating the two would be worthwhile. Or template it..
//TODO: move this code to a cpp file
//TODO: fix the O(N^2) problem in the find_duplicate_map function
//      possibilities: use a locality-sensitive hashing scheme
//                     build the octree first and use that for search
class Mesh3D {
public:
    std::vector<std::array<double,3>> vertices;
    std::vector<std::array<int,3>> faces;
};

bool same_vertex(std::array<double, 3> a, std::array<double, 3> b, double eps) {
    return abs(a[0] - b[0]) < eps && 
           abs(a[1] - b[1]) < eps && 
           abs(a[2] - b[2]) < eps;
}

std::unordered_map<int, int> find_duplicate_map(Mesh3D& mesh, double eps) {
    std::unordered_map<int, int> old_to_new;
    for (unsigned int i = 0; i < mesh.vertices.size(); i++) {
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
        Mesh3D& mesh, std::unordered_map<int,int>& old_to_new) {
    int n_unique = count_unique_vertices(old_to_new);
    std::vector<std::array<double,3>> vertices(n_unique);
    for (unsigned int i = 0; i < mesh.vertices.size(); i++) {
        vertices[old_to_new[i]] = mesh.vertices[i]; 
    }
    return vertices;
}
            
Mesh3D clean_mesh(Mesh3D& mesh, double vertex_smear = 1e-5) {
    //Find duplicate vertices
    auto old_to_new = find_duplicate_map(mesh, vertex_smear);
    auto new_vertices = unique_vertices(mesh, old_to_new);

    std::vector<std::array<int,3>> new_faces;
    for (unsigned int i = 0; i < mesh.faces.size(); i++) {
        std::array<int,3> face;
        for (int v = 0; v < 3; v++) {
            face[v] = old_to_new[mesh.faces[i][v]];
        }
        new_faces.push_back(face);
    }

    Mesh3D retval = {new_vertices, new_faces};
    return retval;
}

std::array<double,3> midpt(std::array<double,3> a, std::array<double,3> b) {
    return {(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2])/ 2.0};
}

/* Produces 4 new triangles from an initial triangle by adding a
 * vertex to each edge of the triangle
 *         /\          /\
 *        /  \   ->   /__\
 *       /    \      /\  /\
 *      /______\    /__\/__\
 *      (My first ever ASCII art, a masterpiece that will
 *      stand the test of time!)
 */
void refine_face(Mesh3D& new_mesh, std::array<int, 3> face) {
    // Find the new vertex and faces.
    const auto v0 = new_mesh.vertices[face[0]];
    const auto v1 = new_mesh.vertices[face[1]];
    const auto v2 = new_mesh.vertices[face[2]];

    // Calculate the midpoints of each edge of the triangle. These
    // are used to refine the triangle
    const auto midpt01 = midpt(v0, v1);
    const auto midpt12 = midpt(v1, v2);
    const auto midpt20 = midpt(v2, v0);

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

Mesh3D refine_mesh(Mesh3D& m, std::vector<int> refine_these) {
    Mesh3D new_mesh;
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

void draw_mesh(Mesh3D& msh) {
    glBegin(GL_TRIANGLES);
    for (unsigned int i = 0; i < msh.faces.size(); i++) {
        for (int v = 0; v < 3; v++) {
            int vert = msh.faces[i][v];
            glVertex3f(msh.vertices[vert][0],
                       msh.vertices[vert][1],
                       msh.vertices[vert][2]);
        }
    }
    glEnd();
}

// chunks to write:
// 3D: 
// mesh cleaning (DONE)
//
// 2D: 
// mesh cleaning 
// subregion determination 
// boundary conditions (PARTIAL)
// refine func (DONE)
// summation func (DONE)
// richardson extrapolation quadrature (NEEDSTEST)
// mappings from reference to real space (DONE)
// basis (DONE)
// constraints --> boundary conditions, non-singular traction BCs
// which kernels to use for the different boundary integral equations
// the kernels
// evaluate solution on the surface after calculation -- just use the 
//     integral equation again
// interior meshing and evaluation
// adaptivity?
// that old list of good things to do!
// transcribe all those old sheets of thoughts
// look into UFL from Fenics as the 

#endif
