#ifndef __MESH_3D_H
#define __MESH_3D_H
#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>
#include <GL/gl.h>

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
