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
class Mesh {
public:
    std::vector<std::array<double,3>> vertices;
    std::vector<std::array<int,3>> faces;
    bool has_refine_mod;
    std::function<std::array<double,3>(std::array<double,3>)> refine_mod;
};

inline bool same_vertex(std::array<double, 3> a,
                        std::array<double, 3> b,
                        double eps) {
    return std::fabs(a[0] - b[0]) < eps && 
           std::fabs(a[1] - b[1]) < eps && 
           std::fabs(a[2] - b[2]) < eps;
}
inline std::array<double,3> midpt(std::array<double,3> a, std::array<double,3> b) {
    return {(a[0] + b[0]) / 2.0,
            (a[1] + b[1]) / 2.0,
            (a[2] + b[2]) / 2.0};
}

std::unordered_map<int, int> find_duplicate_map(const Mesh& mesh, double eps);

int count_unique_vertices(std::unordered_map<int,int>& old_to_new);

std::vector<std::array<double,3>> unique_vertices(
        const Mesh& mesh, std::unordered_map<int,int>& old_to_new);

Mesh clean_mesh(const Mesh& mesh, double vertex_smear = 1e-6);


/* Produces 4 new triangles from an initial triangle by adding a
 * vertex to each edge of the triangle
 *         /\          /\
 *        /  \   ->   /__\
 *       /    \      /\  /\
 *      /______\    /__\/__\
 *      (My first ever ASCII art, a masterpiece that will
 *      stand the test of time!)
 */
void refine_face(Mesh& new_mesh, std::array<int, 3> face);

Mesh refine_mesh(const Mesh& m, std::vector<int> refine_these);
Mesh refine_mesh(const Mesh& m);

/* MESH GENERATION FUNCTION */

//TODO :Eventually move these to another file.

Mesh cube_mesh();


Mesh sphere_mesh(std::array<double,3> center, double radius, bool interior = true);

#endif
