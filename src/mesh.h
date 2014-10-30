#ifndef __MESH_3D_H
#define __MESH_3D_H
#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <GL/gl.h>
#include "vec.h"

class Mesh {
public:
    std::vector<Vec3<double>> vertices;
    std::vector<std::array<int,3>> faces;
    bool has_refine_mod;
    std::function<Vec3<double>(Vec3<double>)> refine_mod;
};

inline bool same_vertex(const Vec3<double>& a,
                        const Vec3<double>& b,
                        double eps) {
    return std::fabs(a[0] - b[0]) < eps && 
           std::fabs(a[1] - b[1]) < eps && 
           std::fabs(a[2] - b[2]) < eps;
}

std::unordered_map<int, int> find_duplicate_map(const Mesh& mesh, double eps);

int count_unique_vertices(const std::unordered_map<int,int>& old_to_new);

std::vector<Vec3<double>> unique_vertices(
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

/* Refine each triangle specified by the refine_these list into 4 
 * subtriangles.
 */
Mesh refine_mesh(const Mesh& m, std::vector<int> refine_these);

/* Refines all the triangles.
 */
Mesh refine_mesh(const Mesh& m);
#endif
