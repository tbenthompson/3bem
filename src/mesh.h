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

inline bool same_vertex(std::array<double, 3> a,
                        std::array<double, 3> b,
                        double eps) {
    return abs(a[0] - b[0]) < eps && 
           abs(a[1] - b[1]) < eps && 
           abs(a[2] - b[2]) < eps;
}
inline std::array<double,3> midpt(std::array<double,3> a, std::array<double,3> b) {
    return {(a[0] + b[0]) / 2.0,
            (a[1] + b[1]) / 2.0,
            (a[2] + b[2]) / 2.0};
}

std::unordered_map<int, int> find_duplicate_map(Mesh3D& mesh, double eps);

int count_unique_vertices(std::unordered_map<int,int>& old_to_new);

std::vector<std::array<double,3>> unique_vertices(
        Mesh3D& mesh, std::unordered_map<int,int>& old_to_new);

Mesh3D clean_mesh(Mesh3D& mesh, double vertex_smear = 1e-5);


/* Produces 4 new triangles from an initial triangle by adding a
 * vertex to each edge of the triangle
 *         /\          /\
 *        /  \   ->   /__\
 *       /    \      /\  /\
 *      /______\    /__\/__\
 *      (My first ever ASCII art, a masterpiece that will
 *      stand the test of time!)
 */
void refine_face(Mesh3D& new_mesh, std::array<int, 3> face);

Mesh3D refine_mesh(Mesh3D& m, std::vector<int> refine_these);

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
