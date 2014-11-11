#ifndef __MESH_3D_H
#define __MESH_3D_H
#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <GL/gl.h>
#include "vec.h"

template <int dim>
class Mesh {
public:
    std::vector<Vec<double,dim>> vertices;
    std::vector<std::array<int,dim>> faces;
    bool has_refine_mod;
    std::function<Vec<double,dim>(Vec<double,dim>)> refine_mod;
};

typedef Mesh<3> Mesh3D;
typedef Mesh<2> Mesh2D;

inline bool same_vertex(const Vec3<double>& a,
                        const Vec3<double>& b,
                        double eps) {
    return std::fabs(a[0] - b[0]) < eps && 
           std::fabs(a[1] - b[1]) < eps && 
           std::fabs(a[2] - b[2]) < eps;
}

std::unordered_map<int, int> find_duplicate_map(const Mesh3D& mesh, double eps);

int count_unique_vertices(const std::unordered_map<int,int>& old_to_new);

std::vector<Vec3<double>> unique_vertices(
        const Mesh3D& mesh, std::unordered_map<int,int>& old_to_new);

Mesh3D clean_mesh(const Mesh3D& mesh, double vertex_smear = 1e-6);


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

/* Refine each triangle specified by the refine_these list into 4 
 * subtriangles.
 */
Mesh3D refine_mesh(const Mesh3D& m, std::vector<int> refine_these);

/* Refines all the triangles.
 */
Mesh3D refine_mesh(const Mesh3D& m);

/* Refine a mesh multiple times and then clean it.
 */
Mesh3D refine_clean(const Mesh3D& m, unsigned int times);

template <int dim>
struct Facet {
    const Vec<Vec<double,dim>,dim> vertices;
};

template <int dim>
struct NewMesh {
    typedef std::function<Vec<double,dim>(Vec<double,dim>)> RefineFnc;
    const std::vector<Facet<dim>> facets;
    const bool has_refine_mod;
    const RefineFnc refine_mod;

    NewMesh<dim> refine(const std::vector<int>& refine_these) const;
    NewMesh<dim> refine() const;
    NewMesh<dim> refine_repeatedly(unsigned int times) const;

    static NewMesh<dim> from_vertices_faces(const std::vector<Vec<double,dim>>& vertices,
                         const std::vector<std::array<int,dim>>& facets,
                         bool has_refine_mod = false,
                         const typename NewMesh<dim>::RefineFnc& refine_mod = nullptr);
};


#endif
