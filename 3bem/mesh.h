#ifndef __MESH_3D_H
#define __MESH_3D_H
#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <GL/gl.h>
#include "vec.h"

namespace tbem {

/* There's a pretty interesting correspondence between the interpolation
 * of a field like position between the topological corners of a triangle
 * and the interpolation of a field like displacment, traction, or heat flux
 * between the same corners. I could generalize the mesh structure into 
 * something like a InterpolatedField structure.
 */

template <typename T, int dim>
struct FacetField {
    const Vec<T,dim> vertices;
};

template <int dim>
using Facet = FacetField<Vec<double,dim>,dim>;

template <typename T, int dim>
struct MeshField {
    typedef std::function<T(T)> RefineFnc;
    const std::vector<FacetField<T,dim>> facets;
    const bool has_refine_mod;
    const RefineFnc refine_mod;

    MeshField<T,dim> refine(const std::vector<int>& refine_these) const;
    MeshField<T,dim> refine() const;
    MeshField<T,dim> refine_repeatedly(unsigned int times) const;
};

template <int dim>
using Mesh = MeshField<Vec<double,dim>,dim>;

template <int dim>
Mesh<dim> mesh_from_vertices_faces(const std::vector<Vec<double,dim>>& vertices,
                         const std::vector<std::array<int,dim>>& facets,
                         bool has_refine_mod = false,
                         const typename Mesh<dim>::RefineFnc& refine_mod = nullptr);

} //END NAMESPACE tbem

#endif
