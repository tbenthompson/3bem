#ifndef __GGggGGggTTFDSSDf_MESH_H
#define __GGggGGggTTFDSSDf_MESH_H
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
 * something like a InterpolatedField structure (or MeshField!).
 */

template <typename T, size_t dim>
struct FacetField {
    const Vec<T,dim> vertices;
};

template <size_t dim>
using Facet = FacetField<Vec<double,dim>,dim>;

template <typename T, size_t dim>
struct FacetCornerIterator;

template <typename T, size_t dim>
struct MeshField {
    typedef std::function<T(T)> RefineFnc;
    const std::vector<FacetField<T,dim>> facets;

    FacetCornerIterator<T,dim> begin() const;
    FacetCornerIterator<T,dim> end() const;

    MeshField<T,dim> refine(const std::vector<int>& refine_these) const;
    MeshField<T,dim> refine() const;
    MeshField<T,dim> refine_repeatedly(unsigned int times) const;

    static MeshField<T,dim> form_union(const std::vector<MeshField<T,dim>>& others);

    static MeshField<T,dim> from_vertices_faces(const std::vector<T>& vertices,
        const std::vector<std::array<int,dim>>& facets_by_vert_idx);
};

template <size_t dim>
using Mesh = MeshField<Vec<double,dim>,dim>;

} //END NAMESPACE tbem

#endif
