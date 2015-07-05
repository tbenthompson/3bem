#ifndef TBEM_MESH_PREPROCESS_H
#define TBEM_MESH_PREPROCESS_H

#include "vec.h"
#include <vector>

namespace tbem {

template <size_t dim>
struct FacetIntersection {
    size_t facet_idx_A;
    size_t facet_idx_B;
    std::vector<Vec<double,dim>> pts;
};


template <size_t dim>
struct MeshPreprocessor {

    /* Find points where the facets of two meshes intersect. This function
     * returns points that are either interior to the facet or on an endpoint.
     */
    std::vector<FacetIntersection<dim>> find_intersections(
        const std::vector<Vec<Vec<double,dim>,dim>>& facetsA,
        const std::vector<Vec<Vec<double,dim>,dim>>& facetsB);

    /* Takes intersections returned by find_intersections and splits the
     * facets of the first mesh and the points where they intersect with
     * the second mesh.
     */
    std::vector<Vec<Vec<double,dim>,dim>> split_facets_at_intersections(
        const std::vector<Vec<Vec<double,dim>,dim>>& facetsA,
        const std::vector<FacetIntersection<dim>>& intersections);
};



} // end namespace tbem

#endif
