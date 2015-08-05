#ifndef ASDJLKJKLQJJQQJQJQJQ_LIMIT_DIRECTION_H
#define ASDJLKJKLQJJQQJQJQJQ_LIMIT_DIRECTION_H

#include <cstdlib>
#include <vector>
#include "vec.h"
#include "nearest_neighbors.h"
#include "facet_info.h"

namespace tbem {

template <size_t dim> struct Ball;

template <size_t dim>
struct NearfieldFacets {
    const std::vector<size_t> facet_indices;
    const std::vector<FacetInfo<dim>> facet_infos;
    const size_t nearest_facet_idx;
    const Vec<double,dim> pt;
    const double distance;
};

template <size_t dim>
struct NearfieldFacetFinder {
    const NearestNeighborData<dim> nn_data;
    double far_threshold;

    NearfieldFacetFinder(const std::vector<Vec<Vec<double,dim>,dim>>& facets,
        double far_threshold);

    NearfieldFacets<dim> find(const Vec<double,dim>& pt) const;

    size_t n_dofs() const 
    {
        return nn_data.facets.size() * dim;
    }
};


// A safety factor of 0.4 should work for all cases.
template <size_t dim>
Vec<double,dim> decide_limit_dir(const Vec<double,dim>& pt,
    const NearfieldFacets<dim>& nearest_facets, 
    const std::vector<Vec<Vec<double,dim>,dim>>& facets, 
    double safety_factor,
    double epsilon = 1e-12);

}

#endif
