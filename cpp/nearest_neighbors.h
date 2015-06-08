#ifndef TBEMQWEKLHHHHAHAHHA_NEAREST_NEIGHBORS_H
#define TBEMQWEKLHHHHAHAHHA_NEAREST_NEIGHBORS_H

#include <cassert>
#include <limits>
#include "octree.h"
#include "vec_ops.h"
#include "geometry.h"
#include "gte_wrapper.h"
#include "numerics.h"

namespace tbem {

template <size_t dim>
struct NearestNeighbor {
    size_t idx;
    Vec<double,dim> pt;
    double distance;
};

template <size_t dim>
NearestNeighbor<dim> nearest_facet_brute_force(const Vec<double,dim>& pt,
    const std::vector<Vec<Vec<double,dim>,dim>>& facets,
    const std::vector<size_t>& indices)
{
    double dist2_to_closest_facet = std::numeric_limits<double>::max();
    size_t closest_facet_idx = 0;
    auto nearest_pt = pt;
    for (auto facet_idx: indices) {
        auto ref_pt = closest_pt_facet(pt, facets[facet_idx]);
        auto mesh_pt = ref_to_real(ref_pt, facets[facet_idx]);
        auto dist2_to_mesh = dist2(pt, mesh_pt);
        if (dist2_to_mesh < dist2_to_closest_facet) {
            dist2_to_closest_facet = dist2_to_mesh;
            closest_facet_idx = facet_idx;
            nearest_pt = mesh_pt;
        }
    }
    return {closest_facet_idx, nearest_pt, std::sqrt(dist2_to_closest_facet)};
}

template <size_t dim>
NearestNeighbor<dim> nearest_facet_brute_force(const Vec<double,dim>& pt,
    const std::vector<Vec<Vec<double,dim>,dim>>& facets)
{
    return nearest_facet_brute_force(pt, facets, range(facets.size()));
}


template <size_t dim>
std::vector<Ball<dim>> make_facet_balls(const std::vector<Vec<Vec<double,dim>,dim>>& f)
{
    std::vector<Ball<dim>> out(f.size());
    for (size_t i = 0; i < f.size(); i++) {
        out[i] = facet_ball(f[i]);
    }
    return out;
}

template <size_t dim>
struct NearestNeighborData {
    const static size_t n_facets_per_leaf = 20;
    const std::vector<Vec<Vec<double,dim>,dim>> facets;
    const std::vector<Ball<dim>> facet_balls;
    const Octree<dim> oct;

    NearestNeighborData(const std::vector<Vec<Vec<double,dim>,dim>>& facets):
        facets(facets),
        facet_balls(make_facet_balls(facets)),
        oct(make_octree<dim>(facet_balls, n_facets_per_leaf))
    {
        assert(facets.size() > 0);
    }
};


template <size_t dim>
NearestNeighbor<dim> nearest_facet_helper(const Vec<double,dim>& pt,
    const NearestNeighborData<dim>& nn_data, const Octree<dim>& cell)
{
    if (cell.is_leaf()) {
        return nearest_facet_brute_force(pt, nn_data.facets, cell.indices);
    } else {
        auto closest_child = cell.find_containing_child(pt);
        auto recurse_closest = nearest_facet_helper(
            pt, nn_data, *cell.children[closest_child]
        );
        return recurse_closest;
    }
}

template <size_t dim>
NearestNeighbor<dim> nearest_facet(const Vec<double,dim>& pt,
    const NearestNeighborData<dim>& nn_data)
{
    return nearest_facet_helper(pt, nn_data, nn_data.oct);
}

} //end namespace tbem

#endif
