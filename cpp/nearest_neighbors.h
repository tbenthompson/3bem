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
    const std::vector<Vec<Vec<double,dim>,dim>>& facets)
{
    double dist2_to_closest_facet = std::numeric_limits<double>::max();
    size_t closest_facet_idx = 0;
    auto nearest_pt = pt;
    for (size_t facet_idx = 0; facet_idx < facets.size(); facet_idx++) {
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
size_t find_closest_child(const Octree<dim>& oct, const Vec<double,dim>& pt) 
{
    size_t closest_child = 0;
    double dist2_to_closest_child = std::numeric_limits<double>::max();
    for (size_t c_idx = 0; c_idx < Octree<dim>::split; c_idx++) {
        auto c = oct.children[c];
        if (c == nullptr) {
            continue;
        }
        auto center = c->data.bounds.center;
        auto dist2_child = dist2(pt, center);
        if (dist2_child < dist2_to_closest_child) {
            closest_child = c_idx;
            dist2_to_closest_child = dist2_child;
        }
    }
    return closest_child;
}

// 
// //find_children
// 
// template <size_t dim>
// size_t nearest_facet_helper(const Octree<dim>& oct,
//     const Vec<double,dim>& pt, const std::vector<Vec<Vec<double,dim>,dim>>& facets)
// {
//     if (oct.is_leaf()) {
//         // found the bottom, now find the closest facet and then go back up!
//         return nearest_facet_brute_force(oct, pt, facets);
//     } else {
//         auto closest_child = find_closest_child(oct, pt);
//         auto current_closest = nearest_facet_helper(
//             *oct.children[closest_child], pt, facets
//         );
// 
//     }
// }
// 
// template <size_t dim>
// NearestFacets<dim> nearest_facet(const Vec<double,dim>& pt,
//     const std::vector<Vec<Vec<double,dim>,dim>>& facets) 
// {
//     assert(facets.size() > 0);
//     std::vector<Ball<dim>> balls(facets.size());
//     for (size_t i = 0; i < facets.size(); i++) {
//         balls[i] = facet_ball(facets[i]);
//     }
//     auto oct = make_octree<dim>(balls, 1);
//     return {{}, pt, 0};
// }

// template <size_t dim>
// NearestFacets<dim> nearest_facets(const Vec<double,dim>& pt,
//     const std::vector<Vec<Vec<double,dim>,dim>>& facets) 
// {
//     assert(facets.size() > 0);
//     std::vector<Vec<Vec<double,dim>,dim>> closest_facets;
//     Vec<double,dim> closest_pt;
//     auto min_dist2 = std::numeric_limits<double>::max();
//     for (size_t facet_idx = 0; facet_idx < facets.size(); facet_idx++) {
//         auto ref_pt = closest_pt_facet(pt, facets[facet_idx]);
//         auto mesh_pt = ref_to_real(ref_pt, facets[facet_idx]);
//         auto dist2_to_mesh = dist2(pt, mesh_pt);
//         const double nearby_factor = 0.2;
//         if (dist2_to_mesh == min_dist2 || 
//             std::fabs(dist2_to_mesh - min_dist2) < nearby_factor * min_dist2) 
//         {
//             closest_facets.push_back(facets[facet_idx]);
//         } else if (dist2_to_mesh < min_dist2) {
//             closest_facets.clear();
//             closest_facets.push_back(facets[facet_idx]);
//             min_dist2 = dist2_to_mesh;
//             closest_pt = mesh_pt;
//         } 
//     }
//     return {closest_facets, closest_pt, std::sqrt(min_dist2)};
// }

} //end namespace tbem

#endif
