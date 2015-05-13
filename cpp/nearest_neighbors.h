#ifndef __QWEKLHHHHAHAHHA_NEAREST_NEIGHBORS_H
#define __QWEKLHHHHAHAHHA_NEAREST_NEIGHBORS_H

#include <cassert>
#include <limits>
#include "octree.h"
#include "vec_ops.h"
#include "geometry.h"
#include "closest_pt.h"
#include "numerics.h"

namespace tbem {

template <size_t dim>
std::vector<std::pair<size_t,size_t>> nearby_points_all_pairs(
    const std::vector<Vec<double,dim>>& ptsA,
    const std::vector<Vec<double,dim>>& ptsB,
    const Octree<dim>& cellA,
    const Octree<dim>& cellB, 
    double r)
{
    auto rA = hypot(cellA.data.bounds.half_width);
    auto rB = hypot(cellB.data.bounds.half_width);
    auto sep = hypot(cellB.data.bounds.center - cellA.data.bounds.center);
    if (sep > rA + rB + 2 * r) {
        // the cells are too far apart, ignore them
        return {};
    }

    if (cellA.is_leaf() && cellB.is_leaf()) {
        // the cells are intersecting leaf cells. perform a direct intersection
        std::vector<std::pair<size_t,size_t>> out_pairs;
        for (const auto& idxA: cellA.data.indices) {
            for (const auto& idxB: cellB.data.indices) {
                bool close_pts = all(fabs(ptsA[idxA] - ptsB[idxB]) <= r);
                if (close_pts) {
                    out_pairs.push_back({idxA, idxB});
                }
            }
        }
        return out_pairs;
    }

    // recurse
    std::vector<std::pair<size_t,size_t>> out_pairs;
    bool B_is_shallower = cellA.data.level > cellB.data.level;
    bool splitB = (B_is_shallower && !cellB.is_leaf()) || cellA.is_leaf();
    if (splitB) {
        //split B because it is shallower
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (cellB.children[c] == nullptr) {
                continue;
            }
            auto child_pairs = nearby_points_all_pairs(
                ptsA, ptsB, cellA, *cellB.children[c], r
            );
            out_pairs.insert(out_pairs.end(), child_pairs.begin(), child_pairs.end());
        }
    } else {
        //split A because it is shallower
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (cellA.children[c] == nullptr) {
                continue;
            }
            auto child_pairs = nearby_points_all_pairs(
                ptsA, ptsB, *cellA.children[c], cellB, r
            );
            out_pairs.insert(out_pairs.end(), child_pairs.begin(), child_pairs.end());
        }
    }

    return out_pairs;
}

template <size_t dim>
std::vector<size_t>
nearby_points(const Vec<double,dim>& ptA, 
    const std::vector<Vec<double,dim>>& ptsB,
    const Octree<dim>& octB, double r)
{
    auto octA = build_octree<dim>({ptA}, 1);
    auto pairs = nearby_points_all_pairs({ptA}, ptsB, octA, octB, r);
    std::vector<size_t> out_indices(pairs.size());
    for (size_t i = 0; i < pairs.size(); i++) {
        assert(pairs[i].first == 0);
        out_indices[i] = pairs[i].second;
    }
    return out_indices;
}

template <size_t dim>
struct NearestFacets {
    const std::vector<Vec<Vec<double,dim>,dim>> facets;
    const Vec<double,dim> pt;
    const double distance;
};

//TODO: Use the bounding balls to speed this up as still as O(n^2) algorithm
//TODO: Do this using an octree to get faster asymptotic times
template <size_t dim>
NearestFacets<dim> nearest_facets(const Vec<double,dim>& pt,
    const std::vector<Vec<Vec<double,dim>,dim>>& facets) 
{
    std::vector<Vec<Vec<double,dim>,dim>> closest_facets;
    Vec<double,dim> closest_pt;
    auto min_dist2 = std::numeric_limits<double>::max();
    for (size_t facet_idx = 0; facet_idx < facets.size(); facet_idx++) {
        auto ref_pt = closest_pt_facet(pt, facets[facet_idx]);
        auto mesh_pt = ref_to_real(ref_pt, facets[facet_idx]);
        auto dist2_to_mesh = dist2(pt, mesh_pt);
        if (dist2_to_mesh < min_dist2) {
            closest_facets.clear();
            closest_facets.push_back(facets[facet_idx]);
            min_dist2 = dist2_to_mesh;
            closest_pt = mesh_pt;
        } else if (dist2_to_mesh == min_dist2) {
            closest_facets.push_back(facets[facet_idx]);
        }
    }
    return {closest_facets, closest_pt, std::sqrt(min_dist2)};
}

template <size_t dim>
Vec<double,dim> away_from_facet_dir(const Vec<Vec<double,dim>,dim>& f,
    const Vec<double,dim>& pt) 
{
    auto facet_normal = unscaled_normal(f); 
    auto which_side = which_side_point(f, pt);
    if (which_side == INTERSECT) {
        return facet_normal;
    } else if (which_side == FRONT) {
        return facet_normal;
    } else {
        return -facet_normal;
    }
}

template <size_t dim>
Vec<double,dim> decide_richardson_dir(const Vec<double,dim>& pt,
    const NearestFacets<dim>& nearest_facets) 
{
    auto direction = zeros<Vec<double,dim>>::make();
    for (size_t i = 0; i < nearest_facets.facets.size(); i++) {
        direction += away_from_facet_dir(nearest_facets.facets[i], pt);
    }
    if (all(direction == 0.0)) {
        direction = away_from_facet_dir(nearest_facets.facets[0], pt);
    }
    direction /= (double)nearest_facets.facets.size();
    return direction;
}

} //end namespace tbem

#endif
