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
std::vector<std::pair<size_t,size_t>> nearby_points_all_pairs(
    const std::vector<Vec<double,dim>>& ptsA,
    const std::vector<Vec<double,dim>>& ptsB,
    const Octree<dim>& cellA,
    const Octree<dim>& cellB, 
    double r)
{
    auto rA = hypot(cellA.data.true_bounds.half_width);
    auto rB = hypot(cellB.data.true_bounds.half_width);
    auto sep = hypot(cellB.data.true_bounds.center - 
        cellA.data.true_bounds.center);
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
    auto octA = make_octree<dim>({ptA}, 1);
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
    const Vec<Vec<double,dim>,dim> nearest_facet;
    const Vec<double,dim> pt;
    const double distance;
};

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

template <size_t dim>
size_t nearest_facet_brute_force(const std::vector<size_t>& indices,
    const Vec<double,dim>& pt, const std::vector<Vec<Vec<double,dim>,dim>>& facets)
{
    double dist2_to_closest_facet = 0.0;
    size_t closest_facet_idx = std::numeric_limits<double>::max();
    for (auto facet_idx: indices) {
        auto ref_pt = closest_pt_facet(pt, facets[facet_idx]);
        auto mesh_pt = ref_to_real(ref_pt, facets[facet_idx]);
        auto dist2_to_mesh = dist2(pt, mesh_pt);
        if (dist2_to_mesh < dist2_to_closest_facet) {
            dist2_to_closest_facet = dist2_to_mesh;
            closest_facet_idx = facet_idx;
        }
    }
    return closest_facet_idx;
}

//find_children

template <size_t dim>
size_t nearest_facet_helper(const Octree<dim>& oct,
    const Vec<double,dim>& pt, const std::vector<Vec<Vec<double,dim>,dim>>& facets)
{
    if (oct.is_leaf()) {
        // found the bottom, now find the closest facet and then go back up!
        return nearest_facet_brute_force(oct, pt, facets);
    } else {
        auto closest_child = find_closest_child(oct, pt);
        auto current_closest = nearest_facet_helper(
            *oct.children[closest_child], pt, facets
        );

    }
}

template <size_t dim>
NearestFacets<dim> nearest_facet(const Vec<double,dim>& pt,
    const std::vector<Vec<Vec<double,dim>,dim>>& facets) 
{
    assert(facets.size() > 0);
    std::vector<Ball<dim>> balls(facets.size());
    for (size_t i = 0; i < facets.size(); i++) {
        balls[i] = facet_ball(facets[i]);
    }
    auto oct = make_octree<dim>(balls, 1);
    return {{}, pt, 0};
}

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

template <size_t dim>
NearestFacets<dim> nearest_facets(const Vec<double,dim>& pt,
    const std::vector<Vec<Vec<double,dim>,dim>>& facets) 
{
    assert(facets.size() > 0);
    size_t closest_facet_idx = 0;
    Vec<double,dim> closest_pt;
    auto min_dist2 = std::numeric_limits<double>::max();
    for (size_t facet_idx = 0; facet_idx < facets.size(); facet_idx++) {
        auto ref_pt = closest_pt_facet(pt, facets[facet_idx]);
        auto mesh_pt = ref_to_real(ref_pt, facets[facet_idx]);
        auto dist2_to_mesh = dist2(pt, mesh_pt);
        if (dist2_to_mesh < min_dist2) {
            closest_facet_idx = facet_idx;
            min_dist2 = dist2_to_mesh;
            closest_pt = mesh_pt;
        } 
    }

    auto closest_facet = facets[closest_facet_idx];
    auto search_ball = facet_ball(closest_facet);
    auto search_r2 = std::pow(search_ball.radius, 2);
    if (search_r2 < min_dist2) {
        return {{}, closest_facet, closest_pt, std::sqrt(min_dist2)};
    }

    std::vector<Vec<Vec<double,dim>,dim>> close_facets;
    for (size_t facet_idx = 0; facet_idx < facets.size(); facet_idx++) {
        if (facet_idx == closest_facet_idx) {
            continue;
        }
        auto ref_pt = closest_pt_facet(pt, facets[facet_idx]);
        auto mesh_pt = ref_to_real(ref_pt, facets[facet_idx]);
        auto dist2_to_mesh = dist2(pt, mesh_pt);
        if (dist2_to_mesh <= search_r2) {
            close_facets.push_back(facets[facet_idx]);
        } 
    }

    return {close_facets, closest_facet, closest_pt, std::sqrt(min_dist2)};
}

} //end namespace tbem

#endif
