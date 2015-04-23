#ifndef __QWEKLHHHHAHAHHA_NEAREST_NEIGHBORS_H
#define __QWEKLHHHHAHAHHA_NEAREST_NEIGHBORS_H

#include <cassert>
#include "octree.h"
#include "vec_ops.h"

namespace tbem {

template <size_t dim>
std::vector<std::pair<size_t,size_t>> identical_points_all_pairs(
    const std::vector<Vec<double,dim>>& ptsA,
    const std::vector<Vec<double,dim>>& ptsB,
    const Octree<dim>& cellA,
    const Octree<dim>& cellB)
{
    auto rA = hypot(cellA.data.bounds.half_width);
    auto rB = hypot(cellB.data.bounds.half_width);
    auto sep = hypot(cellB.data.bounds.center - cellA.data.bounds.center);
    if (sep > rA + rB) {
        // the cells don't intersect. ignore them!
        return {};
    }

    if (cellA.is_leaf() && cellB.is_leaf()) {
        // the cells are intersecting leaf cells. perform a direct intersection
        std::vector<std::pair<size_t,size_t>> out_pairs;
        for (const auto& idxA: cellA.data.indices) {
            for (const auto& idxB: cellB.data.indices) {
                if (ptsA[idxA] == ptsB[idxB]) {
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
            auto child_pairs = identical_points_all_pairs(
                ptsA, ptsB, cellA, *cellB.children[c]
            );
            out_pairs.insert(out_pairs.end(), child_pairs.begin(), child_pairs.end());
        }
    } else {
        //split A because it is shallower
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (cellA.children[c] == nullptr) {
                continue;
            }
            auto child_pairs = identical_points_all_pairs(
                ptsA, ptsB, *cellA.children[c], cellB
            );
            out_pairs.insert(out_pairs.end(), child_pairs.begin(), child_pairs.end());
        }
    }

    return out_pairs;
}

template <size_t dim>
std::vector<size_t>
identical_points(const Vec<double,dim>& ptA, 
    const std::vector<Vec<double,dim>>& ptsB,
    const Octree<dim>& octB)
{
    auto octA = build_octree<dim>({ptA}, 1);
    auto pairs = identical_points_all_pairs({ptA}, ptsB, octA, octB);
    std::vector<size_t> out_indices(pairs.size());
    for (size_t i = 0; i < pairs.size(); i++) {
        assert(pairs[i].first == 0);
        out_indices[i] = pairs[i].second;
    }
    return out_indices;
}



} //end namespace tbem

#endif
