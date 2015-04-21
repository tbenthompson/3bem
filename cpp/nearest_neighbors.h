#ifndef __QWEKLHHHHAHAHHA_NEAREST_NEIGHBORS_H
#define __QWEKLHHHHAHAHHA_NEAREST_NEIGHBORS_H

#include "octree.h"

namespace tbem {

template <size_t dim>
std::vector<size_t>
identical_points(const Vec<double,dim>& pt, 
    const std::vector<Vec<double,dim>>& target_pts,
    const Octree<dim>& oct)
{
    if (!oct.data.bounds.in_box_inclusive(pt)) {
        return {};
    }
    if (oct.is_leaf()) {
        std::vector<size_t> out_indices;
        for (const auto& t_idx: oct.data.indices) {
            if (pt == target_pts[t_idx]) {
                out_indices.push_back(t_idx);
            }
        }
        return out_indices;
    }
    std::vector<size_t> out_indices;
    for (const auto& c: oct.children) {
        if (c == nullptr) {
            continue;
        }
        auto c_indices = identical_points(pt, target_pts, *c);
        out_indices.insert(out_indices.end(), c_indices.begin(), c_indices.end()); 
    }
    return out_indices;
}

} //end namespace tbem

#endif
