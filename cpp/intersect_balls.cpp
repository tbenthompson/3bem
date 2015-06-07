#include <cassert>
#include "intersect_balls.h"
#include "geometry.h"
#include "vec_ops.h"

namespace tbem {

template <size_t dim>
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
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
            auto child_pairs = intersect_balls_all_pairs(
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
            auto child_pairs = intersect_balls_all_pairs(
                ptsA, ptsB, *cellA.children[c], cellB, r
            );
            out_pairs.insert(out_pairs.end(), child_pairs.begin(), child_pairs.end());
        }
    }

    return out_pairs;
}

template 
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
    const std::vector<Vec<double,2>>& ptsA,
    const std::vector<Vec<double,2>>& ptsB,
    const Octree<2>& cellA,
    const Octree<2>& cellB, 
    double r);

template 
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
    const std::vector<Vec<double,3>>& ptsA,
    const std::vector<Vec<double,3>>& ptsB,
    const Octree<3>& cellA,
    const Octree<3>& cellB, 
    double r);

template <size_t dim>
std::vector<size_t> intersect_balls(const Vec<double,dim>& ptA, 
    const std::vector<Vec<double,dim>>& ptsB,
    const Octree<dim>& octB, double r)
{
    auto octA = make_octree<dim>({ptA}, 1);
    auto pairs = intersect_balls_all_pairs({ptA}, ptsB, octA, octB, r);
    std::vector<size_t> out_indices(pairs.size());
    for (size_t i = 0; i < pairs.size(); i++) {
        assert(pairs[i].first == 0);
        out_indices[i] = pairs[i].second;
    }
    return out_indices;
}
template 
std::vector<size_t> intersect_balls(const Vec<double,2>& ptA, 
    const std::vector<Vec<double,2>>& ptsB,
    const Octree<2>& octB, double r);

template 
std::vector<size_t> intersect_balls(const Vec<double,3>& ptA, 
    const std::vector<Vec<double,3>>& ptsB,
    const Octree<3>& octB, double r);

} // end namespace tbem
