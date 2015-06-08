#include <cassert>
#include "intersect_balls.h"
#include "geometry.h"
#include "vec_ops.h"
#include "octree.h"

namespace tbem {


template <size_t dim>
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
    const std::vector<Ball<dim>>& ptsA,
    const std::vector<Ball<dim>>& ptsB,
    const Octree<dim>& cellA,
    const Octree<dim>& cellB)
{
    auto rA = hypot(cellA.true_bounds.half_width);
    auto rB = hypot(cellB.true_bounds.half_width);
    auto sep = hypot(cellB.true_bounds.center - 
        cellA.true_bounds.center);
    if (sep > rA + rB) {
        // the cells are too far apart, ignore them
        return {};
    }

    if (cellA.is_leaf() && cellB.is_leaf()) {
        // the cells are intersecting leaf cells. perform a direct intersection
        std::vector<std::pair<size_t,size_t>> out_pairs;
        for (const auto& idxA: cellA.indices) {
            for (const auto& idxB: cellB.indices) {
                if (balls_intersect(ptsA[idxA], ptsB[idxB])) {
                    out_pairs.push_back({idxA, idxB});
                }
            }
        }
        return out_pairs;
    }

    // recurse
    std::vector<std::pair<size_t,size_t>> out_pairs;
    bool B_is_shallower = cellA.level > cellB.level;
    bool splitB = (B_is_shallower && !cellB.is_leaf()) || cellA.is_leaf();
    if (splitB) {
        //split B because it is shallower
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (cellB.children[c] == nullptr) {
                continue;
            }
            auto child_pairs = intersect_balls_all_pairs(
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
            auto child_pairs = intersect_balls_all_pairs(
                ptsA, ptsB, *cellA.children[c], cellB
            );
            out_pairs.insert(out_pairs.end(), child_pairs.begin(), child_pairs.end());
        }
    }

    return out_pairs;
}

template <size_t dim>
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
    const std::vector<Ball<dim>>& ptsA,
    const std::vector<Ball<dim>>& ptsB)
{
    auto octA = make_octree(ptsA, 20);
    auto octB = make_octree(ptsB, 20);
    return intersect_balls_all_pairs(ptsA, ptsB, octA, octB);
}

template 
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
    const std::vector<Ball<2>>& ptsA,
    const std::vector<Ball<2>>& ptsB);

template 
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
    const std::vector<Ball<3>>& ptsA,
    const std::vector<Ball<3>>& ptsB);



template <size_t dim>
std::vector<size_t> intersect_balls(const Ball<dim>& ptA, 
    const std::vector<Ball<dim>>& ptsB, const Octree<dim>& octB)
{
    auto octA = make_octree<dim>({ptA}, 1);
    auto pairs = intersect_balls_all_pairs({ptA}, ptsB, octA, octB);
    std::vector<size_t> out_indices(pairs.size());
    for (size_t i = 0; i < pairs.size(); i++) {
        assert(pairs[i].first == 0);
        out_indices[i] = pairs[i].second;
    }
    return out_indices;
}
template 
std::vector<size_t> intersect_balls(const Ball<2>& ptA, 
    const std::vector<Ball<2>>& ptsB, const Octree<2>& octB);

template 
std::vector<size_t> intersect_balls(const Ball<3>& ptA, 
    const std::vector<Ball<3>>& ptsB, const Octree<3>& octB);

} // end namespace tbem
