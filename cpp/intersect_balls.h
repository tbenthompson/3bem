#ifndef TBEM6781829929292_INTERSECT_BALLS_H
#define TBEM6781829929292_INTERSECT_BALLS_H
#include <cstdlib>
#include <cassert>
#include <vector>
#include "vec.h"

namespace tbem {

template <size_t dim> struct Ball;

template <size_t dim>
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
    const std::vector<Ball<dim>>& ptsA,
    const std::vector<Ball<dim>>& ptsB);

template <size_t dim> struct Octree;
template <size_t dim>
std::vector<size_t>
intersect_balls(const Ball<dim>& ptA, const std::vector<Ball<dim>>& ptsB,
    const Octree<dim>& octB);

} // end namespace tbem
#endif
