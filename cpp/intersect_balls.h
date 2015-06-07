#ifndef TBEM6781829929292_INTERSECT_BALLS_H
#define TBEM6781829929292_INTERSECT_BALLS_H
#include <cstdlib>
#include "octree.h"
#include "vec.h"

namespace tbem {

template <size_t dim>
std::vector<std::pair<size_t,size_t>> intersect_balls_all_pairs(
    const std::vector<Vec<double,dim>>& ptsA,
    const std::vector<Vec<double,dim>>& ptsB,
    const Octree<dim>& cellA,
    const Octree<dim>& cellB, 
    double r);

template <size_t dim>
std::vector<size_t>
intersect_balls(const Vec<double,dim>& ptA, 
    const std::vector<Vec<double,dim>>& ptsB,
    const Octree<dim>& octB, double r);

} // end namespace tbem
#endif
