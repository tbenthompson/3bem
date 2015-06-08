#ifndef TBEM_12312312_OCTREE_H
#define TBEM_12312312_OCTREE_H

#include <array>
#include <vector>
#include <memory>
#include "vec.h"

namespace tbem { 

template <size_t dim> struct Ball;
template <size_t dim> struct Box;

template <size_t dim>
struct Octree {
    static const size_t split = 2<<(dim-1);
    typedef std::array<std::unique_ptr<Octree<dim>>,split> ChildrenType;
    const Box<dim> bounds;
    const Box<dim> true_bounds;
    const std::vector<size_t> indices;
    const size_t level;
    const size_t index;
    ChildrenType children;

    size_t n_immediate_children() const; 
    size_t n_children() const; 
    bool is_leaf() const; 
    size_t find_containing_child(const Vec<double,dim>& pt) const;
};

template <size_t dim>
Vec<size_t,dim> make_child_idx(size_t i);

template <size_t dim>
Octree<dim>
make_octree(const std::vector<Ball<dim>>& pts, size_t min_pts_per_cell);

template <size_t dim>
Octree<dim> 
make_octree(const std::vector<Vec<double,dim>>& pts, size_t min_pts_per_cell);

} // END namespace tbem
#endif
