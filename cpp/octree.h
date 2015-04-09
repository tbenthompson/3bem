#ifndef _OCTREE2_H
#define _OCTREE2_H

#include <array>
#include <vector>
#include <memory>
#include "vec.h"
#include "numbers.h"

namespace tbem { 

/* A box class defined by its center and half_width. These are used as
 * bounding boxes for the nodes in the octree hierarchy.
 */
template <size_t dim>
class Box {
public:
    const Vec<double,dim> center;
    const Vec<double,dim> half_width;
    Box<dim> get_subcell(const Vec<size_t,dim>& idx) const;
    bool in_box(const Vec<double,dim>& pt) const;

    static Box<dim> bounding_box(const std::vector<Vec<double,dim>>& x);
};

template <size_t dim>
struct Octree {
    static const size_t split = 2<<(dim-1);
    typedef std::array<std::unique_ptr<Octree<dim>>,split> ChildrenType;
    const Box<dim> bounds;
    const std::vector<int> indices;
    const ChildrenType children;
    const size_t level;
};

template <size_t dim>
Vec<size_t,dim> make_child_idx(size_t i);

template <size_t dim>
std::unique_ptr<Octree<dim>> build_octree(const std::vector<Vec<double,dim>>& pts,
    size_t min_pts_per_cell);

} // END namespace tbem
#endif
