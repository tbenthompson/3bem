#include "octree.h"
#include "numerics.h"
#include "vec_ops.h"
#include "util.h"
#include <assert.h>

namespace tbem {

template <size_t dim>
std::ostream& operator<<(std::ostream& os, const Box<dim>& obj)
{
    os << "{Box center=" << obj.center << ", half_width=" << obj.half_width << "}";
    return os;
}
    
template <size_t dim>
Box<dim> Box<dim>::get_subcell(const Vec<size_t,dim>& idx) const
{
    auto new_halfwidth = half_width / 2.0;
    auto new_center = center;
    for (size_t d = 0; d < dim; d++) {
        new_center[d] += ((static_cast<double>(idx[d]) * 2) - 1) * new_halfwidth[d];
    }
    return {new_center, new_halfwidth};
}

template <>
bool Box<2>::in_box(const Vec<double,2>& pt) const
{
    return fabs(pt[0] - center[0]) < half_width[0] && 
           fabs(pt[1] - center[1]) < half_width[1];
}
template <>
bool Box<3>::in_box(const Vec<double,3>& pt) const
{
    return fabs(pt[0] - center[0]) < half_width[0] && 
           fabs(pt[1] - center[1]) < half_width[1] &&
           fabs(pt[2] - center[2]) < half_width[2];
}

template <size_t dim>
Box<dim> Box<dim>::bounding_box(const std::vector<Vec<double,dim>>& x)
{
    auto min_corner = x[0];
    auto max_corner = x[0];
    for (int i = 0; i < x.size(); ++i) {
        for (unsigned int d = 0; d < dim; ++d) {
            min_corner[d] = std::min(min_corner[d], x[i][d]);
            max_corner[d] = std::max(max_corner[d], x[i][d]);
        }
    }
    auto center = (max_corner + min_corner) / 2.0;
    auto half_width = (max_corner - min_corner) / 2.0;
    return {center, half_width};
}
template class Box<2>;
template class Box<3>;

template <size_t dim>
Vec<size_t,dim> make_child_idx(size_t i) 
{
    Vec<size_t,dim> child_idx;
    for (int d = dim - 1; d >= 0; d--) {
        auto idx = i % 2;
        i = (i - idx) / 2;
        child_idx[d] = idx;
    }
    return child_idx;
}
template 
Vec<size_t,2> make_child_idx<2>(size_t i);
template 
Vec<size_t,3> make_child_idx<3>(size_t i);

template <size_t dim>
std::unique_ptr<Octree<dim>> 
make_child(const Vec<size_t,dim>& idx, const Box<dim>& bounds,
    const std::vector<Vec<double,dim>>& pts, const std::vector<int> indices,
    size_t level, size_t min_pts_per_cell) 
{
    auto child_bounds = bounds.get_subcell(idx);
    std::vector<int> child_indices;
    for (size_t i = 0; i < indices.size(); i++) {
        if (child_bounds.in_box(pts[indices[i]])) {
            child_indices.push_back(indices[i]);
        }
    }
    if (child_indices.size() == 0) {
        return nullptr;
    }
    auto child_level = level + 1;

    auto children = build_children(
        child_bounds, pts, child_indices, child_level, min_pts_per_cell
    );
    return std::unique_ptr<Octree<dim>>(new Octree<dim>{
        {child_bounds, child_indices, child_level},
        std::move(children)
    });
}

template <size_t dim>
typename Octree<dim>::ChildrenType
build_children(const Box<dim>& bounds,
    const std::vector<Vec<double,dim>>& pts, const std::vector<int> indices,
    size_t level, size_t min_pts_per_cell) 
{
    typename Octree<dim>::ChildrenType children;
    if (indices.size() <= min_pts_per_cell) {
        return std::move(children);
    }
#pragma omp parallel for if(level == 0)
    for (size_t i = 0; i < Octree<dim>::split; i++) {
        auto idx = make_child_idx<dim>(i);
        children[i] = make_child(idx, bounds, pts, indices, level, min_pts_per_cell);
    }
    return std::move(children);
}

template <size_t dim>
Octree<dim> build_octree(const std::vector<Vec<double,dim>>& pts,
    size_t min_pts_per_cell)
{
    auto true_box = Box<dim>::bounding_box(pts);
    const double expand_factor = 1.0 + 1e-5;
    Box<dim> expanded_box = {true_box.center, true_box.half_width * expand_factor};
    auto all_indices = range(pts.size());
    auto children = build_children(
        expanded_box, pts, all_indices, 0, min_pts_per_cell
    );
    return Octree<dim>{{expanded_box, all_indices, 0}, std::move(children)};
}

template 
Octree<2> build_octree(const std::vector<Vec<double,2>>& pts,
    size_t min_pts_per_cell);
template 
Octree<3> build_octree(const std::vector<Vec<double,3>>& pts,
    size_t min_pts_per_cell);

} // END namespace tbem
