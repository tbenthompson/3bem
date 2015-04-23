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

template <size_t dim>
bool Box<dim>::in_box(const Vec<double,dim>& pt, const Vec<bool,dim>& inclusive) const
{
    bool in = true;
    for (size_t d = 0; d < dim; d++) {
        auto sep = fabs(pt[d] - center[d]);
        in = in && (
            (sep < half_width[d]) ||
            (inclusive[d] && (sep == half_width[d])));
    }
    return in;
}

template <size_t dim>
bool Box<dim>::in_box_inclusive(const Vec<double,dim>& pt) const
{
    Vec<bool,dim> inclusive;
    for (size_t d = 0; d < dim; d++) {
        inclusive[d] = true;
    }
    return in_box(pt, inclusive);
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
        i = i >> 1;
        child_idx[d] = idx;
    }
    return child_idx;
}
template 
Vec<size_t,2> make_child_idx<2>(size_t i);
template 
Vec<size_t,3> make_child_idx<3>(size_t i);

template <size_t dim>
typename Octree<dim>::ChildrenType
build_children(const Box<dim>& bounds, 
    const std::vector<Vec<double,dim>>& pts, const std::vector<int>& indices,
    size_t level, size_t min_pts_per_cell) 
{
    typename Octree<dim>::ChildrenType children;
    if (indices.size() <= min_pts_per_cell) {
        return std::move(children);
    }


    std::array<std::vector<int>,Octree<dim>::split> child_indices;

    for (size_t i = 0; i < indices.size(); i++) {
        auto pt = pts[indices[i]];
        size_t child_idx = 0;
        for (size_t d = 0; d < dim; d++) {
            if (pt[d] > bounds.center[d]) {
                child_idx++; 
            }
            if (d < dim - 1) {
                child_idx = child_idx << 1;
            }
        }

        child_indices[child_idx].push_back(indices[i]);
    }

    auto child_level = level + 1;
#pragma omp parallel for if(level == 0)
    for (size_t i = 0; i < Octree<dim>::split; i++) {
        if (child_indices[i].size() == 0) {
            continue;
        }
        auto idx = make_child_idx<dim>(i);
        auto child_bounds = bounds.get_subcell(idx);

        std::vector<Vec<double,dim>> tight_pts(child_indices[i].size());
        for (size_t pt_idx = 0; pt_idx < child_indices[i].size(); pt_idx++) {
            tight_pts[pt_idx] = pts[child_indices[i][pt_idx]];
        }
        auto tight_bounds = Box<dim>::bounding_box(tight_pts);

        typename Octree<dim>::ChildrenType sub_children;
        if (any(tight_bounds.half_width != 0.0)) {
            sub_children = build_children(
                child_bounds, pts, child_indices[i], child_level, min_pts_per_cell
            );
        } 

        children[i] = std::unique_ptr<Octree<dim>>(new Octree<dim>{
            {child_bounds, tight_bounds, child_indices[i], child_level},
            std::move(sub_children)
        });
    }
    return std::move(children);
}

template <size_t dim>
Octree<dim> build_octree(const std::vector<Vec<double,dim>>& pts,
    size_t min_pts_per_cell)
{
    auto box = Box<dim>::bounding_box(pts);
    auto half_width = box.half_width;
    auto max_axis = max(half_width);
    half_width += 1e-5 * max_axis;
    Box<dim> expanded_box{box.center, half_width};
    auto all_indices = range(pts.size());
    auto children = build_children(
        box, pts, all_indices, 0, min_pts_per_cell
    );
    return Octree<dim>{
        {expanded_box, expanded_box, all_indices, 0},
        std::move(children)
    };
}

template 
Octree<2> build_octree(const std::vector<Vec<double,2>>& pts,
    size_t min_pts_per_cell);
template 
Octree<3> build_octree(const std::vector<Vec<double,3>>& pts,
    size_t min_pts_per_cell);

} // END namespace tbem
