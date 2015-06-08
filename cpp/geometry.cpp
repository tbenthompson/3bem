#include "geometry.h"

namespace tbem {

template <size_t dim>
std::ostream& operator<<(std::ostream& os, const Box<dim>& obj)
{
    os << "{Box center=" << obj.center << ", half_width=" << obj.half_width << "}";
    return os;
}

template <size_t dim>
Box<dim> Box<dim>::expand_by_max_axis_multiple(double factor) const
{
    auto max_axis = max(half_width);
    return {center, half_width + factor * max_axis};
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
bool Box<dim>::in_box(const Ball<dim>& b) const
{
    return all(center - half_width <= b.center - b.radius) &&
        all(center + half_width >= b.center + b.radius);
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
Box<dim> Box<dim>::bounding_box(const std::vector<Ball<dim>>& pts)
{
    if (pts.size() == 0) {
        auto z = zeros<Vec<double,dim>>::make();
        return {z, z};
    }
    auto min_corner = pts[0].center - pts[0].radius;
    auto max_corner = pts[0].center + pts[0].radius;
    for (int i = 1; i < pts.size(); ++i) {
        auto c = pts[i].center;
        auto r = pts[i].radius;
        for (unsigned int d = 0; d < dim; ++d) {
            min_corner[d] = std::min(min_corner[d], c[d] - r);
            max_corner[d] = std::max(max_corner[d], c[d] + r);
        }
    }
    auto center = (max_corner + min_corner) / 2.0;
    auto half_width = (max_corner - min_corner) / 2.0;
    return {center, half_width};
}

template <size_t dim>
size_t Box<dim>::find_containing_subcell(const Vec<double,dim>& pt) const
{
    size_t child_idx = 0;
    for (size_t d = 0; d < dim; d++) {
        if (pt[d] > center[d]) {
            child_idx++; 
        }
        if (d < dim - 1) {
            child_idx = child_idx << 1;
        }
    }
    return child_idx;
}

template class Box<2>;
template class Box<3>;

} // end namespace tbem
