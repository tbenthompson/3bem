#ifndef TBEMQQQQQQQ2222222_GEOMETRY_H
#define TBEMQQQQQQQ2222222_GEOMETRY_H
#include "vec_ops.h"
#include <cassert>

namespace tbem {

template <size_t dim> struct Ball;

/* A box class defined by its center and half_width. These are mainly 
 * used as bounding boxes for the nodes in the octree hierarchy.
 */
template <size_t dim>
struct Box {
    const Vec<double,dim> center;
    const Vec<double,dim> half_width;

    Box<dim> expand_by_max_axis_multiple(double factor) const;
    Box<dim> get_subcell(const Vec<size_t,dim>& idx) const;
    bool in_box(const Ball<dim>& b) const;
    bool in_box(const Vec<double,dim>& pt, const Vec<bool,dim>& inclusive) const;
    bool in_box_inclusive(const Vec<double,dim>& pt) const;
    size_t find_containing_subcell(const Vec<double,dim>& pt) const;

    static Box<dim> bounding_box(const std::vector<Ball<dim>>& x);
};

template <typename T>
Vec<T,3> outer_product(const Vec<double,3>& a, const T& b) 
{
    return {b * a[0], b * a[1], b * a[2]};
}

template <typename T>
Vec<T,2> outer_product(const Vec<double,2>& a, const T& b) 
{
    return {b * a[0], b * a[1]};
}

template <typename T>
Vec3<T> cross(const Vec3<T>& x, const Vec3<T>& y) 
{
    return {
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0]
    };
}

template <typename T, size_t dim>
T dot_product(const Vec<T,dim>& x, const Vec<T,dim>& y) 
{
    return sum(x * y);
}

template <typename T, size_t dim>
T hypot2(const Vec<T,dim>& v) 
{
    return dot_product(v, v);
}

template <typename T, size_t dim>
T hypot(const Vec<T,dim>& v) 
{
    return std::sqrt(hypot2(v));
}

template <typename T, size_t dim>
void normalize(Vec<T,dim>& v) 
{
    v /= hypot(v);
}

template <typename T, size_t dim>
Vec<T,dim> normalized(const Vec<T,dim>& v) 
{
    Vec<T,dim> res = v;
    normalize(res);
    return res;
}

template <typename T, size_t dim>
inline T dist2(const Vec<T,dim>& v0, const Vec<T,dim>& v1) 
{
    return hypot2(v1 - v0);
}

template <typename T, size_t dim>
inline T dist(const Vec<T,dim>& v0, const Vec<T,dim>& v1) 
{
    return hypot(v1 - v0);
}

template <typename T, size_t dim>
inline Vec<T,dim> unit(const int k);

template <>
inline Vec3<double> unit<double,3>(const int k) 
{
    Vec3<double> e_k = {0.0, 0.0, 0.0};
    e_k[k] = 1.0;
    return e_k;
}

template <>
inline Vec2<double> unit<double,2>(const int k) 
{
    Vec2<double> e_k = {0.0, 0.0};
    e_k[k] = 1.0;
    return e_k;
}

inline Vec3<double> 
unscaled_normal(const Vec<Vec3<double>,3>& corners) 
{
    return cross(corners[2] - corners[0], corners[2] - corners[1]);
}

inline Vec2<double> 
unscaled_normal(const Vec<Vec2<double>,2>& corners) 
{
    return {
        -(corners[1][1] - corners[0][1]),
        corners[1][0] - corners[0][0]
    };
}

template <size_t dim>
Vec<double,dim> facet_normal(const Vec<Vec<double,dim>,dim>& corners) 
{
    auto unscaled = unscaled_normal(corners);
    return normalized(unscaled);
}

inline double tri_area(const Vec3<double>& unscaled_normal) 
{
    return 0.5 * hypot(unscaled_normal);
}

inline double tri_area(const Vec<Vec3<double>,3>& corners) 
{
    return tri_area(unscaled_normal(corners));
}

enum Side {FRONT, INTERSECT, BEHIND};

/* Determine which side of the plane/line defined by triangle/segment 
 * the provided point is on.
 */
template <size_t dim>
Side which_side_point(const Vec<Vec<double,dim>,dim>& face,
                const Vec<double,dim>& pt) 
{
    auto normal = unscaled_normal(face);
    double dot_val = dot_product(pt - face[0], normal);
    if (dot_val > 0) { return FRONT; }
    else if (dot_val < 0) { return BEHIND; }
    else { return INTERSECT; }
}

/* Returns the side of a plane that a triangle/segment is on. */
template <size_t dim>
Side facet_side(const std::array<Side,dim>& s);

template <>
inline Side facet_side<2>(const std::array<Side,2>& s) 
{
    if (s[0] == s[1]) { return s[0]; } 
    else if(s[0] == INTERSECT) { return s[1]; }
    else if(s[1] == INTERSECT) { return s[0]; }
    else { return INTERSECT; }
}

template <>
inline Side facet_side<3>(const std::array<Side,3>& s) 
{
    auto edge0 = facet_side<2>({s[0], s[1]});
    auto edge1 = facet_side<2>({s[0], s[2]});
    auto edge2 = facet_side<2>({s[1], s[2]});
    if (edge0 == INTERSECT && edge1 == edge2) {
        return edge1;
    }
    if (edge1 == INTERSECT && edge2 == edge0) {
        return edge2;
    }
    if (edge2 == INTERSECT && edge0 == edge1) {
        return edge0;
    }
    return edge0;
}


/* Determine the side of the plane/line defined by triangle/segment
 * that the given triangle/segment is on
 */
template <size_t dim>
Side which_side_facet(const Vec<Vec<double,dim>,dim>& plane,
    const Vec<Vec<double,dim>,dim>& face) 
{
    std::array<Side,dim> sides;
    for (size_t d = 0; d < dim; d++) {
        sides[d] = which_side_point<dim>(plane, face[d]);
    }
    return facet_side<dim>(sides);
}

template <size_t dim>
Vec<double,dim> centroid(const Vec<Vec<double,dim>,dim>& f) 
{
    Vec<double,dim> centroid = zeros<Vec<double,dim>>::make(); 
    for (size_t d = 0; d < dim; d++) {
        centroid += f[d];
    }
    centroid /= static_cast<double>(dim);
    return centroid;
}

template <size_t dim>
struct Ball 
{
    Vec<double,dim> center;
    double radius;

    bool inside(const Vec<double,dim>& pt) const
    {
        return dist2(pt, center) < std::pow(radius, 2);
    }

    friend std::ostream& operator<<(std::ostream& os, const Ball<dim>& b) {
        os << "(" << b.center << ", " << b.radius << ")";
        return os;
    }
};

template <size_t dim>
bool balls_intersect(const Ball<dim>& a, const Ball<dim>& b) 
{
    return all(fabs(a.center - b.center) <= a.radius + b.radius);
}

template <size_t dim>
std::vector<Ball<dim>> balls_from_centers_radii(
    const std::vector<Vec<double,dim>>& centers, const std::vector<double> radii)
{
    assert(centers.size() == radii.size());
    std::vector<Ball<dim>> out(centers.size());
    for (size_t i = 0; i < centers.size(); i++) {
        out[i] = {centers[i], radii[i]};
    }
    return out;
}

//TODO: In the 3D case, this is not the minimum bounding ball. 
//see here: http://realtimecollisiondetection.net/blog/?p=20
//To test for minimum bounding ball, I could check that either
//2 or 3 points touch the bounding ball in 3D. If only one 
//touches, then the ball can be improved.
template <size_t dim>
Ball<dim> facet_ball(const Vec<Vec<double,dim>,dim>& f)
{
    auto c = centroid(f);
    double r2 = dist2(f[0], c);
    for (size_t d = 1; d < dim; d++) {
        r2 = std::max(r2, dist2(f[d], c));
    }
    return {c, std::sqrt(r2)};
}


} // end namespace tbem

#endif
