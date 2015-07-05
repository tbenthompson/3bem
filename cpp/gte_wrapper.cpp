#include "gte_wrapper.h"
#include <GTEngine.h>
#include "geometry.h"

namespace tbem {

template <size_t dim>
gte::Vector<dim,double> gte_vec(const Vec<double,dim>& in)
{
    gte::Vector<dim,double> out;
    for (size_t d = 0; d < dim; d++) {
        out[d] = in[d];
    }
    return out;
}

template <int dim>
Vec<double,dim> tbem_vec(const gte::Vector<dim,double>& in)
{
    Vec<double,dim> out;
    for (size_t d = 0; d < dim; d++) {
        out[d] = in[d];
    }
    return out;
}

NearestPoint<2> closest_pt_seg(const Vec<double,2>& pt, const Vec<Vec<double,2>,2> seg)
{
    gte::Segment<2,double> gte_seg(gte_vec(seg[0]), gte_vec(seg[1]));
    auto gte_pt = gte_vec(pt);
    gte::DCPQuery<double,gte::Vector<2,double>,gte::Segment<2,double>> Q;
    auto result = Q(gte_pt, gte_seg);
    return {
        {2 * result.segmentParameter - 1},
        tbem_vec(result.segmentClosest),
        result.distance
    };
}

NearestPoint<3> closest_pt_tri(const Vec<double,3>& pt, const Vec<Vec<double,3>,3> tri) 
{
    gte::Triangle3<double> gte_tri(gte_vec(tri[0]), gte_vec(tri[1]), gte_vec(tri[2]));
    auto gte_pt = gte_vec(pt);
    gte::DCPQuery<double,gte::Vector<3,double>,gte::Triangle<3,double>> Q;
    auto result = Q(gte_pt, gte_tri);
    return {
        {result.parameter[1], result.parameter[2]},
        tbem_vec(result.closest),
        result.distance
    };
}


template <>
NearestPoint<2> closest_pt_facet(const Vec<double,2>& pt,
    const Vec<Vec<double,2>,2> tri) 
{
    return closest_pt_seg(pt, tri);
}

template <>
NearestPoint<3> closest_pt_facet(const Vec<double,3>& pt,
    const Vec<Vec<double,3>,3> tri) 
{
    return closest_pt_tri(pt, tri);
}

std::vector<Vec<double,2>> seg_seg_intersection(const Vec<Vec<double,2>,2>& A,
    const Vec<Vec<double,2>,2>& B)
{
    gte::Segment<2,double> seg0(gte_vec(A[0]), gte_vec(A[1]));
    gte::Segment<2,double> seg1(gte_vec(B[0]), gte_vec(B[1]));
    gte::FIQuery<double,gte::Segment2<double>,gte::Segment2<double>> Q;
    auto result = Q(seg0, seg1);

    std::vector<Vec<double,2>> out(result.numIntersections);
    for (size_t i = 0; i < static_cast<size_t>(result.numIntersections); i++) {
        out[i] = tbem_vec(result.point[i]);
    }
    return out;
}

Vec<Vec<double,3>,2> choose_planar_axes(const Vec<double,3>& plane_normal)
{
    Vec<double,3> first_axis{1, 0, 0};
    auto first_axis_dot = dot_product(plane_normal, first_axis);
    if (first_axis_dot > 0.95) {
        first_axis = {0, 1, 0}; 
        first_axis_dot = dot_product(plane_normal, first_axis);
    }
    first_axis = normalized(first_axis - first_axis_dot * plane_normal);
    auto second_axis = cross(plane_normal, first_axis);
    return {first_axis, second_axis};
}


Vec<double,2> to_planar_coords(const Vec<double,3>& origin,
    const Vec<Vec<double,3>,2>& axes, const Vec<double,3>& pt) 
{
    return {
        dot_product(pt - origin, axes[0]),
        dot_product(pt - origin, axes[1])
    };
}

Vec<double,3> from_planar_coords(const Vec<double,3>& origin,
    const Vec<Vec<double,3>,2>& axes, const Vec<double,2>& pt) 
{
    return pt[0] * axes[0] + pt[1] * axes[1] + origin;
}

bool point_already_in_set(const std::vector<Vec<double,3>>& set,
    const Vec<double,3>& pt, double threshold) 
{
    for (size_t j = 0; j < set.size(); j++) {
        if (dist(set[j], pt) < threshold) {
            return true; 
        }
    }
    return false;
}

void add_seg_seg_intersections(const Vec<double,3>& origin,
    const Vec<Vec<double,3>,2>& axes, const Vec<Vec<double,2>,2>& seg0,
    const Vec<Vec<double,2>,2>& seg1, std::vector<Vec<double,3>>& intersections)
{
    auto new_pts = seg_seg_intersection(seg0, seg1);
    double threshold = 1e-14 * facet_ball(seg1).radius;
    for (size_t i = 0; i < new_pts.size(); i++) {
        auto new_pt = from_planar_coords(origin, axes, new_pts[i]);
        if (point_already_in_set(intersections, new_pt, threshold)) {
            continue;
        }
        intersections.push_back(new_pt); 
    }
}

std::vector<Vec<double,3>> seg_tri_intersection(const Vec<Vec<double,3>,3>& A,
    const Vec<Vec<double,3>,2>& B)
{
    gte::Triangle3<double> tri(gte_vec(A[0]), gte_vec(A[1]), gte_vec(A[2]));
    gte::Segment<3,double> seg(gte_vec(B[0]), gte_vec(B[1]));
    gte::FIQuery<double,gte::Segment3<double>,gte::Triangle3<double>> Q;
    auto result = Q(seg, tri);
    if (result.intersect) {
        return {tbem_vec(result.point)};
    } 
    
    // GTE unfortunately does not handle the case in which the segment is 
    // coplanar with the triangle
    // To check manually, project everything into the plane of the triangle,
    // then perform a 2D triangle, line segment intersection test.
    auto plane_normal = normalized(cross(A[2] - A[0], A[1] - A[0]));

    auto dist_to_plane0 = dot_product(plane_normal, B[0] - A[0]);
    auto dist_to_plane1 = dot_product(plane_normal, B[1] - A[0]);

    // If not coplanar, just return now.
    double threshold = 1e-14 * facet_ball(A).radius;
    if (std::fabs(dist_to_plane0) > threshold ||
        std::fabs(dist_to_plane1) > threshold) 
    {
        return {}; 
    }

    auto axes = choose_planar_axes(plane_normal);
    auto origin = A[0]; 

    auto seg0 = to_planar_coords(origin, axes, B[0]);
    auto seg1 = to_planar_coords(origin, axes, B[1]);
    auto tri0 = to_planar_coords(origin, axes, A[0]);
    auto tri1 = to_planar_coords(origin, axes, A[1]);
    auto tri2 = to_planar_coords(origin, axes, A[2]);

    std::vector<Vec<double,3>> intersections;
    add_seg_seg_intersections(
        origin, axes, {seg0, seg1}, {tri0, tri1}, intersections
    );
    add_seg_seg_intersections(
        origin, axes, {seg0, seg1}, {tri1, tri2}, intersections
    );
    add_seg_seg_intersections(
        origin, axes, {seg0, seg1}, {tri2, tri0}, intersections
    );

    if (closest_pt_tri(B[0], A).distance < threshold &&
        !point_already_in_set(intersections, B[0], threshold)) 
    {
        intersections.push_back(B[0]);
    }

    if (closest_pt_tri(B[1], A).distance < threshold &&
        !point_already_in_set(intersections, B[1], threshold)) 
    {
        intersections.push_back(B[1]);
    }

    return intersections;
}

void add_intersections(const Vec<Vec<double,3>,3>& tri,
    const Vec<Vec<double,3>,2>& seg, std::vector<Vec<double,3>>& intersections)
{
    auto new_pts = seg_tri_intersection(tri, seg);
    double threshold = 1e-14 * facet_ball(tri).radius;
    for (size_t i = 0; i < new_pts.size(); i++) {
        if (point_already_in_set(intersections, new_pts[i], threshold)) {
            continue;
        }
        intersections.push_back(new_pts[i]); 
    }
}

std::vector<Vec<double,3>> tri_tri_intersection(const Vec<Vec<double,3>,3>& A,
    const Vec<Vec<double,3>,3>& B)
{
    std::vector<Vec<double,3>> intersections;
    add_intersections(A, {B[0], B[1]}, intersections);
    add_intersections(A, {B[1], B[2]}, intersections);
    add_intersections(A, {B[2], B[0]}, intersections);
    add_intersections(B, {A[0], A[1]}, intersections);
    add_intersections(B, {A[1], A[2]}, intersections);
    add_intersections(B, {A[2], A[0]}, intersections);
    assert(intersections.size() <= 3);
    return intersections;
}

template <> 
std::vector<Vec<double,2>> seg_facet_intersection(
    const Vec<Vec<double,2>,2>& f, const Vec<Vec<double,2>,2>& seg)
{
    return seg_seg_intersection(f, seg);
}

template <> 
std::vector<Vec<double,3>> seg_facet_intersection(
    const Vec<Vec<double,3>,3>& f, const Vec<Vec<double,3>,2>& seg) 
{
    return seg_tri_intersection(f, seg);
}

template <> 
std::vector<Vec<double,2>> facet_facet_intersection(
    const Vec<Vec<double,2>,2>& fA, const Vec<Vec<double,2>,2>& fB)
{
    return seg_seg_intersection(fA, fB);
}

template <> 
std::vector<Vec<double,3>> facet_facet_intersection(
    const Vec<Vec<double,3>,3>& fA, const Vec<Vec<double,3>,3>& fB)
{
    return tri_tri_intersection(fA, fB);
}

bool in_polygon(const std::vector<Vec<double,2>>& poly,
    const Vec<double,2>& pt) 
{
    std::vector<gte::Vector<2,double>> gte_pts;
    for (size_t i = 0; i < poly.size(); i++) {
        gte_pts.push_back({poly[i][0], poly[i][1]});
    }
    gte::PointInPolygon2<double> gte_poly(gte_pts.size(), gte_pts.data());
    return gte_poly.Contains({pt[0], pt[1]});
}

template <size_t dim>
bool is_intersection_box_ball(const Box<dim>& box, const Ball<dim>& ball)
{
    gte::TIQuery<double,gte::OrientedBox<dim,double>,gte::Hypersphere<dim,double>> Q;
    gte::Hypersphere<dim,double> S(gte_vec(ball.center), ball.radius);
    std::array<gte::Vector<dim,double>,dim> axes;
    for (size_t d = 0; d < dim; d++) {
        axes[d] = gte_vec(unit<double,dim>(d));
    }
    gte::OrientedBox<dim,double> B(
        gte_vec(box.center), axes, gte_vec(box.half_width)
    );
    return Q(B, S).intersect;
}

template 
bool is_intersection_box_ball<2>(const Box<2>& box, const Ball<2>& ball);
template 
bool is_intersection_box_ball<3>(const Box<3>& box, const Ball<3>& ball);

} //end namespace tbem
