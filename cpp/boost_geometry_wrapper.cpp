#include "boost_geometry_wrapper.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/segment.hpp> 
#include <boost/geometry/algorithms/intersection.hpp>

namespace tbem {

std::vector<Vec<double,2>> seg_seg_intersection(const Vec<Vec<double,2>,2>& A,
    const Vec<Vec<double,2>,2>& B)
{
    typedef boost::geometry::model::d2::point_xy<double> Point;
    typedef boost::geometry::model::segment<Point> Segment;
    Segment seg0(Point(A[0][0],A[0][1]), Point(A[1][0], A[1][1]));
    Segment seg1(Point(B[0][0],B[0][1]), Point(B[1][0], B[1][1]));
    
    std::vector<Point> intersect_pts; 
    boost::geometry::intersection(seg0, seg1, intersect_pts);

    std::vector<Vec<double,2>> out;
    for (auto p: intersect_pts) {
        out.push_back({p.x(), p.y()});
    }
    return out;
}

std::vector<Vec<double,3>> seg_tri_intersection(const Vec<Vec<double,3>,3>& A,
    const Vec<Vec<double,3>,2>& B)
{
    assert(false);//NOT IMPLEMENTED
    typedef boost::geometry::model::point<
        double,3,boost::geometry::cs::cartesian> Point;
    typedef boost::geometry::model::polygon<Point> Polygon;
    typedef boost::geometry::model::segment<Point> Segment;
    Polygon tri;
    tri.outer().push_back(Point(A[0][0], A[0][1], A[0][2]));
    tri.outer().push_back(Point(A[1][0], A[1][1], A[1][2]));
    tri.outer().push_back(Point(A[2][0], A[2][1], A[2][2]));
    Segment seg(Point(B[0][0], B[0][1], B[0][2]), Point(B[1][0], B[1][1], B[1][2]));
    
    std::vector<Point> intersect_pts; 
    // boost::geometry::intersection(tri, seg, intersect_pts);

    std::vector<Vec<double,3>> out;
    // for (auto p: intersect_pts) {
    //     out.push_back({p.get<0>(), p.get<1>(), p.get<2>()});
    // }
    return out;
}

template <> 
std::vector<Vec<double,2>> seg_facet_intersection(const Vec<Vec<double,2>,2>& f,
    const Vec<Vec<double,2>,2>& seg)
{
    return seg_seg_intersection(f, seg);
}

template <> 
std::vector<Vec<double,3>> seg_facet_intersection(const Vec<Vec<double,3>,3>& f,
    const Vec<Vec<double,3>,2>& seg) 
{
    return seg_tri_intersection(f, seg);
}

bool in_polygon(const std::vector<Vec<double,2>>& poly,
    const Vec<double,2>& pt) 
{
    typedef boost::geometry::model::d2::point_xy<double> point_type;
    typedef boost::geometry::model::polygon<point_type> polygon_type;
    polygon_type boost_geom_poly;
    for (size_t i = 0; i < poly.size(); i++) {
        boost_geom_poly.outer().push_back({poly[i][0], poly[i][1]});
    }
    point_type boost_geom_pt(pt[0], pt[1]);
    return boost::geometry::within(boost_geom_pt, boost_geom_poly);
}

} //end namespace tbem
