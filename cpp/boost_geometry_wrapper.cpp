#include "boost_geometry_wrapper.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/segment.hpp> 
#include <boost/geometry/algorithms/intersection.hpp>
#include <gte/Include/GTEngine.h>

//TODO: Remove the remainder of the boost geometry junk and rename to gte_wrapper
namespace tbem {

std::vector<Vec<double,2>> seg_seg_intersection(const Vec<Vec<double,2>,2>& A,
    const Vec<Vec<double,2>,2>& B)
{
    gte::Segment<2,double> seg0(
        gte::Vector<2,double>({A[0][0], A[0][1]}),
        gte::Vector<2,double>({A[1][0], A[1][1]})
    );
    gte::Segment<2,double> seg1(
        gte::Vector<2,double>({B[0][0], B[0][1]}),
        gte::Vector<2,double>({B[1][0], B[1][1]})
    );
    gte::FIQuery<double,gte::Segment2<double>,gte::Segment2<double>> Q;
    auto result = Q(seg0, seg1);

    if (result.numIntersections == 0) {
        return {};
    } else if (result.numIntersections == 1) {
        return {Vec<double,2>{result.point[0][0], result.point[0][1]}};
    } else if (result.numIntersections == 2) {
        return {
            Vec<double,2>{result.point[0][0], result.point[0][1]},
            Vec<double,2>{result.point[1][0], result.point[1][1]}
        };
    }
}

std::vector<Vec<double,3>> seg_tri_intersection(const Vec<Vec<double,3>,3>& A,
    const Vec<Vec<double,3>,2>& B)
{
    gte::Triangle3<double> tri(
        gte::Vector<3,double>({A[0][0], A[0][1], A[0][2]}),
        gte::Vector<3,double>({A[1][0], A[1][1], A[1][2]}),
        gte::Vector<3,double>({A[2][0], A[2][1], A[2][2]})
    );
    gte::Segment<3,double> seg(
        gte::Vector<3,double>({B[0][0], B[0][1], B[0][2]}),
        gte::Vector<3,double>({B[1][0], B[1][1], B[1][2]})
    );
    gte::FIQuery<double,gte::Segment3<double>,gte::Triangle3<double>> Q;
    auto result = Q(seg, tri);
    if (result.intersect) {
        return {Vec<double,3>{result.point[0], result.point[1], result.point[2]}};
    } else {
        return {};
    }
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
