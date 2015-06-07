#include "gte_wrapper.h"
#include <gte/Include/GTEngine.h>

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
    } else {
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
    std::vector<gte::Vector<2,double>> gte_pts;
    for (size_t i = 0; i < poly.size(); i++) {
        gte_pts.push_back({poly[i][0], poly[i][1]});
    }
    gte::PointInPolygon2<double> gte_poly(gte_pts.size(), gte_pts.data());
    return gte_poly.Contains({pt[0], pt[1]});
}

Vec<double,1> closest_pt_seg(const Vec<double,2>& pt, const Vec<Vec<double,2>,2> seg)
{
    gte::Segment<2,double> gte_seg(
        gte::Vector<2,double>({seg[0][0], seg[0][1]}),
        gte::Vector<2,double>({seg[1][0], seg[1][1]})
    );
    gte::Vector<2,double> gte_pt({pt[0], pt[1]});
    gte::DCPQuery<double,gte::Vector<2,double>,gte::Segment<2,double>> Q;
    auto result = Q(gte_pt, gte_seg);
    return {2 * result.segmentParameter - 1};
}

/* Find the closest reference point on a triangle. I took this from somewhere 
 * online and can't remember where. 
 *
 * A faster approach might be to treat the problem as a quadratic minimization
 * for the triangle reference coordinates. (see Ericson's Real-Time Collision
 * Detection)
 */
Vec<double,2> closest_pt_tri(const Vec<double,3>& pt, const Vec<Vec<double,3>,3> tri) 
{
    gte::Triangle3<double> gte_tri(
        gte::Vector<3,double>({tri[0][0], tri[0][1], tri[0][2]}),
        gte::Vector<3,double>({tri[1][0], tri[1][1], tri[1][2]}),
        gte::Vector<3,double>({tri[2][0], tri[2][1], tri[2][2]})
    );
    gte::Vector<3,double> gte_pt({pt[0], pt[1], pt[2]});
    gte::DCPQuery<double,gte::Vector<3,double>,gte::Triangle<3,double>> Q;
    auto result = Q(gte_pt, gte_tri);
    return {result.parameter[1], result.parameter[2]};
}

template <>
Vec<double,1> closest_pt_facet(const Vec<double,2>& pt,
    const Vec<Vec<double,2>,2> tri) 
{
    return closest_pt_seg(pt, tri);
}

template <>
Vec<double,2> closest_pt_facet(const Vec<double,3>& pt,
    const Vec<Vec<double,3>,3> tri) 
{
    return closest_pt_tri(pt, tri);
}

} //end namespace tbem
