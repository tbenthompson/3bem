#include "closest_pt.h"
#include "vec_ops.h"
#include "numerics.h"

namespace tbem {

Vec<double,1> closest_pt_seg(const Vec<double,2>& pt, const Vec<Vec<double,2>,2> seg)
{
     auto v = seg[1] - seg[0];
     auto w = pt - seg[0];

     double c1 = dot_product(w, v);
     if (c1 <= 0) {
          return {-1.0};
     }

     double c2 = dot_product(v, v);
     if (c2 <= c1) {
          return {1.0};
     }

     double b = c1 / c2;
     return {b * 2 - 1};
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
    auto diff = pt - tri[0];
    auto edge0 = tri[1] - tri[0];
    auto edge1 = tri[2] - tri[0];
    auto a00 = dot_product(edge0, edge0);
    auto a01 = dot_product(edge0, edge1);
    auto a11 = dot_product(edge1, edge1);
    auto b0 = -dot_product(diff, edge0);
    auto b1 = -dot_product(diff, edge1);
    auto det = a00 * a11 - a01 * a01;
    auto t0 = a01 * b1 - a11 * b0;
    auto t1 = a01 * b0 - a00 * b1;

    if (t0 + t1 <= det)
    {
        if (t0 < 0)
        {
            if (t1 < 0)  // region 4
            {
                if (b0 < 0)
                {
                    t1 = 0;
                    if (-b0 >= a00)  // V0
                    {
                        t0 = 1;
                    }
                    else  // E01
                    {
                        t0 = -b0 / a00;
                    }
                }
                else
                {
                    t0 = 0;
                    if (b1 >= 0)  // V0
                    {
                        t1 = 0;
                    }
                    else if (-b1 >= a11)  // V2
                    {
                        t1 = 1;
                    }
                    else  // E20
                    {
                        t1 = -b1 / a11;
                    }
                }
            }
            else  // region 3
            {
                t0 = 0;
                if (b1 >= 0)  // V0
                {
                    t1 = 0;
                }
                else if (-b1 >= a11)  // V2
                {
                    t1 = 1;
                }
                else  // E20
                {
                    t1 = -b1 / a11;
                }
            }
        }
        else if (t1 < 0)  // region 5
        {
            t1 = 0;
            if (b0 >= 0)  // V0
            {
                t0 = 0;
            }
            else if (-b0 >= a00)  // V1
            {
                t0 = 1;
            }
            else  // E01
            {
                t0 = -b0 / a00;
            }
        }
        else  // region 0, interior
        {
            auto invDet = 1 / det;
            t0 *= invDet;
            t1 *= invDet;
        }
    }
    else
    {
        double tmp0, tmp1, numer, denom;

        if (t0 < 0)  // region 2
        {
            tmp0 = a01 + b0;
            tmp1 = a11 + b1;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - a01 * 2.0 + a11;
                if (numer >= denom)  // V1
                {
                    t0 = 1;
                    t1 = 0;
                }
                else  // E12
                {
                    t0 = numer / denom;
                    t1 = 1 - t0;
                }
            }
            else
            {
                t0 = 0;
                if (tmp1 <= 0)  // V2
                {
                    t1 = 1;
                }
                else if (b1 >= 0)  // V0
                {
                    t1 = 0;
                }
                else  // E20
                {
                    t1 = -b1 / a11;
                }
            }
        }
        else if (t1 < 0)  // region 6
        {
            tmp0 = a01 + b1;
            tmp1 = a00 + b0;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - a01 * 2.0 + a11;
                if (numer >= denom)  // V2
                {
                    t1 = 1;
                    t0 = 0;
                }
                else  // E12
                {
                    t1 = numer / denom;
                    t0 = 1 - t1;
                }
            }
            else
            {
                t1 = 0;
                if (tmp1 <= 0)  // V1
                {
                    t0 = 1;
                }
                else if (b0 >= 0)  // V0
                {
                    t0 = 0;
                }
                else  // E01
                {
                    t0 = -b0 / a00;
                }
            }
        }
        else  // region 1
        {
            numer = a11 + b1 - a01 - b0;
            if (numer <= 0)  // V2
            {
                t0 = 0;
                t1 = 1;
            }
            else
            {
                denom = a00 - a01 * 2.0 + a11;
                if (numer >= denom)  // V1
                {
                    t0 = 1;
                    t1 = 0;
                }
                else  // 12
                {
                    t0 = numer / denom;
                    t1 = 1 - t0;
                }
            }
        }
    }

    return {t0, t1};
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

} // end namespace tbem
