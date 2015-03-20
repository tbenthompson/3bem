/* Heavily modified from original version obtained from
 * https://github.com/nblinov/Adaptive-Integrator
 *
 * The original license:
 *
 * The MIT License (MIT)
 * 
 * Copyright (c) 2013 Nikita Blinov
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifndef __QWkJELHFUUASASDK_ADAPTIVE_QUAD_H
#define __QWkJELHFUUASASDK_ADAPTIVE_QUAD_H

#include <cmath>
#include <limits>
#include <iostream>
#include "vec_ops.h"
#include "numerics.h"

namespace tbem {
using std::fabs;

const double lobatto_alpha = std::sqrt(2./3.);
const double lobatto_beta = 1./std::sqrt(5.);
const double lobatto_x1 = .94288241569547971905635175843185720232;
const double lobatto_x2 = .64185334234578130578123554132903188354;
const double lobatto_x3 = .23638319966214988028222377349205292599;

static int adaptive_quad_n = 0;
template <typename T, typename F>
T adaptlobstp(const F& f, const double a, const double b, 
              const T& fa, const T& fb, const T& is )
{
    // std::cout << a << " " << b << std::endl;
    adaptive_quad_n++;
    double m = (a + b) / 2.; 
    double h = (b - a) / 2.;

    double ah = lobatto_alpha * h;
    double bh = lobatto_beta * h;
    double mll = m - ah;
    double ml = m - bh;
    double mr = m + bh;
    double mrr = m + ah;

    T fmll = f(mll);
    T fml = f(ml);
    T fm = f(m);
    T fmr = f(mr);
    T fmrr = f(mrr);

    T i2 = (h / 6.) * (fa + fb + 5.0 * (fml + fmr));    
    T i1 = (h / 1470.) * (
            77.0 * (fa + fb) + 
            432.0 * (fmll + fmrr) +
            625.0 * (fml + fmr) +
            672.0 * fm);

    if (all(is + (i1 - i2) == is)) {
        return i1;
    } else if (mll <= a or b <= mrr) {
        return i1;
    } else {
        return adaptlobstp(f, a, mll, fa, fmll, is)
             + adaptlobstp(f, mll, ml, fmll, fml, is)
             + adaptlobstp(f, ml, m, fml, fm, is)
             + adaptlobstp(f, m, mr, fm, fmr, is)
             + adaptlobstp(f, mr, mrr, fmr, fmrr, is)
             + adaptlobstp(f, mrr, b, fmrr, fb, is);
    }
}

inline double get_error_is(double p_tol, double erri1, double erri2, double is, 
                    double a, double b) {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const double R = ((erri2 == 0.0) ? 1.0 : (erri1 / erri2));
    const double tol = ((R > 0.0 and R < 1.0) ? (p_tol / R) : p_tol);

    double retval = fabs(is) * tol / eps;
    if (fabs(retval) < p_tol) {
        retval = b - a;
    }
    return retval;
}

template <typename T>
Vec1<T> get_error_is(double p_tol, Vec1<T> erri1, Vec1<T> erri2,
                          Vec1<T> is, double a, double b) {
    return {get_error_is(p_tol, erri1[0], erri2[0], is[0], a, b)};
}

template <typename T>
Vec2<T> get_error_is(double p_tol, Vec2<T> erri1, Vec2<T> erri2,
                          Vec2<T> is, double a, double b) {
    return {
        get_error_is(p_tol, erri1[0], erri2[0], is[0], a, b),
        get_error_is(p_tol, erri1[1], erri2[1], is[1], a, b)
    };
}

template <typename T>
Vec3<T> get_error_is(double p_tol, Vec3<T> erri1, Vec3<T> erri2,
                          Vec3<T> is, double a, double b) {
    return {
        get_error_is(p_tol, erri1[0], erri2[0], is[0], a, b),
        get_error_is(p_tol, erri1[1], erri2[1], is[1], a, b),
        get_error_is(p_tol, erri1[2], erri2[2], is[2], a, b)
    };
}

template <typename T, typename F>
T adaptive_integrate(const F& f, double a,
                     double b, double p_tol) {
    double m = (a + b) / 2.; 
    double h = (b - a) / 2.;

    const T y[13] = {
        f(a),
        f(m - lobatto_x1 * h),
        f(m - lobatto_alpha * h),
        f(m - lobatto_x2 * h),
        f(m - lobatto_beta * h),
        f(m - lobatto_x3 * h),
        f(m),
        f(m + lobatto_x3 * h),
        f(m + lobatto_beta * h),
        f(m + lobatto_x2 * h),
        f(m + lobatto_alpha * h),
        f(m + lobatto_x1 * h),
        f(b)
    };
    
    const T fa = y[0];
    const T fb = y[12];
    
    const T i2 = (h / 6.0) * (y[0] + y[12] + 5.0 * (y[4] + y[8]));
    const T i1 = (h / 1470.0) * (
            77.0 * (y[0] + y[12]) +
            432.0 * (y[2] + y[10]) +
            625.0 * (y[4] + y[8]) +
            672.0 * y[6]);
    const T is = h * (
        0.0158271919734802 * (y[0] + y[12]) + 
        0.0942738402188500 * (y[1] + y[11]) + 
        0.155071987336585  * (y[2] + y[10]) +
        0.188821573960182  * (y[3] + y[9]) + 
        0.199773405226859  * (y[4] + y[8]) +
        0.224926465333340  * (y[5] + y[7]) + 
        0.242611071901408  * y[6]);    
   
    const T erri1 = fabs(i1 - is);
    const T erri2 = fabs(i2 - is);
    
    const T err_is = get_error_is(p_tol, erri1, erri2, is, a, b);

    return adaptlobstp(f, a, b, fa, fb, err_is);
}

template <typename T, typename F>
T sinh_sigmoidal_adaptive_integral2d(double x0, double b,
    double tol, const F& f)
{
    auto mu_0 = 0.5 * (std::asinh((1.0 + x0) / b) + std::asinh((1.0 - x0) / b));
    auto eta_0 = 0.5 * (std::asinh((1.0 + x0) / b) - std::asinh((1.0 - x0) / b));
    return adaptive_integrate<T>([&](double s) {
            double x_hat = x0 + b * std::sinh(mu_0 * s - eta_0); 
            double jacobian = b * mu_0 * std::cosh(mu_0 * s - eta_0);
            return jacobian * f(x_hat);
        }, -1.0, 1.0, tol);
}

template <typename T, typename F>
T sinh_sigmoidal_single_tri(double b, double tol,
    const Vec<Vec<double,3>,3>& tri, const F& f) 
{
    double tri_area_jacobian = 2.0 * tri_area(tri);
    if (tri_area_jacobian == 0.0) {
        return zeros<T>::make();
    }
    constexpr double alpha = M_PI / 4.0;
    return tri_area_jacobian * adaptive_integrate<T>([&](double sigma) {
            double theta = (M_PI / 2.0) * sigma;
            double theta_jacobian = (M_PI / 2.0);
            double R_theta = std::sin(alpha) / std::sin(theta + alpha);
            double mu_1 = 0.5 * std::asinh(R_theta / b);
            return theta_jacobian * adaptive_integrate<T>([=](double s) {
                double r = b * std::sinh(mu_1 * (s + 1));
                double x = r * std::cos(theta);
                double y = r * std::sin(theta);
                auto pt = ref_to_real({x, y}, tri);

                double polar_to_cartes_jacobian = r;
                double sinh_jacobian = b * mu_1 * std::cosh(mu_1 * (s + 1));
                return polar_to_cartes_jacobian * sinh_jacobian * f({pt[0], pt[1]});
            }, -1.0, 1.0, tol);
        }, 0.0, 1.0, tol);
}

template <typename T, typename F>
T sinh_sigmoidal_adaptive_integral3d(double x0, double y0, double b,
    double tol, const F& f)
{
    Vec<double,3> pt0{0, 0, 0};
    Vec<double,3> pt1{1, 0, 0};
    Vec<double,3> pt2{0, 1, 0};
    Vec<double,3> singular_pt{x0, y0, 0};

    // upper left tri
    auto eval1 = sinh_sigmoidal_single_tri<T>(b, tol, {singular_pt, pt2, pt0}, f); 
    // upper right tri                                                       
    auto eval2 = sinh_sigmoidal_single_tri<T>(b, tol, {singular_pt, pt1, pt2}, f); 
    // lower                                                                 
    auto eval3 = sinh_sigmoidal_single_tri<T>(b, tol, {singular_pt, pt0, pt1}, f); 

    return eval1 + eval2 + eval3;
}

} //END NAMESPACE tbem
#endif
