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
using std::fabs;


template <typename T>
inline Vec3<bool> operator==(const Vec3<T>& t, const T& rhs) {
    return {t[0] == rhs, t[1] == rhs, t[2] == rhs};
}
template <typename T>
inline Vec3<bool> operator!=(const Vec3<T>& t, const T& rhs) {
    return {t[0] != rhs, t[1] != rhs, t[2] != rhs};
}
inline bool any(bool a) {return a;}
inline bool any(Vec3<bool> v) {return v[0] || v[1] || v[2];}
inline bool all(bool a) {return a;}
inline bool all(Vec3<bool> v) {return v[0] && v[1] && v[2];}
template <typename T>
inline T ones();
template <>
inline double ones<double>() {return 1.0;}
template <>
inline Vec3<double> ones<Vec3<double>>() {return {1.0, 1.0, 1.0};}
inline double max(double x) {return x;}
inline double max(Vec3<double> x) {return std::max(x[0], std::max(x[1], x[2]));}

const double alpha = std::sqrt(2./3.);
const double beta = 1./std::sqrt(5.);
const double x1 = .94288241569547971905635175843185720232;
const double x2 = .64185334234578130578123554132903188354;
const double x3 = .23638319966214988028222377349205292599;

template <typename T>
T adaptlobstp(const std::function<T(double)>& f, const double a, const double b, 
              const T& fa, const T& fb, const T& is )
{
    // std::cout << a << " " << b << std::endl;
    double m = (a + b) / 2.; 
    double h = (b - a) / 2.;

    double ah = alpha * h;
    double bh = beta * h;
    double mll = m - ah;
    double ml = m - bh;
    double mr = m + bh;
    double mrr = m + ah;

    T fmll, fml, fm, fmr, fmrr;
    fmll = f(mll); fml = f(ml); fm = f(m); fmr = f(mr); fmrr = f(mrr);

    T i1, i2;
    i2 = (h / 6.) * (fa + fb + 5.0 * (fml + fmr));    
    i1 = (h / 1470.) * (
            77.0 * (fa + fb) + 
            432.0 * (fmll + fmrr) +
            625.0 * (fml + fmr) +
            672.0 * fm);

    if (all(is + (i1 - i2) == is)){
        return i1;
    } else if (mll <= a or b <= mrr) {
        std::cout << "YIKES!" << std::endl;
        return i1;
    }
    else
    {
        return adaptlobstp(f, a, mll, fa, fmll, is)
             + adaptlobstp(f, mll, ml, fmll, fml, is)
             + adaptlobstp(f, ml, m, fml, fm, is)
             + adaptlobstp(f, m, mr, fm, fmr, is)
             + adaptlobstp(f, mr, mrr, fmr, fmrr, is)
             + adaptlobstp(f, mrr, b, fmrr, fb, is);
    }
}

//TODO: Refactor the shit out of this! Super ugly.
template <typename T>
T adaptive_integrate(std::function<T(double)> f, double a, double b, double p_tol)
{
    // std::cout << "START" << std::endl;

    double eps = std::numeric_limits<double>::epsilon();
    // double tol = (p_tol < eps) ?  eps : p_tol;

    double m = (a + b) / 2.; 
    double h = (b - a) / 2.;

    T y[13] = {
        f(a),f(m-x1*h),f(m-alpha*h),f(m-x2*h),f(m-beta*h),
        f(m-x3*h),f(m), f(m+x3*h),f(m+beta*h),f(m+x2*h),
        f(m+alpha*h),f(m+x1*h),f(b)
    };
    
    T fa = y[0];
    T fb = y[12];
    
    T i2 = (h / 6.0) * (y[0] + y[12] + 5.0 *(y[4] + y[8]));
    T i1 = (h / 1470.0) * (
            77.0 * (y[0] + y[12]) +
            432.0 * (y[2] + y[10]) +
            625.0 * (y[4] + y[8]) +
            672.0 * y[6]);
    T is = h * (
        0.0158271919734802 * (y[0] + y[12]) + 
        0.0942738402188500 * (y[1] + y[11]) + 
        0.155071987336585  * (y[2] + y[10]) +
        0.188821573960182  * (y[3] + y[9]) + 
        0.199773405226859  * (y[4] + y[8]) +
        0.224926465333340  * (y[5] + y[7]) + 
        0.242611071901408  * y[6]);    
   
    T erri1 = fabs(i1 - is);
    T erri2 = fabs(i2 - is);
    double R = (all(erri2 != 0.0)) ? max(erri1 / erri2) : 1.0;
    std::cout << R << std::endl;
    
    double tol = (R > 0. and R < 1.) ? p_tol / R : p_tol; 
    is = fabs(is) * tol / eps;
    if (any(is == 0.0)) {
        is = ones<T>() * (b - a);
    }

    return adaptlobstp(f, a, b, fa, fb, is);
}

#endif
