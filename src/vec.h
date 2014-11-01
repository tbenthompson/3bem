#ifndef __ASDKJWEORIOQ_VEC_ARRAYS_H 
#define __ASDKJWEORIOQ_VEC_ARRAYS_H 
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

/* Some methods for efficient indexing of std::vector by a std::array<int,n>
 */
template <typename T>
std::array<T, 3> index3(const std::vector<T>& x, const std::array<int, 3>& indices) {
    std::array<T, 3> out;
    for (std::size_t i = 0; i < 3; i++) {
        out[i] = x[indices[i]];
    }
    return out;
}

/* [Note: Justification of overloading std::array operations]
 * Having a 3-vector class is useful because it reduces instances where
 * the same code has to be replicated 3 times for each element.
 *
 * Writing a new 3-vector class is not hard at all, but it can't hurt
 * to rely on an existing standard class. 
 *
 * This avoids large libraries like Eigen just to have a 3-vector class.
 *
 * This can later be extended to work on std::vector or arbitrary length 
 * std::array.
 *
 * C++11 templated typedefs allow this Vec3 "class" to be templated on float
 * vs double.
 *
 * It would be interesting to extend this to use expression templates for
 * lazy evaluation of results.
 */

// Cool! C++11 templated typedef
template <typename T>
using Vec3 = std::array<T,3>;

template <typename T>
void operator+=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
}

template <typename T>
void operator-=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] -= b[0]; a[1] -= b[1]; a[2] -= b[2];
}

template <typename T>
void operator*=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] *= b[0]; a[1] *= b[1]; a[2] *= b[2];
}

template <typename T>
void operator/=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] /= b[0]; a[1] /= b[1]; a[2] /= b[2];
}

template <typename T>
void operator*=(Vec3<T>& a, const T s) {
    a[0] *= s; a[1] *= s; a[2] *= s;
}

template <typename T>
void operator/=(Vec3<T>& a, const T s) {
    a[0] /= s; a[1] /= s; a[2] /= s;
}

template <typename T>
Vec3<T> operator+(const Vec3<T>& a, const Vec3<T>& b) {
    Vec3<T> res = a; res += b; return res;
}

template <typename T>
Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b) {
    Vec3<T> res = a; res -= b; return res;
}

template <typename T>
Vec3<T> operator*(const Vec3<T>& a, const Vec3<T>& b) {
    Vec3<T> res = a; res *= b; return res;
}

template <typename T>
Vec3<T> operator/(const Vec3<T>& a, const Vec3<T>& b) {
    Vec3<T> res = a; res /= b; return res;
}

template <typename T>
Vec3<T> operator*(const Vec3<T>& a, const T s) {
    Vec3<T> res = a; res *= s; return res;
}

template <typename T>
Vec3<T> operator*(const T s, const Vec3<T>& a) {
    Vec3<T> res = a; res *= s; return res;
}

template <typename T>
Vec3<T> operator/(const Vec3<T>& a, const T s) {
    Vec3<T> res = a; res /= s; return res;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vec3<T>& a) {
    os << "(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
    return os;
}

template <typename T>
Vec3<T> operator-(const Vec3<T>& v) {
    return {-v[0], -v[1], -v[2]};
}


template <typename T>
T sum(const Vec3<T>& a) {
    return a[0] + a[1] + a[2];
}

template <typename T>
Vec3<T> cross(const Vec3<T>& x, const Vec3<T>& y) {
    return {
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0]
    };
}

template <typename T>
T dot(const Vec3<T>& x, const Vec3<T>& y) {
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

template <typename T>
T hypot2(const Vec3<T>& v) {
    return dot(v, v);
}

template <typename T>
T hypot(const Vec3<T>& v) {
    return std::sqrt(hypot2(v));
}

template <typename T>
void normalize(Vec3<T>& v) {
    v /= hypot(v);
}

template <typename T>
Vec3<T> normalized(const Vec3<T>& v) {
    Vec3<T> res = v;
    normalize(res);
    return res;
}

template <typename T>
inline double dist2(const Vec3<T>& v0, 
                    const Vec3<T>& v1) {
    return hypot2(v1 - v0);
}

template <typename T>
Vec3<T> unit(const int k) {
    Vec3<T> e_k = {0.0, 0.0, 0.0};
    e_k[k] = 1.0;
    return e_k;
}

#endif
