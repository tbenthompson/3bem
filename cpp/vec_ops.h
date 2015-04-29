#ifndef __QWEQWEQWEQEWQWEQWE_VEC_OPS_H
#define __QWEQWEQWEQEWQWEQWE_VEC_OPS_H

#include "vec.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "numbers.h"

namespace tbem {

template <typename T>
void operator+=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
}
template <typename T>
void operator+=(Vec2<T>& a, const Vec2<T>& b) {
    a[0] += b[0]; a[1] += b[1];
}
template <typename T>
void operator+=(Vec1<T>& a, const Vec1<T>& b) {
    a[0] += b[0];
}

template <typename T, typename F>
void operator+=(Vec3<T>& a, const F& s) {
    a[0] += s; a[1] += s; a[2] += s;
}
template <typename T, typename F>
void operator+=(Vec2<T>& a, const F& s) {
    a[0] += s; a[1] += s;
}
template <typename T, typename F>
void operator+=(Vec1<T>& a, const F& s) {
    a[0] += s; 
}

template <typename T>
void operator-=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] -= b[0]; a[1] -= b[1]; a[2] -= b[2];
}
template <typename T>
void operator-=(Vec2<T>& a, const Vec2<T>& b) {
    a[0] -= b[0]; a[1] -= b[1];
}
template <typename T>
void operator-=(Vec1<T>& a, const Vec1<T>& b) {
    a[0] -= b[0];
}

template <typename T, typename F>
void operator-=(Vec3<T>& a, const F& s) {
    a[0] -= s; a[1] -= s; a[2] -= s;
}
template <typename T, typename F>
void operator-=(Vec2<T>& a, const F& s) {
    a[0] -= s; a[1] -= s;
}
template <typename T, typename F>
void operator-=(Vec1<T>& a, const F& s) {
    a[0] -= s; 
}

template <typename T>
void operator*=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] *= b[0]; a[1] *= b[1]; a[2] *= b[2];
}
template <typename T>
void operator*=(Vec2<T>& a, const Vec2<T>& b) {
    a[0] *= b[0]; a[1] *= b[1];
}
template <typename T>
void operator*=(Vec1<T>& a, const Vec1<T>& b) {
    a[0] *= b[0]; 
}

template <typename T>
void operator/=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] /= b[0]; a[1] /= b[1]; a[2] /= b[2];
}
template <typename T>
void operator/=(Vec2<T>& a, const Vec2<T>& b) {
    a[0] /= b[0]; a[1] /= b[1]; 
}
template <typename T>
void operator/=(Vec1<T>& a, const Vec1<T>& b) {
    a[0] /= b[0]; 
}

template <typename T, typename F>
void operator*=(Vec3<T>& a, const F& s) {
    a[0] *= s; a[1] *= s; a[2] *= s;
}
template <typename T, typename F>
void operator*=(Vec2<T>& a, const F& s) {
    a[0] *= s; a[1] *= s;
}
template <typename T, typename F>
void operator*=(Vec1<T>& a, const F& s) {
    a[0] *= s; 
}

template <typename T, size_t dim>
void operator/=(Vec<T,dim>& a, const T& s) {
    double inv_s = 1 / s;
    a *= inv_s;
}

template <typename T, size_t dim>
Vec<T,dim> operator+(const Vec<T,dim>& a, const Vec<T,dim>& b) {
    Vec<T,dim> res = a; res += b; return res;
}

template <typename T, typename F, size_t dim>
Vec<T,dim> operator+(const Vec<T,dim>& a, const F& s) {
    Vec<T,dim> res = a; res += s; return res;
}

template <typename T, typename F, size_t dim>
Vec<T,dim> operator+(const F& s, const Vec<T,dim>& a) {
    Vec<T,dim> res = a; res += s; return res;
}

template <typename T, size_t dim>
Vec<T,dim> operator-(const Vec<T,dim>& a, const Vec<T,dim>& b) {
    Vec<T,dim> res = a; res -= b; return res;
}

template <typename T, typename F, size_t dim>
Vec<T,dim> operator-(const Vec<T,dim>& a, const F& s) {
    Vec<T,dim> res = a; res -= s; return res;
}

template <typename T, typename F, size_t dim>
Vec<T,dim> operator-(const F& s, const Vec<T,dim>& a) {
    Vec<T,dim> res = a; res -= s; return res;
}

template <typename T, size_t dim>
Vec<T,dim> operator*(const Vec<T,dim>& a, const Vec<T,dim>& b) {
    Vec<T,dim> res = a; res *= b; return res;
}

template <typename T, size_t dim>
Vec<T,dim> operator/(const Vec<T,dim>& a, const Vec<T,dim>& b) {
    Vec<T,dim> res = a; res /= b; return res;
}

template <typename T, typename F, size_t dim>
Vec<T,dim> operator*(const Vec<T,dim>& a, const F& s) {
    Vec<T,dim> res = a; res *= s; return res;
}

template <typename T, typename F, size_t dim>
Vec<T,dim> operator*(const F& s, const Vec<T,dim>& a) {
    Vec<T,dim> res = a; res *= s; return res;
}

template <typename T, size_t dim>
Vec<T,dim> operator/(const Vec<T,dim>& a, const T& s) {
    Vec<T,dim> res = a; res /= s; return res;
}

template <typename T>
Vec<T,1> operator/(const T& s, const Vec<T,1>& a) {
    return {s / a[0]};
}

template <typename T>
Vec<T,2> operator/(const T& s, const Vec<T,2>& a) {
    return {s / a[0], s / a[1]};
}

template <typename T>
Vec<T,3> operator/(const T& s, const Vec<T,3>& a) {
    return {s / a[0], s / a[1], s / a[2]};
}


template <typename T, size_t dim>
std::ostream& operator<<(std::ostream& os, const Vec<T,dim>& a) {
    os << "(";
    if (dim > 0) {
        for (std::size_t i = 0; i < dim - 1; i++) {
            os << a[i] << ", ";
        }
        os << a[dim - 1];
    }
    os << ")";
    return os;
}

template <typename T, size_t dim>
Vec<T,dim> operator-(const Vec<T,dim>& v) {
    return zeros<Vec<T,dim>>::make() - v;
}


template <typename T>
T sum(const Vec3<T>& a) {
    return a[0] + a[1] + a[2];
}

template <typename T>
T sum(const Vec2<T>& a) {
    return a[0] + a[1];
}

using std::fabs;

template <typename T>
Vec1<T> fabs(const Vec1<T>& v) {
    return {fabs(v[0])};
}
template <typename T>
Vec2<T> fabs(const Vec2<T>& v) {
    return {fabs(v[0]), fabs(v[1])};
}
template <typename T>
Vec3<T> fabs(const Vec3<T>& v) {
    return {fabs(v[0]), fabs(v[1]), fabs(v[2])};
}

template <typename T>
inline Vec3<bool> operator==(const Vec3<T>& t, const T& rhs) {
    return {t[0] == rhs, t[1] == rhs, t[2] == rhs};
}
template <typename T>
inline Vec2<bool> operator==(const Vec2<T>& t, const T& rhs) {
    return {t[0] == rhs, t[1] == rhs};
}
template <typename T>
inline Vec1<bool> operator==(const Vec1<T>& t, const T& rhs) {
    return {t[0] == rhs};
}

template <typename T>
inline Vec3<bool> operator!=(const Vec3<T>& t, const T& rhs) {
    return {t[0] != rhs, t[1] != rhs, t[2] != rhs};
}
template <typename T>
inline Vec2<bool> operator!=(const Vec2<T>& t, const T& rhs) {
    return {t[0] != rhs, t[1] != rhs};
}
template <typename T>
inline Vec1<bool> operator!=(const Vec1<T>& t, const T& rhs) {
    return {t[0] != rhs};
}

template <typename T>
inline Vec3<bool> operator<(const Vec3<T>& t, const T& rhs) {
    return {t[0] < rhs, t[1] < rhs, t[2] < rhs};
}
template <typename T>
inline Vec2<bool> operator<(const Vec2<T>& t, const T& rhs) {
    return {t[0] < rhs, t[1] < rhs};
}
template <typename T>
inline Vec1<bool> operator<(const Vec1<T>& t, const T& rhs) {
    return {t[0] < rhs};
}

template <typename T>
inline Vec3<bool> operator>(const Vec3<T>& t, const T& rhs) {
    return {t[0] > rhs, t[1] > rhs, t[2] > rhs};
}
template <typename T>
inline Vec2<bool> operator>(const Vec2<T>& t, const T& rhs) {
    return {t[0] > rhs, t[1] > rhs};
}
template <typename T>
inline Vec1<bool> operator>(const Vec1<T>& t, const T& rhs) {
    return {t[0] > rhs};
}

inline double max(double x) {return x;}
inline double max(Vec2<double> x) {return std::max(x[0], x[1]);}
inline double max(Vec3<double> x) {return std::max(x[0], std::max(x[1], x[2]));}

inline double min(double x) {return x;}
inline double min(Vec3<double> x) {return std::min(x[0], std::min(x[1], x[2]));}

inline bool any(bool a) {return a;}
inline bool any(Vec3<bool> v) {return v[0] || v[1] || v[2];}
inline bool any(Vec2<bool> v) {return v[0] || v[1];}

inline bool all(bool a) {return a;}
inline bool all(Vec3<bool> v) {return v[0] && v[1] && v[2];}
inline bool all(Vec2<bool> v) {return v[0] && v[1];}

template <typename F>
struct constant<Vec1<F>> {
    static Vec1<F> make(double val) { 
        return {constant<F>::make(val)};
    }
};

template <typename F>
struct constant<Vec2<F>> {
    static Vec2<F> make(double val) { 
        return {
            constant<F>::make(val), constant<F>::make(val)
        };
    }
};

template <typename F>
struct constant<Vec3<F>> {
    static Vec3<F> make(double val) { 
        return {
            constant<F>::make(val), constant<F>::make(val), constant<F>::make(val) 
        };
    }
};

} // end namespace tbem

#endif
