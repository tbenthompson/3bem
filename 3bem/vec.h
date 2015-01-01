#ifndef __ASDKJWEORIOQ_VEC_ARRAYS_H 
#define __ASDKJWEORIOQ_VEC_ARRAYS_H 
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace tbem {

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
 *
 * Extending this to use recursion to do 1D,2D,3D all in one.
 */

// Cool! C++11 templated typedef
template <typename T, unsigned long dim>
using Vec = std::array<T,dim>;
template <typename T>
using Vec3 = Vec<T,3>;
template <typename T>
using Vec2 = Vec<T,2>;

template <typename T>
void operator+=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
}
template <typename T>
void operator+=(Vec2<T>& a, const Vec2<T>& b) {
    a[0] += b[0]; a[1] += b[1];
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
void operator*=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] *= b[0]; a[1] *= b[1]; a[2] *= b[2];
}
template <typename T>
void operator*=(Vec2<T>& a, const Vec2<T>& b) {
    a[0] *= b[0]; a[1] *= b[1];
}

template <typename T>
void operator/=(Vec3<T>& a, const Vec3<T>& b) {
    a[0] /= b[0]; a[1] /= b[1]; a[2] /= b[2];
}
template <typename T>
void operator/=(Vec2<T>& a, const Vec2<T>& b) {
    a[0] /= b[0]; a[1] /= b[1]; 
}

template <typename T, typename F>
void operator*=(Vec3<T>& a, const F& s) {
    a[0] *= s; a[1] *= s; a[2] *= s;
}
template <typename T, typename F>
void operator*=(Vec2<T>& a, const F& s) {
    a[0] *= s; a[1] *= s;
}

template <typename T, unsigned long dim>
void operator/=(Vec<T,dim>& a, const T& s) {
    double inv_s = 1 / s;
    a *= inv_s;
}

template <typename T, unsigned long dim>
Vec<T,dim> operator+(const Vec<T,dim>& a, const Vec<T,dim>& b) {
    Vec<T,dim> res = a; res += b; return res;
}

template <typename T, unsigned long dim>
Vec<T,dim> operator-(const Vec<T,dim>& a, const Vec<T,dim>& b) {
    Vec<T,dim> res = a; res -= b; return res;
}

template <typename T, unsigned long dim>
Vec<T,dim> operator*(const Vec<T,dim>& a, const Vec<T,dim>& b) {
    Vec<T,dim> res = a; res *= b; return res;
}

template <typename T, unsigned long dim>
Vec<T,dim> operator/(const Vec<T,dim>& a, const Vec<T,dim>& b) {
    Vec<T,dim> res = a; res /= b; return res;
}

template <typename T, typename F, unsigned long dim>
Vec<T,dim> operator*(const Vec<T,dim>& a, const F& s) {
    Vec<T,dim> res = a; res *= s; return res;
}

template <typename T, typename F, unsigned long dim>
Vec<T,dim> operator*(const F& s, const Vec<T,dim>& a) {
    Vec<T,dim> res = a; res *= s; return res;
}

template <typename T, unsigned long dim>
Vec<T,dim> operator/(const Vec<T,dim>& a, const T& s) {
    Vec<T,dim> res = a; res /= s; return res;
}

template <typename T, unsigned long dim>
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

template <typename T>
Vec3<T> operator-(const Vec3<T>& v) {
    return {-v[0], -v[1], -v[2]};
}


template <typename T>
T sum(const Vec3<T>& a) {
    return a[0] + a[1] + a[2];
}

template <typename T>
T sum(const Vec2<T>& a) {
    return a[0] + a[1];
}

template <typename T>
Vec<T,3> outer_product(const Vec<double,3>& a, const T& b) {
    return {b * a[0], b * a[1], b * a[2]};
}

template <typename T>
Vec<T,2> outer_product(const Vec<double,2>& a, const T& b) {
    return {b * a[0], b * a[1]};
}

template <typename T>
Vec3<T> cross(const Vec3<T>& x, const Vec3<T>& y) {
    return {
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0]
    };
}

inline double dot_product(const double& x, const double& y) {
    return x * y;
}

template <typename T>
T dot_product(const Vec<double,3>& x, const Vec<T,3>& y) {
    return sum(Vec<T,3>{{x[0] * y[0], x[1] * y[1], x[2] * y[2]}});
}

template <typename T>
T dot_product(const Vec<double,2>& x, const Vec<T,2>& y) {
    return sum(Vec<T,2>{{x[0] * y[0], x[1] * y[1]}});
}

template <typename T, unsigned long dim>
T hypot2(const Vec<T,dim>& v) {
    return dot_product(v, v);
}

template <typename T, unsigned long dim>
T hypot(const Vec<T,dim>& v) {
    return std::sqrt(hypot2(v));
}

template <typename T, unsigned long dim>
void normalize(Vec<T,dim>& v) {
    v /= hypot(v);
}

template <typename T, unsigned long dim>
Vec<T,dim> normalized(const Vec<T,dim>& v) {
    Vec<T,dim> res = v;
    normalize(res);
    return res;
}

template <typename T, unsigned long dim>
inline T dist2(const Vec<T,dim>& v0, const Vec<T,dim>& v1) {
    return hypot2(v1 - v0);
}

template <typename T, unsigned long dim>
inline T dist(const Vec<T,dim>& v0, const Vec<T,dim>& v1) {
    return hypot(v1 - v0);
}

template <typename T, unsigned long dim>
inline Vec<T,dim> unit(const int k);

template <>
inline Vec3<double> unit<double,3>(const int k) {
    Vec3<double> e_k = {0.0, 0.0, 0.0};
    e_k[k] = 1.0;
    return e_k;
}

template <>
inline Vec2<double> unit<double,2>(const int k) {
    Vec2<double> e_k = {0.0, 0.0};
    e_k[k] = 1.0;
    return e_k;
}

using std::fabs;
template <typename T>
Vec3<T> fabs(const Vec3<T>& v) {
    return {fabs(v[0]), fabs(v[1]), fabs(v[2])};
}
template <typename T>
Vec2<T> fabs(const Vec2<T>& v) {
    return {fabs(v[0]), fabs(v[1])};
}

inline Vec3<double> 
unscaled_normal(const std::array<Vec3<double>,3>& corners) {
    return cross(corners[2] - corners[0], corners[2] - corners[1]);
}

inline Vec2<double> 
unscaled_normal(const std::array<Vec2<double>,2>& corners) {
    return {
        -(corners[1][1] - corners[0][1]),
        corners[1][0] - corners[0][0]
    };
}

inline Vec3<double> 
tri_normal(const std::array<Vec3<double>,3>& corners) {
    auto unscaled = unscaled_normal(corners);
    return normalized(unscaled);
}

inline double tri_area(const Vec3<double>& unscaled_normal) {
    return 0.5 * hypot(unscaled_normal);
}

inline double tri_area(const std::array<Vec3<double>,3>& corners) {
    return tri_area(unscaled_normal(corners));
}

enum Side {FRONT, INTERSECT, BEHIND};

/* Determine which side of the plane/line defined by triangle/segment 
 * the provided point is on.
 */
template <int dim>
Side which_side_point(const std::array<Vec<double,dim>,dim>& face,
                const Vec<double,dim>& pt) {
    auto normal = unscaled_normal(face);
    double dot_val = dot_product(pt - face[0], normal);
    if (dot_val > 0) { return FRONT; }
    else if (dot_val < 0) { return BEHIND; }
    else { return INTERSECT; }
}

/* Returns the side of a plane that a triangle/segment is on. */
template <int dim>
Side facet_side(std::array<Side,dim> s);

template <>
inline Side facet_side<2>(std::array<Side,2> s) {
    if (s[0] == s[1]) { return s[0]; } 
    else if(s[0] == INTERSECT) { return s[1]; }
    else if(s[1] == INTERSECT) { return s[0]; }
    else { return INTERSECT; }
}

template <>
inline Side facet_side<3>(std::array<Side,3> s) {
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
template <int dim>
Side which_side_facet(const std::array<Vec<double,dim>,dim>& plane,
                const std::array<Vec<double,dim>,dim>& face) {
    std::array<Side,dim> sides;
    for (int d = 0; d < dim; d++) {
        sides[d] = which_side_point<dim>(plane, face[d]);
    }
    return facet_side<dim>(sides);
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
inline Vec3<bool> operator!=(const Vec3<T>& t, const T& rhs) {
    return {t[0] != rhs, t[1] != rhs, t[2] != rhs};
}
template <typename T>
inline Vec2<bool> operator!=(const Vec2<T>& t, const T& rhs) {
    return {t[0] != rhs, t[1] != rhs};
}

template <typename T>
inline Vec3<bool> operator<(const Vec3<T>& t, const T& rhs) {
    return {t[0] < rhs, t[1] < rhs, t[2] < rhs};
}
template <typename T>
inline Vec2<bool> operator<(const Vec2<T>& t, const T& rhs) {
    return {t[0] < rhs, t[1] < rhs};
}

inline double max(double x) {return x;}
inline double max(Vec3<double> x) {return std::max(x[0], std::max(x[1], x[2]));}

inline double min(double x) {return x;}
inline double min(Vec3<double> x) {return std::min(x[0], std::min(x[1], x[2]));}

inline bool any(bool a) {return a;}
inline bool any(Vec3<bool> v) {return v[0] || v[1] || v[2];}

inline bool all(bool a) {return a;}
inline bool all(Vec3<bool> v) {return v[0] && v[1] && v[2];}
inline bool all(Vec2<bool> v) {return v[0] && v[1];}

template <typename T, typename F = void>
struct constant;

template <>
struct constant<double> {
    static double make(double val) { return val; }
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

template <typename T, typename F = void>
struct ones {
    static T make() {
        return constant<T,F>::make(1.0);
    }
};

template <typename T, typename F = void>
struct zeros {
    static T make() {
        return constant<T,F>::make(0.0);
    }
};

template <typename T, typename F>
std::vector<T> reinterpret_vector(const std::vector<F>& A) {
    std::vector<T> out((T*)A.data(), (T*)(A.data() + A.size()));
    return out;
}

} // END namespace tbem

#endif
