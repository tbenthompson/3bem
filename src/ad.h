#ifndef __ASDQWEQHWJKH_TAYLOR_H
#define __ASDQWEQHWJKH_TAYLOR_H
#include <array>
#include <iostream>

/* Tower-based automatic differentiation for univariate functions using a 
 * recursive template-based approach. This partially follows the exposition
 * in
 * http://conal.net/papers/beautiful-differentiation/beautiful-differentiation-long.pdf
 */

namespace tower_ad {

/* The recursive tower series class.
 * T = float, double, etc -- anything that implements standard operations
 * M = the maximum-degree for this tower series tower.
 * Nd = the degree of this object.
 */
template<typename T, int M, int Nd = 0>
struct D {
    typedef D<T,M,Nd> type;
    typedef D<T,M,Nd+1> d_type;

    //Allows implicit conversion from the underlying type to a constant
    //tower series
    D(T x): v(x), d(d_type::D_from_const(0.0)) {}

    // Default constructor allows initializer lists, which are otherwise
    // disallowed by the conversion constructor.
    D(T p_v, d_type p_d): v(p_v), d(p_d) {}

    const T v;
    const d_type d;

    type operator*(const type& b) const {
        // std::cout << v << " " << d.v << " " << b.v << " " << b.d.v << std::endl;
        return type{v * b.v, b.d * shorter(*this) + d * shorter(b)};
    }

    type operator+(const type& b) const {
        return type{v + b.v, d + b.d};
    }

    type operator-() const {
        return type{-v, -d};
    }

    type operator-(const type& b) const {
        return (*this + (-b));
    }

    friend std::ostream& operator<<(std::ostream& out, 
                                    const type& d) {
        out << d.v << " " << d.d;
        return out;
    }

    static type D_from_const(T x) {
        return type{
            x, 
            d_type::D_from_const(0)
        };
    }

    static type D_from_id(T x) {
        return type{
            x, 
            d_type::D_from_const(1)
        };
    }
};

template <typename T, int M, int Nd>
typename std::enable_if<M==Nd+1,D<T,M,Nd+1>>::type shorter(const D<T,M,Nd>& b) {
    return D<T,M,Nd+1>{b.v};
}

template <typename T, int M, int Nd>
typename std::enable_if<M!=Nd+1,D<T,M,Nd+1>>::type shorter(const D<T,M,Nd>& b) {
    return D<T,M,Nd+1>{b.v, shorter(b.d)};
}

/* The base case for the recursive data structure, when
 * degree = max degree (M == Nd)
 */
template<typename T, int M>
struct D<T,M,M> {
    typedef D<T,M,M> type;

    const T v;

    //Allows implicit conversion from the underlying type to a constant
    //tower series
    D(T x): v(x) {}

    type operator*(const type& b) const {
        return type{v * b.v};
    }

    type operator-() const {
        return type{-v};
    }

    type operator+(const type& b) const {
        return type{v + b.v};
    }

    friend std::ostream& operator<<(std::ostream& out, 
                                    const type& d) {
        out << d.v;
        return out;
    }

    static type D_from_const(T x) {
        return type{x};
    }

    static type D_from_id(T x) {
        return D_from_const(x);
    }
};

// A convenience typedef for a derivative tower for doubles.
template <int M>
using D_d = D<double,M,0>;

}
#endif
