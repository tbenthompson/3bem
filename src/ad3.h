#ifndef __ASDQWEQHWJKH_AD_H
#define __ASDQWEQHWJKH_AD_H
#include <array>
#include <iostream>
#include <functional>
#include <cmath>

/* Tower-based automatic differentiation for univariate functions using a 
 * recursive template-based approach. This partially follows the exposition
 * in
 * http://conal.net/papers/beautiful-differentiation/beautiful-differentiation-long.pdf
 */

namespace ad {

/* The recursive tower derivatives class.
 * T = float, double, etc -- anything that implements standard operations
 * M = the maximum-degree for this derivatves tower.
 * Nd = the degree of this object.
 */
template<typename T, int M, int Nd = 0>
struct D {
    typedef D<T,M,Nd> type;
    typedef D<T,M,Nd+1> d_type;

    //Allows implicit conversion from the underlying type to a constant
    //tower derivatives
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

    void get_all(std::array<T,M+1>& out) const {
        out[Nd] = v;
        d.get_all(out);
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

    static type chain(const type& a,
                      const std::function<T(T)>& f,
                      const std::function<d_type(const d_type&)>& fp) {
        return type{f(a.v), a.d * fp(shorter(a))};
    }
};


/* The base case for the recursive data structure, when
 * degree = max degree (M == Nd)
 */
template<typename T, int M>
struct D<T,M,M> {
    typedef D<T,M,M> type;

    const T v;

    //Allows implicit conversion from the underlying type to a constant
    //tower derivatives
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

    void get_all(std::array<T,M+1>& out) const {
        out[M] = v;
    }

    static type D_from_const(T x) {
        return type{x};
    }

    static type D_from_id(T x) {
        return D_from_const(x);
    }

    static type chain(const type& a,
                      const std::function<T(T)>& f) {
        return type{f(a.v)};
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

template <typename T, int M, int Nd> 
struct helper {
    static D<T,M,Nd> help(const D<T,M,Nd>& a) {
        auto f = exp<T,M,Nd+1>;
        return D<T,M,Nd>::chain(a, (double(*)(double))&std::exp, &f);
    }
};

template <typename T, int M> 
struct helper<T,M,M> {
    static D<T,M,M> help(const D<T,M,M>& a) {
        return D<T,M,M>::chain(a, (double(*)(double))&std::exp);
    }
};

template <typename T, int M, int Nd>
D<T,M,Nd> exp(const D<T,M,Nd>& a) {
    return helper<T,M,Nd>::help(a);
}

// A convenience typedef for a derivative tower for doubles.
template <int M>
using D_d = D<double,M,0>;
}
#endif
