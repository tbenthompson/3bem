#ifndef QWEKJHWELR_AD2_H
#define QWEKJHWELR_AD2_H
#include <array>
#include <functional>
#include <cmath>
#include <iostream>
#include <cassert>

template <typename T, int M>
std::array<T,M> _constant(const T& x) {
    std::array<T,M> v{};
    v[0] = x;
    return v;
}

template <typename T, int M>
struct D {
    static_assert(M >= 1, "M must be >= 1");

    typedef D<T,M> type;
    typedef D<T,M-1> d_type;

    D(const T& x): v(_constant<T,M>(x)), degree(M) {}
    D(const std::array<T,M>& vals): v(vals), degree(M) {}
    D(const std::array<T,M>& vals, int degree): v(vals), degree(degree) {}

    const std::array<T, M> v;
    int degree;

    type init() const {
        std::array<T,M> out{};
        for (int i = 0; i < degree - 1; i++) {
            out[i] = v[i];
        }
        return type{out, degree - 1};
    }

    type tail() const {
        std::array<T,M> out{};
        for (int i = 0; i < degree - 1; i++) {
            out[i] = v[i + 1];
        }
        return type{out, degree - 1};
    }

    type chain(std::function<double(double)> f,
               std::function<type(const type&)> fp) const {
        std::array<T,M> res{};    
        res[0] = f(v[0]);
        if (degree == 0) {
            return type{res, 0};
        }
        auto c = fp(init()) * tail();
        for (int i = 0; i < c.degree; i++) {
            res[i + 1] += c.v[i];
        }
        return type{res, degree};
    }

    type operator*(const type& b) const {
        assert(b.degree == degree);
        std::array<T,M> res;
        res[0] = v[0] * b.v[0];
        if (degree == 0) {
            return type{res, 0};
        }
        auto pr = init() * b.tail() + tail() * b.init();
        for(int i = 0; i < pr.degree; i++) {
            res[i + 1] = pr.v[i];
        }
        return type{res, degree};
    }

    type operator+(const type& b) const {
        assert(b.degree == degree);
        std::array<T,M> vals = b.v;
        for (int i = 0; i < degree; i++) {
            vals[i] += v[i];
        }
        return type{vals};
    }

    type operator-() const {
        std::array<T,M> vals{};
        for (int i = 0; i < degree; i++) {
            vals[i] = -v[i];
        }
        return type{vals};
    }

    type operator-(const type& b) const {
        return (*this) + (-b);
    }

    friend std::ostream& operator<<(std::ostream& os, const type& a) {
        os << "(";
        for (int i = 0; i < M - 1; i++) {
            os << a.v[i] << ", ";
        }
        os << a.v[M-1] << ")";
        return os;
    }
};

template <typename T, int M>
D<T,M> var(const T& val) {
    std::array<T,M> v{};
    v[0] = val;
    v[1] = 1.0;
    return D<T,M>{v};
}

template <typename T, int M>
D<T,M> constant(const T& x) {
    return D<T,M>{_constant<T,M>(x)};
}

template <typename T, int M>
D<T,M> exp(const D<T,M>& a) {
    return a.chain((double(*)(double))&std::exp, &exp<T,M>);
}

template <typename T, int M>
D<T,M> sin(const D<T,M>& a);

template <typename T, int M>
D<T,M> cos(const D<T,M>& a) {
    return a.chain((double(*)(double))&std::cos, [](const D<T,M>& x) {
        return -sin<T,M>(x);
    });
}

template <typename T, int M>
D<T,M> sin(const D<T,M>& a) {
    return a.chain((double(*)(double))&std::sin, &cos<T,M>);
}

template <typename T>
T recip(const T& x) {return 1.0 / x;}

template <typename T, int M>
D<T,M> recip(const D<T,M>& a) {
    return a.chain(recip<T>, [](const D<T,M>& x) {
        auto r = recip(x); 
        return (-r) * r;
    });
}

template <typename T, int M>
D<T,M> log(const D<T,M>& a) {
    return a.chain((double(*)(double))&std::log, &recip<T,M>);
}

template <typename T, int M>
D<T,M> sqrt(const D<T,M>& a) {
    return a.chain((double(*)(double))&std::sqrt, [](const D<T,M>& x) {
        return recip(sqrt(x) * 2.0);
    });
}

template <typename T, int M>
D<T,M> asin(const D<T,M>& a) {
    return a.chain((double(*)(double))&std::asin, [](const D<T,M>& x) {
        return recip(sqrt(-(x * x - 1.0)));
    });
}


template <int M>
using D_d = D<double, M>;

#endif
