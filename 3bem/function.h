#ifndef __rH1JH11LHLHAHAh_FUNCTION_H
#define __rH1JH11LHLHAHAh_FUNCTION_H

#include <vector>
#include <cassert>
#include <iostream>

namespace tbem {

/* typedef std::vector<double> Function; */

template <typename T>
struct InternalFunction {
    std::vector<T> _data;

    InternalFunction() {}

    InternalFunction(size_t n_elements):
        _data(n_elements)
    {}

    InternalFunction(size_t n_elements, const T& value):
        _data(n_elements, value)
    {}

    InternalFunction(const std::vector<T>& data):
        _data(data)
    {}

    InternalFunction(std::initializer_list<T> s):
        _data(s)
    {}

    void resize(size_t new_size) {
        _data.resize(new_size);
    }

    T& operator[] (size_t idx) {
        return _data[idx];
    }

    const T& operator[] (size_t idx) const {
        return _data[idx];
    }

    typename std::vector<T>::iterator begin() {
        return _data.begin();
    }

    typename std::vector<T>::const_iterator begin() const {
        return _data.begin();
    }

    typename std::vector<T>::iterator end() {
        return _data.end();
    }

    typename std::vector<T>::const_iterator end() const {
        return _data.end();
    }

    double* data() {
        return _data.data();
    }

    const double* data() const {
        return _data.data();
    }

    size_t size() const {
        return _data.size();
    }
};
typedef InternalFunction<double> Function;
typedef InternalFunction<Function> BlockFunction;

inline std::ostream& operator<<(std::ostream& os, const BlockFunction& a) {
    os << "(";
    for (size_t d = 0; d < a.size(); d++) {
        os << "[";
        for (size_t i = 0; i < a[d].size(); i++) {
            os << a[d][i] << ", ";
        }
        os << "], ";
    }
    os << ")";
    return os;
}

//TODO: Maybe get rid of this. Interesting macro experiment. What about python template?
#define __FunctionLoopBegin \
    assert(a.size() == b.size());\
    for (size_t d = 0; d < a.size(); d++) {\
        assert(a[d].size() == b[d].size());\
        for (size_t i = 0; i < a[d].size(); i++) {

#define __FunctionLoopEnd \
        } \
    }

template <typename T>
InternalFunction<T>& operator+=(InternalFunction<T>& a, const InternalFunction<T>& b) {
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
    return a;
}

inline BlockFunction& operator+=(BlockFunction& a, const BlockFunction& b) {
    assert(a.size() == b.size());
    for (size_t d = 0; d < a.size(); d++) {
        a[d] += b[d];
    }
    return a;
}


inline BlockFunction& operator-=(BlockFunction& a, const BlockFunction& b) {
    __FunctionLoopBegin;
    a[d][i] -= b[d][i];
    __FunctionLoopEnd;
    return a;
}

inline BlockFunction& operator*=(BlockFunction& a, const double b) {
    for (size_t d = 0; d < a.size(); d++) {\
        for (size_t i = 0; i < a[d].size(); i++) {
            a[d][i] *= b;
        }
    }
    return a;
}

inline BlockFunction operator*(const BlockFunction& a, double b) {
    BlockFunction out = a;
    out *= b;
    return out;
}

inline BlockFunction constant_function(size_t components, size_t dofs, double value) {
    return BlockFunction(components, Function(dofs, value));
}

inline BlockFunction operator-(const BlockFunction& a) {
    BlockFunction out(a.size());
    for (size_t d = 0; d < a.size(); d++) {
        out[d].resize(a[d].size());
        for (size_t i = 0; i < a[d].size(); i++) {
            out[d][i] = -a[d][i];  
        }
    }
    return out;
}

inline bool operator==(const BlockFunction& a, const BlockFunction& b) {
    __FunctionLoopBegin;
    if (a[d][i] != b[d][i]) {
        return false;
    }
    __FunctionLoopEnd;
    return true;
}

} // end namespace tbem

#endif
