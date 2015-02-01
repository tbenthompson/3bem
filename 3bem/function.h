#ifndef __rH1JH11LHLHAHAh_FUNCTION_H
#define __rH1JH11LHLHAHAh_FUNCTION_H

#include <vector>
#include <cassert>
#include <iostream>

namespace tbem {

typedef std::vector<double> Function;
typedef std::vector<std::vector<double>> BlockFunction;

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

inline BlockFunction& operator+=(BlockFunction& a, const BlockFunction& b) {
    __FunctionLoopBegin;
    a[d][i] += b[d][i];
    __FunctionLoopEnd;
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
    return BlockFunction(components, std::vector<double>(dofs, value));
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

inline bool operator==(BlockFunction& a, const BlockFunction& b) {
    __FunctionLoopBegin;
    if (a[d][i] != b[d][i]) {
        return false;
    }
    __FunctionLoopEnd;
    return true;
}

struct ConcatenatedFunction 
{
    const size_t components;
    const std::vector<double> data;
    const std::vector<size_t> component_lengths;
};


ConcatenatedFunction concatenate(const BlockFunction& fncs); 

BlockFunction expand(const ConcatenatedFunction& block_fnc,
    const std::vector<double>& replacement_data);

BlockFunction expand(const ConcatenatedFunction& block_fnc); 

} // end namespace tbem

#endif
