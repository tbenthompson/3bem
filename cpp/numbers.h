#ifndef TBEMRRRRJJJJJJZZZZZZEEEAKSDJLKJ_CONSTANTS_H
#define TBEMRRRRJJJJJJZZZZZZEEEAKSDJLKJ_CONSTANTS_H
#include <vector>
#include <cstddef>
#include <numeric>

namespace tbem {

/* Equivalent of python's range, returns all n for which min <= n < max
 */
template <typename T = int>
std::vector<T> range(T min, T max) {
    std::vector<T> indices(max - min);
    std::iota(indices.begin(), indices.end(), min);
    return indices;
}

/* range(0, max)
 */
template <typename T = int>
std::vector<T> range(T max) {
    return range<T>(0, max);
}

std::vector<double> linspace(const double& a, const double& b, size_t count);

template <typename T, typename F = void>
struct constant;

template <>
struct constant<double> {
    static double make(double val) { return val; }
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

}//END namespace tbem

#endif
