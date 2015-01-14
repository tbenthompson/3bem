#ifndef __RRRRJJJJJJZZZZZZEEEAKSDJLKJ_CONSTANTS_H
#define __RRRRJJJJJJZZZZZZEEEAKSDJLKJ_CONSTANTS_H
#include <vector>
#include <cstddef>

namespace tbem {

/* Equivalent of python's range, returns all n for which min <= n < max
 */
std::vector<int> range(int min, int max);

/* range(0, max)
 */
std::vector<int> range(int max);

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
