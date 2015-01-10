#ifndef __RRRRJJJJJJZZZZZZEEEAKSDJLKJ_CONSTANTS_H
#define __RRRRJJJJJJZZZZZZEEEAKSDJLKJ_CONSTANTS_H
#include <vector>

namespace tbem {

std::vector<int> integers(int min, int max);
std::vector<int> integers(int max);

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
