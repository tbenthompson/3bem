#include "numbers.h"
#include <numeric>

namespace tbem {

/* equivalent to range(min, max) in python */
std::vector<int> integers(int min, int max) {
    std::vector<int> indices(max);
    std::iota(indices.begin(), indices.end(), min);
    return indices;
}

/* equivalent to range(0, max) in python */
std::vector<int> integers(int max) {
    return integers(0, max);
}

} //END namespace tbem
