#include "numbers.h"
#include <numeric>

namespace tbem {

std::vector<int> range(int min, int max) {
    std::vector<int> indices(max - min);
    std::iota(indices.begin(), indices.end(), min);
    return indices;
}

std::vector<int> range(int max) {
    return range(0, max);
}

std::vector<double> linspace(const double& a, const double& b, size_t count)
{
    std::vector<double> out(count);
    for (size_t i = 0; i < count; i++) {
        out[i] = a + i * (b - a) / ((double)count - 1);
    }
    return out;
}

} //END namespace tbem
