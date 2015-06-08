#include "numbers.h"

namespace tbem {

std::vector<double> linspace(const double& a, const double& b, size_t count)
{
    std::vector<double> out(count);
    for (size_t i = 0; i < count; i++) {
        out[i] = a + i * (b - a) / ((double)count - 1);
    }
    return out;
}

} //END namespace tbem
