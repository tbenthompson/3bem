#include "tanh_sinh.h"
#include <cmath>
//TODO: REMOVE
#include <iostream>

namespace tbem {

QuadRule<1> tanh_sinh(size_t n)
{
    const double log_offset = 1.50;
    double h = std::exp(-std::log(static_cast<double>(n)) + log_offset);
    std::cout << h << std::endl;
    int ni = static_cast<int>(n);
    QuadRule<1> retval;
    for (int i = -ni; i <= ni; i++) {
        auto sinhterm = 0.5 * M_PI * std::sinh(i * h);
        auto coshsinh = std::cosh(sinhterm);
        auto x = std::tanh(sinhterm);
        auto w = h * 0.5 * M_PI * std::cosh(i * h) / (coshsinh * coshsinh);
        retval.push_back(QuadPt<1>{x, w});
    }
    return retval;
}

} //end namespace tbem
