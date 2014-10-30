#ifndef __YUIQWEOIUQWE_KERNELS_H
#define __YUIQWEOIUQWE_KERNELS_H

#include <array>
#include <functional>

#include "vec.h"

inline double one(double r2, 
                  Vec3<double> delta,
                  Vec3<double> nsrc,
                  Vec3<double> nobs) {
    return 1.0;
}

inline double laplace_single(double r2,
                             Vec3<double> delta,
                             Vec3<double> nsrc,
                             Vec3<double> nobs) {
    return 1.0 / (4.0 * M_PI * std::sqrt(r2));
}

inline double laplace_double(double r2,
                             Vec3<double> delta,
                             Vec3<double> nsrc,
                             Vec3<double> nobs) {
    return -(nsrc[0] * delta[0] + nsrc[1] * delta[1] + nsrc[2] * delta[2]) / 
           (4 * M_PI * pow(r2, 1.5));
}

#endif
