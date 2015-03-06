#ifndef __QWELKASHNCN_LAPLACE_SOLNS_H
#define __QWELKASHNCN_LAPLACE_SOLNS_H

#include "laplace.h"

double log_u(Vec2<double> x) {
    return std::log(std::sqrt(x[0] * x[0] + x[1] * x[1]));
}

double theta_u(Vec2<double> x) {
    return std::atan2(x[1], x[0]);
}

struct LogDudn {
    const Vec2<double> center;
    double operator()(Vec2<double> loc) const {
        auto n = normalized(center - loc);
        return dot_product(n, loc) / hypot2(loc);
    }
};
    
struct ThetaDudn {
    const Vec2<double> center;
    double operator()(Vec2<double> loc) const {
        auto n = normalized(center - loc);
        double x = loc[0];
        double y = loc[1];
        double dy = 1.0 / (x * (1 + (y * y / (x * x))));
        double dx = (-y / x) * dy;
        return dot_product(n, Vec2<double>{dx, dy});
    }
};
#endif
