#include "mesh_gen.h"
#include "laplace.h"
#include <iostream>

double log_u(Vec2<double> x) {
    return std::log(std::sqrt(x[0] * x[0] + x[1] * x[1]));
}

double theta_u(Vec2<double> x) {
    return std::atan(x[1] / x[0]);
}

double log_dudn(Vec2<double> loc, Vec2<double> n) {
    return dot(n, loc) / hypot2(loc);
}
    
double theta_dudn(Vec2<double> loc, Vec2<double> n) {
    double x = loc[0];
    double y = loc[1];
    double dy = 1.0 / (x * (1 + (y * y / (x * x))));
    double dx = (-y / x) * dy;
    return dot(n, {dx, dy});
}

int main() {
    const Vec2<double> center = {5, 0};
    double r = 3.0;
    double obs_radius = 2.9;
    int refine_level = 10;
    int n_test_pts = 100;
    auto circle = circle_mesh(center, r).refine_repeatedly(refine_level);
    std::vector<Vec2<double>> test_pts(n_test_pts);
    for (int i = 0; i < n_test_pts; i++) {
        test_pts[i] = random_pt_sphere<2>(center, random_val() * obs_radius);
    }
    run_laplace_test<2>(circle, test_pts, log_u, log_dudn);
    run_laplace_test<2>(circle, test_pts, theta_u, theta_dudn);
}
