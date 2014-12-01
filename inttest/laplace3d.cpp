#include "laplace.h"
#include "mesh_gen.h"
#include <iostream>

double harmonic_u(Vec3<double> x) {
    return 1.0 / hypot(x);
}

double harmonic_dudn(Vec3<double> x, Vec3<double> n) {
    return dot(x, n) / pow(hypot(x), 3);
}

int main() {
    //TODO: Something is seriously wrong when I use obs_quad_pts = 3
    const Vec3<double> center = {5, 0, 0};
    double r = 3.0;
    double obs_radius = 2.7;
    int refine_level = 3;
    int n_test_pts = 100;
    auto sphere = sphere_mesh(center, r).refine_repeatedly(refine_level);
    std::vector<Vec3<double>> test_pts(n_test_pts);
    for (int i = 0; i < n_test_pts; i++) {
        test_pts[i] = random_pt_sphere<3>(center, random_val() * obs_radius);
    }
    run_laplace_test<3>(sphere, test_pts, harmonic_u, harmonic_dudn);
}
