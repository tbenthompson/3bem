#include "laplace.h"

Vec3<double> center = {5, 0, 0};
double r = 3.0;

double harmonic_u(Vec3<double> loc) {
    return 1.0 / hypot(loc);
}

struct InvDudn {
    double operator()(Vec3<double> loc) const {
        auto n = normalized(center - loc);
        return dot_product(loc, n) / pow(hypot(loc), 3);
    }
};

int main() {
    double obs_radius = 2.9;
    int refine_level = 3;
    int n_test_pts = 100;
    auto sphere = sphere_mesh(center, r, refine_level);
    std::vector<Vec3<double>> test_pts(n_test_pts);
    for (int i = 0; i < n_test_pts; i++) {
        test_pts[i] = random_pt_sphere<3>(center, random_val() * obs_radius);
    }
    dirichlet_laplace_test<3>(sphere, test_pts, harmonic_u, InvDudn());
}
