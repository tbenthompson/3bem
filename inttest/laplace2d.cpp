#include "laplace_solns.h"
#include "laplace.h"

Vec2<double> center = {5, 0};
double r = 3.0;

int main() {
    double obs_radius = 2.9;
    int refine_level = 4;
    int n_test_pts = 100;
    auto circle = circle_mesh(center, r, refine_level);
    std::vector<Vec2<double>> test_pts(n_test_pts);
    for (int i = 0; i < n_test_pts; i++) {
        test_pts[i] = random_pt_sphere(center, random_val() * obs_radius);
    }
    dirichlet_laplace_test<2>(circle, test_pts, log_u, LogDudn{center});
    dirichlet_laplace_test<2>(circle, test_pts, theta_u, ThetaDudn{center});
}
