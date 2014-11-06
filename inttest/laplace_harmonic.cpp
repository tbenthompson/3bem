#include "kernels.h"
#include "bem.h"
#include "quadrature.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "util.h"
#include "petsc_interface.h"
#include <iostream>

double harmonic_u(Vec3<double> x) {
    return 1.0 / hypot(x);
}

double harmonic_dudn(Vec3<double> x, Vec3<double> center) {
    return -(sum(x * normalized(x - center))) / pow(hypot(x), 3);
}

int main() {
    //THIS IS HOT!
    //TODO: Solving for dudn here should not have equality constraints
    //imposed
    const Vec3<double> center = {5, 0, 0};
    double r = 3.0;
    double obs_radius = 2.7;
    double far_threshold = 2.0;
    int refine_level = 3;
    int near_quad_pts = 5;
    int adjacent_quad_pts = 5;
    int near_steps = 5;
    int src_quad_pts = 2;
    int obs_quad_pts = 2;

    auto sphere = refine_clean(sphere_mesh(center, r), refine_level);

    auto q_src = tri_gauss(src_quad_pts);
    auto q_obs = tri_gauss(obs_quad_pts);
    auto K = laplace_single<double>;
    auto Kdn = laplace_double<double>;
    NearEval ne(near_quad_pts, adjacent_quad_pts, near_steps, q_obs);

    int n_verts = sphere.vertices.size();
    std::vector<double> u(n_verts);
    for (int i = 0; i < n_verts; i++) {
        u[i] = harmonic_u(sphere.vertices[i]);
    }

    std::vector<double> dudn(n_verts);
    for (int i = 0; i < n_verts; i++) {
        dudn[i] = harmonic_dudn(sphere.vertices[i], center);
    }

    TIC
    auto rhs = direct_interact(sphere, sphere, q_src, q_obs, 
                               Kdn, u, ne, far_threshold);
    TOC("RHS Eval")

    auto rhs_mass = mass_term(sphere, q_src, u);

    for (unsigned int i = 0; i < rhs.size(); i++){
        //TODO: I think the signs here are wrong. Where is there a sign flip?
        rhs[i] = -rhs[i] + rhs_mass[i];
    }

    TIC2
    auto matrix = interact_matrix(sphere, sphere, q_src, q_obs,
                                  K, ne, far_threshold);
    TOC("Matrix construct on " + std::to_string(sphere.faces.size()) + " faces");
    int count = 0;
    auto dudn_solved = solve_system(rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            TIC
            for (int i = 0; i < n_verts; i++) {
                y[i] = 0.0;
                for (int j = 0; j < n_verts; j++) {
                    y[i] += matrix[i][j] * x[j];
                }
            }
            TOC("Matrix multiply on " + std::to_string(sphere.faces.size()) + " faces");
        });
    std::cout << error_inf(dudn_solved, dudn) << std::endl;
    hdf_out("laplace.hdf5", sphere, dudn_solved); 

    double obs_len_scale = get_len_scale(sphere, 0, obs_quad_pts);
    for(int i = 0; i < 100; i++) {
        auto obs_pt = random_pt_sphere(center, obs_radius);

        auto obs_normal = normalized(center - obs_pt);

        double u_effect = eval_integral_equation(sphere, q_src, Kdn, ne, -1, obs_pt,
                                               obs_normal, obs_len_scale, u);
        double dudn_effect = eval_integral_equation(sphere, q_src, K, ne, -1, obs_pt,
                                         obs_normal, obs_len_scale, dudn);
        double result = u_effect + dudn_effect;
        double exact = 1.0 / hypot(obs_pt);
        double error = std::fabs(exact - result);
        if (error > 5e-2) {
            std::cout << "Failed with point: " << obs_pt << std::endl;
        }
    }
}
