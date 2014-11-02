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
    const Vec3<double> center = {5, 0, 0};
    double r = 3.0;
    double obs_radius = 2.7;
    double far_threshold = 2.0;
    int refine_level = 2;
    int near_field = 4;
    int src_quad_pts = 2;
    int obs_quad_pts = 3;

    auto sphere = sphere_mesh(center, r);
    for (int i = 0; i < refine_level; i++) {
        sphere = refine_mesh(sphere);
    }
    sphere = clean_mesh(sphere);
    auto q_src = tri_gauss(src_quad_pts);
    auto q_obs = tri_gauss(obs_quad_pts);
    KernelFnc K = laplace_single;
    KernelFnc Kdn = laplace_double;
    NearEval ne(near_field);

    int n_verts = sphere.vertices.size();
    std::vector<double> u(n_verts);
    for (int i = 0; i < n_verts; i++) {
        u[i] = harmonic_u(sphere.vertices[i]);
    }

    std::vector<double> dudn(n_verts);
    for (int i = 0; i < n_verts; i++) {
        dudn[i] = harmonic_dudn(sphere.vertices[i], center);
    }

    auto rhs = direct_interact(sphere, sphere, q_src, q_obs, 
                               Kdn, u, near_field, far_threshold);
    auto rhs_mass = mass_term(sphere, q_src, u);

    for (unsigned int i = 0; i < rhs.size(); i++){
        //TODO: I think the signs here are wrong. Where is there a sign flip?
        rhs[i] = -rhs[i] + rhs_mass[i];
    }

    int count = 0;
    auto dudn_solved = solve_system(rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            TIC
            auto y_temp = direct_interact(sphere, sphere, q_src, q_obs,
                                          K, x, near_field, far_threshold);
            TOC("Direct interact on " + std::to_string(sphere.faces.size()) + " faces");
            std::copy(y_temp.begin(), y_temp.end(), y.begin());
        });
    std::cout << error_inf(dudn_solved, dudn) << std::endl;
    hdf_out("laplace.hdf5", sphere, dudn_solved); 

    double obs_len_scale = get_len_scale(sphere, 0, obs_quad_pts);
    for(int i = 0; i < 100; i++) {
        auto obs_pt = random_pt_sphere(center, obs_radius);

        auto obs_normal = normalized(center - obs_pt);

        double u_effect = eval_integral_equation(sphere, q_src, Kdn, ne, obs_pt,
                                               obs_normal, obs_len_scale, u);
        double dudn_effect = eval_integral_equation(sphere, q_src, K, ne, obs_pt,
                                         obs_normal, obs_len_scale, dudn);
        double result = u_effect + dudn_effect;
        double exact = 1.0 / hypot(obs_pt);
        double error = std::fabs(exact - result);
        if (error > 5e-2) {
            std::cout << "Failed with point: " << obs_pt << std::endl;
        }
    }
}
