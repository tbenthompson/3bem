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
    int near_steps = 5;
    int src_quad_pts = 2;
    int obs_quad_pts = 2;

    QuadStrategy qs(obs_quad_pts, src_quad_pts, near_quad_pts,
                    near_steps, far_threshold);

    auto sphere = refine_clean(sphere_mesh(center, r), refine_level);

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
    Problem p_double = {sphere, sphere, laplace_double, u};
    auto rhs = direct_interact(p_double, qs);
    TOC("RHS Eval")

    Problem p_mass = {sphere, sphere, one, u};
    auto rhs_mass = mass_term(p_mass, qs);

    for (unsigned int i = 0; i < rhs.size(); i++){
        //TODO: I think the signs here are wrong. Where is there a sign flip?
        rhs[i] = -rhs[i] + rhs_mass[i];
    }

    TIC2
    Problem p_single = {sphere, sphere, laplace_single, {}};
    auto matrix = interact_matrix(p_single, qs);
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
    return 0;

    double obs_len_scale = get_len_scale(sphere, 0, obs_quad_pts);
    for(int i = 0; i < 100; i++) {
        auto obs_pt = random_pt_sphere(center, obs_radius);

        auto obs_normal = normalized(center - obs_pt);
        ObsPt obs = {obs_len_scale, obs_pt, obs_normal};
       
        double u_effect = eval_integral_equation(p_double, qs, obs);
        Problem p_s = {sphere, sphere, laplace_single, dudn};
        double dudn_effect = eval_integral_equation(p_s, qs, obs);
        double result = u_effect + dudn_effect;
        double exact = 1.0 / hypot(obs_pt);
        double error = std::fabs(exact - result);
        if (error > 5e-2) {
            std::cout << "Failed with point: " << obs_pt << std::endl;
        }
    }
}
