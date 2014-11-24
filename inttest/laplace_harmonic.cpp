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
    double far_threshold = 3.0;
    int refine_level = 3;
    int near_quad_pts = 3;
    int near_steps = 7;
    int src_quad_pts = 3;
    //TODO: Something is seriously wrong when I use obs_quad_pts = 3
    int obs_quad_pts = 2;
    double tol = 5e-4;

    QuadStrategy qs(obs_quad_pts, src_quad_pts, near_quad_pts,
                    near_steps, far_threshold, tol);

    auto sphere = sphere_mesh(center, r).refine_repeatedly(refine_level);

    auto constraints = ConstraintMatrix::from_constraints(mesh_continuity(sphere));

    int n_dofs = 3 * sphere.facets.size();
    std::vector<double> u(n_dofs);
    for (unsigned int i = 0; i < sphere.facets.size(); i++) {
        for (int d = 0; d < 3; d++) {
            u[3 * i + d] = harmonic_u(sphere.facets[i].vertices[d]);
        }
    }

    std::vector<double> dudn(n_dofs);
    for (unsigned int i = 0; i < sphere.facets.size(); i++) {
        for (int d = 0; d < 3; d++) {
            dudn[3 * i + d] = harmonic_dudn(sphere.facets[i].vertices[d], center);
        }
    }

    TIC
    Problem p_double = {sphere, sphere, laplace_double, u};
    auto rhs_full = direct_interact(p_double, qs);
    TOC("RHS Eval")

    Problem p_mass = {sphere, sphere, one, u};
    auto rhs_mass = mass_term(p_mass, qs);

    for (unsigned int i = 0; i < rhs_full.size(); i++){
        rhs_full[i] = rhs_full[i] - rhs_mass[i];
    }

    auto rhs = constraints.get_reduced(rhs_full);
    for(auto r: rhs) {
        std::cout << r << std::endl;
    }
    std::cout << "N_dofs: " << rhs.size() << std::endl;

    TIC2
    Problem p_single = {sphere, sphere, laplace_single, {}};
    auto matrix = interact_matrix(p_single, qs);
    TOC("Matrix construct on " + std::to_string(sphere.facets.size()) + " facets");
    int count = 0;
    auto dudn_solved_subset = solve_system(rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            // std::cout << "iteration " << count << std::endl;
            count++;
            // TIC
            auto x_full = constraints.get_all(x, n_dofs);
            auto y_mult = bem_mat_mult(matrix, n_dofs, x_full); 
            auto y_temp = constraints.get_reduced(y_mult);
            /* TOC("Matrix multiply on " + std::to_string(sphere.facets.size()) + " faces"); */
            std::copy(y_temp.begin(), y_temp.end(), y.begin());
        });
    auto dudn_solved = constraints.get_all(dudn_solved_subset, n_dofs);
    std::cout << error_inf(dudn_solved, dudn) << std::endl;
    hdf_out("laplace.hdf5", sphere, dudn_solved); 
    return 0;

    double obs_len_scale = get_len_scale(sphere, 0, obs_quad_pts);
    for(int i = 0; i < 100; i++) {
        auto obs_pt = random_pt_sphere(center, obs_radius);

        auto obs_normal = normalized(center - obs_pt);
        ObsPt obs = {obs_len_scale, obs_pt, obs_normal};
       
        double double_layer = eval_integral_equation(p_double, qs, obs);
        Problem p_s = {sphere, sphere, laplace_single, dudn};
        double single_layer = eval_integral_equation(p_s, qs, obs);
        double result = double_layer - single_layer;
        double exact = 1.0 / hypot(obs_pt);
        double error = std::fabs(exact - result);
        if (error > 1e-2) {
            std::cout << "Failed with point: " << obs_pt << std::endl;
            std::cout << result << " " << exact << std::endl;
        }
    }
}
