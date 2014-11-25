#include "kernels.h"
#include "bem.h"
#include "quadrature.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "util.h"
#include "petsc_interface.h"
#include <iostream>

double harmonic_u(Vec2<double> x) {
    return std::atan(x[1] / x[0]);
}

double harmonic_dudn(Vec2<double> loc, Vec2<double> center) {
    auto dir = loc - center;
    double x = loc[0];
    double y = loc[1];
    double dy = 1.0 / (x * (1 + (y * y / (x * x))));
    double dx = (-y / x) * dy;
    return dir[0] * dx + dir[1] * dy;
}

int main() {
    //THIS IS HOT!
    //TODO: Very similar to 3D case. Can they be templated and combined, maybe
    //with the harmonic function part abstracted out?
    const Vec2<double> center = {5, 5};
    double r = 3.0;
    double obs_radius = 2.7;
    double far_threshold = 3.0;
    int refine_level = 5;
    int near_quad_pts = 3;
    int near_steps = 5;
    int src_quad_pts = 2;
    //TODO: Something is seriously wrong when I use obs_quad_pts = 3
    int obs_quad_pts = 2;
    double tol = 1e-4;

    QuadStrategy<2> qs(obs_quad_pts, src_quad_pts, near_quad_pts,
                    near_steps, far_threshold, tol);

    auto circle = circle_mesh(center, r).refine_repeatedly(refine_level);
    std::cout << circle.facets.size() << std::endl;

    auto constraints = ConstraintMatrix::from_constraints(mesh_continuity(circle));

    int n_dofs = 2 * circle.facets.size();
    std::vector<double> u(n_dofs);
    for (unsigned int i = 0; i < circle.facets.size(); i++) {
        for (int d = 0; d < 2; d++) {
            u[2 * i + d] = harmonic_u(circle.facets[i].vertices[d]);
        }
    }

    std::vector<double> dudn(n_dofs);
    for (unsigned int i = 0; i < circle.facets.size(); i++) {
        for (int d = 0; d < 2; d++) {
            dudn[2 * i + d] = harmonic_dudn(circle.facets[i].vertices[d], center);
        }
    }

    TIC
    Problem<2> p_double = {circle, circle, laplace_double2d, u};
    auto rhs_full = direct_interact(p_double, qs);
    TOC("RHS Eval")

    Problem<2> p_mass = {circle, circle, one<2>, u};
    auto rhs_mass = mass_term(p_mass, qs);

    for (unsigned int i = 0; i < rhs_full.size(); i++){
        rhs_full[i] = -rhs_full[i] + rhs_mass[i];
    }

    auto rhs = constraints.get_reduced(rhs_full);
    std::cout << "N_dofs: " << rhs.size() << std::endl;

    TIC2
    Problem<2> p_single = {circle, circle, laplace_single2d, {}};
    auto matrix = interact_matrix(p_single, qs);
    TOC("Matrix construct on " + std::to_string(circle.facets.size()) + " facets");
    int count = 0;
    // return 0;
    auto dudn_solved_subset = solve_system(rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            TIC
            auto x_full = constraints.get_all(x, n_dofs);
            auto y_mult = bem_mat_mult(matrix, n_dofs, x_full); 
            auto y_temp = constraints.get_reduced(y_mult);
            TOC("Matrix multiply on " + std::to_string(circle.facets.size()) + " faces");
            std::copy(y_temp.begin(), y_temp.end(), y.begin());
        });
    auto dudn_solved = constraints.get_all(dudn_solved_subset, n_dofs);
    std::cout << error_inf(dudn_solved, dudn) << std::endl;
    // hdf_out("laplace2d.hdf5", circle, dudn_solved); 
    return 0;

    // double obs_len_scale = get_len_scale(circle, 0, obs_quad_pts);
    // for(int i = 0; i < 100; i++) {
    //     auto obs_pt = random_pt_circle(center, obs_radius);

    //     auto obs_normal = normalized(center - obs_pt);
    //     ObsPt<2> obs = {obs_len_scale, obs_pt, obs_normal};
    //    
    //     double double_layer = eval_integral_equation(p_double, qs, obs);
    //     Problem<2> p_s = {circle, circle, laplace_single, dudn};
    //     double single_layer = eval_integral_equation(p_s, qs, obs);
    //     double result = double_layer - single_layer;
    //     double exact = 1.0 / hypot(obs_pt);
    //     double error = std::fabs(exact - result);
    //     if (error > 1e-2) {
    //         std::cout << "Failed with point: " << obs_pt << std::endl;
    //         std::cout << result << " " << exact << std::endl;
    //     }
    // }
}
