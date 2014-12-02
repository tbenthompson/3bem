#include "vec.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "kernels.h"
#include "quadrature.h"
#include "util.h"
#include "petsc_interface.h"
#include "basis.h"
#include "bem.h"

/* This code:
 * 1. solves for the displacement at y = -0.5 in a full space due to the
 *    slip of a strike-slip fault on the line segment (0,-1)-(0,0).
 * 2. solves for the displacement at y = 0.0 in a half space due to the
 *    slip of a strike-slip fault on the line segment (0,-1)-(0,0).
 */

// Fault mesh.
auto fault = line_mesh({0, -1}, {0, 0}).refine_repeatedly(0);

// Unit slip on the fault plane.
std::vector<double> one_vec(2 * fault.facets.size(), 1.0);

void full_space() {
    // FULL SPACE STRIKE SLIP FAULT IN ANTIPLANE STRAIN

    // Surface mesh at y = -0.5
    auto surface1 = line_mesh({-10, -0.5}, {10, -0.5}).refine_repeatedly(9);

    // Quadrature details -- these parameters basically achieve machine precision
    double far_threshold = 4.0;
    int near_quad_pts = 17;
    int near_steps = 8;
    int src_quad_pts = 6;
    int obs_quad_pts = 2;
    double tol = 1e-13;
    QuadStrategy<2> qs(obs_quad_pts, src_quad_pts, near_quad_pts,
                         near_steps, far_threshold, tol);


    // Interpolate the integral equation onto the y = -0.5 surface
    Problem<2> p_fullspace = {fault, surface1, laplace_double<2>, one_vec};
    auto disp = interpolate(surface1, [&] (Vec2<double> x) {
            ObsPt<2> obs = {0.001, x, {0, 1}};
            double val = eval_integral_equation(p_fullspace, qs, obs);
            double exact = std::atan(0.5 / x[0]) / M_PI;
            if (std::fabs(exact - val) > 1e-13 && x[0] != 0) {
                std::cout << "FAILED Antiplane for x = " << x << std::endl;
            }
            return val;
        });

    hdf_out("antiplane_full_space.hdf5", surface1, disp);
}
    
void half_space() {
    // HALF SPACE STRIKE SLIP FAULT IN ANTIPLANE STRAIN

    QuadStrategy<2> qs(2);

    // Earth's surface
    auto surface2 = line_mesh({-25, 0.0}, {25, 0.0}).refine_repeatedly(10);
    auto raw_constraints = ConstraintMatrix::from_constraints(mesh_continuity(surface2));
    //TODO: generalize the breaking of constraints across a fault-surface intersection
    auto new_c_map = raw_constraints.c_map;
    for (std::size_t i = 0; i < surface2.facets.size(); i++) {
        auto x = surface2.facets[i].vertices[0][0];
        if (std::fabs(x) < 0.1) {
            new_c_map.erase(2 * i);
            new_c_map.erase(2 * i - 1);
        }
    }
    auto constraints = ConstraintMatrix{new_c_map};
    
    TIC
    // The RHS is the effect of the fault on the surface.
    Problem<2> p_rhs_halfspace = {fault, surface2, laplace_hypersingular<2>, one_vec};
    auto rhs_all_dofs = direct_interact(p_rhs_halfspace, qs);
    for (std::size_t i = 0; i < rhs_all_dofs.size(); i++) {
        rhs_all_dofs[i] = -rhs_all_dofs[i];
    }
    auto rhs = constraints.get_reduced(rhs_all_dofs);

    // The LHS is the effect of the surface on the surface.
    Problem<2> p_lhs_halfspace = {surface2, surface2, laplace_hypersingular<2>, one_vec};
    auto lhs = interact_matrix(p_lhs_halfspace, qs);

    // Solve the linear system.
    double linear_solve_tol = 1e-5;
    int n_dofs = 2 * surface2.facets.size();
    auto soln_reduced = solve_system(rhs, linear_solve_tol,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            auto x_full = constraints.get_all(x, n_dofs);
            auto y_mult = bem_mat_mult(lhs, n_dofs, x_full); 
            auto y_temp = constraints.get_reduced(y_mult);
            std::copy(y_temp.begin(), y_temp.end(), y.begin());
        });
    TOC("Solve antiplane half space.");
    auto soln = constraints.get_all(soln_reduced, n_dofs);

    hdf_out("antiplane_half_space.hdf5", surface2, soln);
}

int main() {
    full_space();
    half_space();
}
