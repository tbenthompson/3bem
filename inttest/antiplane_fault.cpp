#include "3bem.h"
#include "laplace_kernels.h"
using namespace tbem;

/* This code:
 * 1. solves for the displacement at y = -0.5 in a full space due to the
 *    slip of a strike-slip fault on the line segment (0,-1)-(0,0).
 * 2. solves for the displacement at y = 0.0 in a half space due to the
 *    slip of a strike-slip fault on the line segment (0,-1)-(0,0).
 */

// Fault mesh.
auto fault = line_mesh({0, -1}, {0, 0}).refine_repeatedly(0);

// Unit slip on the fault plane.
std::vector<double> one_vec(fault.n_dofs(), 1.0);


void full_space() {
    // FULL SPACE STRIKE SLIP FAULT IN ANTIPLANE STRAIN

    // Surface mesh at y = -0.5
    auto surface1 = line_mesh({-10, -0.5}, {10, -0.5}).refine_repeatedly(14);

    auto continuity = mesh_continuity(surface1.begin());
    auto cut_cont = cut_at_intersection(continuity, surface1.begin(), fault.begin());
    auto constraints = convert_to_constraints(cut_cont);
    auto constraint_matrix = ConstraintMatrix::from_constraints(constraints);

    // Quadrature details -- these parameters basically achieve machine precision
    double far_threshold = 4.0;
    int near_steps = 8;
    int src_quad_pts = 6;
    int obs_quad_pts = 2;
    double tol = 1e-13;
    QuadStrategy<2> qs(obs_quad_pts, src_quad_pts, near_steps, far_threshold, tol);


    // Interpolate the integral equation onto the y = -0.5 surface
    TIC
    auto p_fullspace = make_problem<2>(fault, surface1, LaplaceDouble<2>(), one_vec);
    auto disp = constrained_interpolate<2>(surface1, [&] (Vec2<double> x) {
            ObsPt<2> obs = {0.001, x, {0, 1}, {0, -1}};
            double val = eval_integral_equation(p_fullspace, qs, obs);
            double exact = std::atan(0.5 / x[0]) / M_PI;
            if (std::fabs(exact - val) > 1e-13 && x[0] != 0) {
                std::cout << "FAILED Antiplane for x = " << x << std::endl;
            }
            return val;
        }, constraint_matrix);
    TOC("Solve fullspace antiplane strike slip motion")

    auto file = HDFOutputter("test_out/antiplane_full_space.hdf5");
    out_surface<2>(file, surface1, disp, 1);
}
    
void half_space() {
    // HALF SPACE STRIKE SLIP FAULT IN ANTIPLANE STRAIN

    // Simple low accuracy quadrature strategy.
    QuadStrategy<2> qs(2);

    // Earth's surface
    auto surface2 = line_mesh({-50, 0.0}, {50, 0.0}).refine_repeatedly(9);
    
    auto continuity = mesh_continuity(surface2.begin());
    auto cut_cont = cut_at_intersection(continuity, surface2.begin(), fault.begin());
    auto constraints = convert_to_constraints(cut_cont);
    auto constraint_matrix = ConstraintMatrix::from_constraints(constraints);
    
    TIC
    // The RHS is the effect of the fault on the surface.
    auto p_rhs_halfspace =
        make_problem<2>(fault, surface2, LaplaceHypersingular<2>(), one_vec);

    auto rhs_all_dofs = direct_interact(p_rhs_halfspace, qs);
    for (std::size_t i = 0; i < rhs_all_dofs.size(); i++) {
        rhs_all_dofs[i] = -rhs_all_dofs[i];
    }
    auto rhs = constraint_matrix.get_reduced(rhs_all_dofs);

    // The LHS is the effect of the surface on the surface.
    auto hypersingular_kernel = LaplaceHypersingular<2>();
    auto p_lhs_halfspace = make_problem<2>(surface2, surface2,
                                           hypersingular_kernel, one_vec);
    auto lhs = interact_matrix(p_lhs_halfspace, qs);

    // Solve the linear system.
    double linear_solve_tol = 1e-5;
    auto soln_reduced = solve_system(rhs, linear_solve_tol,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            auto x_full = constraint_matrix.get_all(x, surface2.n_dofs());
            auto y_mult = bem_mat_mult(lhs, hypersingular_kernel, surface2.n_dofs(), x_full); 
            auto y_temp = constraint_matrix.get_reduced(y_mult);
            std::copy(y_temp.begin(), y_temp.end(), y.begin());
        });
    TOC("Solve antiplane half space.");
    auto soln = constraint_matrix.get_all(soln_reduced, surface2.n_dofs());

    auto filesurface = HDFOutputter("test_out/antiplane_half_space.hdf5");
    out_surface<2>(filesurface, surface2, soln, 1);


    // Loop over a bunch of interior points and evaluate the displacement
    int nx = 200; int ny = 200;
    int n_tot = nx * ny;
    std::vector<Vec2<double>> interior_pts(n_tot);
    std::vector<double> interior_disp(n_tot);
    std::vector<double> interior_trac[2];
    interior_trac[0].resize(n_tot);
    interior_trac[1].resize(n_tot);

    Vec2<double> x_bounds = {-5, 5};
    Vec2<double> y_bounds = {-5, 0};
    TIC2
#pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        double x = x_bounds[0] + i * (x_bounds[1] - x_bounds[0]) / ((double)nx - 1);
        for (int j = 0; j < ny; j++) {
            double y = y_bounds[0] + j * (y_bounds[1] - y_bounds[0]) / ((double)ny - 1);
            Vec2<double> pt = {x, y};
            interior_pts[i * ny + j] = pt;

            ObsPt<2> obs[] = {
                {0.001, pt, {0, 1}, {0, -1}},
                {0.001, pt, {1, 0}, {0, -1}}
            };

            auto p_disp_fault = make_problem<2>(fault, surface2,
                                             LaplaceDouble<2>(), one_vec);
            auto p_disp_surf = make_problem<2>(surface2, surface2,
                                            LaplaceDouble<2>(), soln);
            double eval_fault = eval_integral_equation(p_disp_fault, qs, obs[0]);
            double eval_surf = eval_integral_equation(p_disp_surf, qs, obs[0]);
            double eval = eval_surf + eval_fault;
            interior_disp[i * ny + j] = eval;

            for (int d = 0; d < 2; d++) {
                auto p_trac_fault = make_problem<2>(fault, surface2,
                                                LaplaceHypersingular<2>(), one_vec);
                auto p_trac_surf = make_problem<2>(surface2, surface2,
                                                LaplaceHypersingular<2>(), soln);
                eval_fault = eval_integral_equation(p_trac_fault, qs, obs[d]);
                eval_surf = eval_integral_equation(p_trac_surf, qs, obs[d]);
                eval = eval_surf + eval_fault;
                interior_trac[d][i * ny + j] = eval;
            }

        }
    }
    auto fileu = HDFOutputter("test_out/antiplane_half_space_volu.hdf5");
    out_volume<2>(fileu, interior_pts, interior_disp, 1);
    auto filetx = HDFOutputter("test_out/antiplane_half_space_voltx.hdf5");
    out_volume<2>(filetx, interior_pts, interior_trac[0], 1);
    auto filety = HDFOutputter("test_out/antiplane_half_space_volty.hdf5");
    out_volume<2>(filety, interior_pts, interior_trac[1], 1);
    TOC("Interior eval.")
}

int main() {
    full_space();
    half_space();
}
