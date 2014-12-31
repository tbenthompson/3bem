#include "3bem.h"

using namespace tbem;

int main() {
    // HALF SPACE THRUST FAULT IN PLANE STRAIN

    // Fault mesh.
    auto fault = line_mesh({-1, -1}, {0, 0}).refine_repeatedly(0);

    // Simple low accuracy quadrature strategy.
    QuadStrategy<2> qs(2);

    // Earth's surface
    auto surface = line_mesh({-25, 0.0}, {25, 0.0}).refine_repeatedly(8);
    auto raw_constraints =
        ConstraintMatrix::from_constraints(mesh_continuity<2>(surface));
    auto constraints = apply_discontinuities<2>(surface, fault, raw_constraints);

    double shear_mod = 30e9;
    double poisson = 0.25;
    ElasticHypersingular<2> hyp(shear_mod, poisson);

    std::size_t n_fault_dofs = 2 * fault.facets.size();
    std::size_t n_surface_dofs = 2 * surface.facets.size();

    double slip = 1;
    std::vector<Vec2<double>> du(n_fault_dofs, slip * ones<Vec2<double>>::make());

    std::vector<Vec2<double>> all_dofs_rhs(n_surface_dofs, zeros<Vec2<double>>::make());

    auto p = make_problem<2>(fault, surface, hyp, du);
    auto res = direct_interact(p, qs);
    for (unsigned int i = 0; i < res.size(); i++) {
        all_dofs_rhs[i] += res[i];
    }

    // std::array<std::vector<double>,2> rhs = {
    //     constraints.get_reduced(all_dofs_rhs[0]),
    //     constraints.get_reduced(all_dofs_rhs[1]),
    // };

    // int n_reduced_surface_dofs = rhs[0].size();

    // std::vector<double> vector_rhs(2 * n_reduced_surface_dofs);
    // for (int d = 0; d < 2; d++) {
    //     std::copy(rhs[d].begin(), rhs[d].end(), vector_rhs.begin() + 
    //                                             d * n_reduced_surface_dofs);
    // }

    // std::array<std::array<std::vector<double>,2>,2> mats;
    // TIC
    // for (int k = 0; k < 2; k++) {
    //     for (int j = 0; j < 2; j++) {
    //         Problem<2> p = {surface, surface, ek.hypersingular_mat[k][j], {}};
    //         mats[k][j] = interact_matrix(p, qs);
    //     }
    // }
    // TOC("Building matrices")

    // int count = 0;
    // auto surface_disp = solve_system(vector_rhs, 1e-5,
    //     [&] (std::vector<double>& x, std::vector<double>& y) {
    //         std::cout << "iteration " << count << std::endl;
    //         count++;
    //         std::array<std::vector<double>,2> x_temp;
    //         for (int i = 0; i < 2; i++) {
    //             auto x_temp_reduced = std::vector<double>(n_surface_dofs);
    //             std::copy(x.begin() + i * n_reduced_surface_dofs,
    //                       x.begin() + (i + 1) * n_reduced_surface_dofs,
    //                       x_temp_reduced.begin());
    //             x_temp[i] = constraints.get_all(x_temp_reduced, n_surface_dofs);
    //         }

    //         std::array<std::vector<double>,2> y_temp;
    //         for (int k = 0; k < 2; k++) {
    //             auto y_temp_full = std::vector<double>(n_surface_dofs, 0.0);
    //             for (int j = 0; j < 2; j++) {
    //                 for (unsigned int mi = 0; mi < n_surface_dofs; mi++) {
    //                     for (unsigned int ni = 0; ni < n_surface_dofs; ni++) {
    //                         y_temp_full[mi] += 
    //                             mats[k][j][mi * n_surface_dofs + ni] * x_temp[j][ni];
    //                     }
    //                 }
    //             }
    //             auto y_temp_reduced = constraints.get_reduced(y_temp_full);
    //             std::copy(y_temp_reduced.begin(), y_temp_reduced.end(),
    //                       y.begin() + k * n_reduced_surface_dofs);
    //         }
    //     });

    // std::array<std::vector<double>,2> soln;
    // for (int i = 0; i < 2; i++) {
    //     auto reduced_soln = std::vector<double>(n_reduced_surface_dofs);
    //     std::copy(surface_disp.begin() + i * n_reduced_surface_dofs,
    //               surface_disp.begin() + (i + 1) * n_reduced_surface_dofs,
    //               reduced_soln.begin());
    //     soln[i] = constraints.get_all(reduced_soln, n_surface_dofs);
    // }
    // hdf_out_surface<2>("2dthrust0.hdf5", surface, {soln[0], soln[1]});
}
