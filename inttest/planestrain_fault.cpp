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
    std::vector<Vec2<double>> du(n_fault_dofs, constant<Vec2<double>>::make(slip));

    std::vector<Vec2<double>> all_dofs_rhs(n_surface_dofs, zeros<Vec2<double>>::make());

    TIC
    auto p_rhs = make_problem<2>(fault, surface, hyp, du);
    auto res = direct_interact(p_rhs, qs);
    for (unsigned int i = 0; i < res.size(); i++) {
        all_dofs_rhs[i] += res[i];
    }
    auto rhs = constraints.get_reduced(all_dofs_rhs);
    TOC("Building RHS");

    int n_reduced_surface_dofs = 2 * rhs.size();

    TIC2
    auto p_lhs = make_problem<2>(surface, surface, hyp, {});
    auto lhs = interact_matrix(p_lhs, qs);
    TOC("Building LHS matrices");

    int count = 0;
    auto disp_reduced = solve_system((double*)rhs.data(),
                                             n_reduced_surface_dofs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            auto x_vec_reduced = reinterpret_vector<Vec2<double>>(x);
            auto x_vec = constraints.get_all(x_vec_reduced, n_surface_dofs);
            auto y_vec = bem_mat_mult(lhs, hyp, n_surface_dofs, x_vec);
            auto y_vec_reduced = constraints.get_reduced(y_vec);
            for (std::size_t i = 0; i < y_vec_reduced.size(); i++) {
                y[2 * i] = y_vec_reduced[i][0];
                y[2 * i + 1] = y_vec_reduced[i][1];
            }
        }
    );

    auto disp_reduced_vec = reinterpret_vector<Vec2<double>>(disp_reduced);
    auto disp_vec = constraints.get_all(disp_reduced_vec, n_surface_dofs);

    hdf_out_surface<2,2>("planestrain_thrust.hdf5", surface, disp_vec);
}
