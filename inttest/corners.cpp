#include "3bem.h"
#include "laplace_kernels.h"
#include "laplace_solns.h"

using namespace tbem;

Vec2<double> center = {5, 0};
double r = 3.0;

template <typename F>
ConstraintEQ bc_constraint(const Mesh<2>& mesh, const F& fnc, 
        const VertexIterator<2>& dof) 
{
    return ConstraintEQ{
        {{(size_t)dof.absolute_index(), 1.0}},
        fnc(*dof)
    };
}

int main() {
    // Constraints on the end points!

    int refine_level = 7;
    double linear_solve_tol = 1e-9;
    QuadStrategy<2> qs(4, 4, 10, 3.0, 1e-5);
    auto surface = circle_mesh(center, r, refine_level);
    auto n_dofs = surface.n_dofs(); 
    size_t swap_bc = 1000;
    if (swap_bc > surface.n_facets()) {
        swap_bc = surface.n_facets();
    }

    auto continuity = mesh_continuity(surface.begin());
    std::vector<ConstraintEQ> potential_constraints = convert_to_constraints(continuity);
    for (auto it = surface.begin(); it != surface.begin() + 2 * swap_bc; ++it) {
        potential_constraints.push_back(bc_constraint(surface, theta_u, it));
    }
    auto potential_cm = from_constraints(potential_constraints);

    std::vector<ConstraintEQ> flux_constraints;
    for (auto it = surface.begin() + 2 * swap_bc;
              it != surface.begin() + surface.n_dofs(); ++it) {
        flux_constraints.push_back(bc_constraint(surface, ThetaDudn{center}, it));
    }
    auto flux_cm = from_constraints(flux_constraints);

    LaplaceSingle<2> single_kernel;
    LaplaceDouble<2> double_kernel;
    IdentityScalar<2> id_scalar;

    auto p_single = make_boundary_integral<2>(surface, surface, single_kernel);
    auto single_op = mesh_to_mesh_operator(p_single, qs);

    auto p_double = make_boundary_integral<2>(surface, surface, double_kernel);
    auto double_op = mesh_to_mesh_operator(p_double, qs);
    
    auto p_mass = make_boundary_integral<2>(surface, surface, id_scalar);
    auto mass_op = mass_operator(p_mass, qs);

    std::vector<double> rhs_potential(n_dofs, 0.0);
    std::vector<double> rhs_flux(n_dofs, 0.0);
    std::vector<VectorX> condensed{
        condense_vector(potential_cm, rhs_potential),
        condense_vector(flux_cm, rhs_flux),
    };
    auto dof_map = block_dof_map_from_functions(condensed);
    auto rhs = concatenate(dof_map, condensed);
    // If the RHS is all zero, PETSc just returns. This forces it to
    // actually solve the problem.
    rhs[0] += 2 * linear_solve_tol;

    int count = 0;
    auto soln_reduced = solve_system(rhs, linear_solve_tol,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            auto both = expand(dof_map, x);
            auto potential = distribute_vector(potential_cm, both[0], n_dofs);
            auto flux = distribute_vector(flux_cm, both[1], n_dofs);
            auto single_eval = single_op.apply({flux})[0];
            auto double_eval = double_op.apply({potential})[0];
            auto mass_eval = mass_op.apply({potential})[0];
            auto out = mass_eval;
            for (size_t i = 0; i < out.size(); i++) {
                out[i] += double_eval[i];
                out[i] -= single_eval[i];
            }
            auto y_temp_potential = condense_vector(potential_cm, out);
            auto y_temp_flux = condense_vector(flux_cm, out);
            for (size_t i = 0; i < y_temp_potential.size(); i++) {
                y[i] = y_temp_potential[i];
            }
            for (size_t i = 0; i < y_temp_flux.size(); i++) {
                y[y_temp_potential.size() + i] = y_temp_flux[i];
            }
        });

    auto both_soln = expand(dof_map, soln_reduced);
    auto potential_soln = distribute_vector(potential_cm, both_soln[0], n_dofs);
    auto flux_soln = distribute_vector(flux_cm, both_soln[1], n_dofs);
    auto single_eval = single_op.apply({flux_soln})[0];
    auto double_eval = double_op.apply({potential_soln})[0];
    auto mass_eval = mass_op.apply({potential_soln})[0];
    for (size_t i = 0;i < single_eval.size(); i++) {
        mass_eval[i] += double_eval[i];
        mass_eval[i] -= single_eval[i];
    }

    out_surface<2>(HDFOutputter("test_out/corners_potential.hdf5"), surface, {potential_soln});
    out_surface<2>(HDFOutputter("test_out/corners_flux.hdf5"), surface, {flux_soln});
    out_surface<2>(HDFOutputter("test_out/corners_single.hdf5"), surface, {single_eval});
    out_surface<2>(HDFOutputter("test_out/corners_double.hdf5"), surface, {double_eval});
    out_surface<2>(HDFOutputter("test_out/corners_all.hdf5"), surface, {mass_eval});
}
