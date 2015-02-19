#include "disloc_shared.h"

using namespace tbem;

int main() {
    // HALF SPACE THRUST FAULT IN PLANE STRAIN

    // Fault mesh.
    auto fault = line_mesh({-1, -1}, {0, 0}).refine_repeatedly(0);

    // Simple low accuracy quadrature strategy.
    QuadStrategy<2> qs(3);

    // Earth's surface
    auto surface = line_mesh({-100, 0.0}, {100, 0.0}).refine_repeatedly(9);

    auto constraint_matrix = surf_fault_constraints(surface.begin(), fault.begin());

    double shear_mod = 30e9;
    double poisson = 0.25;
    ElasticHypersingular<2> hyp(shear_mod, poisson);

    double slip = 1;
    Function duxy(fault.n_dofs(), slip);
    BlockFunction du{duxy, duxy};

    TIC
    auto p_rhs = make_problem<2>(surface, fault, hyp);
    auto rhs_op = mesh_to_mesh_operator(p_rhs, qs);
    auto res = apply_operator(rhs_op, du);
    auto all_dofs_rhs = apply_operator(rhs_op, du);
    BlockFunction condensed{
        condense_vector(constraint_matrix, all_dofs_rhs[0]),
        condense_vector(constraint_matrix, all_dofs_rhs[1])
    };
    auto dof_map = block_dof_map_from_functions(condensed);
    auto rhs = concatenate(dof_map, condensed);
    TOC("Building RHS");

    TIC2
    auto p_lhs = make_problem<2>(surface, surface, hyp);
    auto lhs = mesh_to_mesh_operator(p_lhs, qs);
    TOC("Building LHS matrices");

    auto disp_reduced = solve(lhs, rhs, dof_map, surface, constraint_matrix);

    auto disp_reduced_vec = expand(dof_map, disp_reduced);
    BlockFunction soln{
        distribute_vector(constraint_matrix, disp_reduced_vec[0], surface.n_dofs()),
        distribute_vector(constraint_matrix, disp_reduced_vec[1], surface.n_dofs())
    };

    auto file = HDFOutputter("test_out/planestrain_u.hdf5");
    out_surface(file, surface, soln);
}
