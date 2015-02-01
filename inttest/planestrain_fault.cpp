#include "disloc_shared.h"
#include "armadillo_facade.h"
#include <armadillo>

using namespace tbem;

int main() {
    // HALF SPACE THRUST FAULT IN PLANE STRAIN

    // Fault mesh.
    auto fault = line_mesh({-1, -1}, {0, 0}).refine_repeatedly(0);

    // Simple low accuracy quadrature strategy.
    QuadStrategy<2> qs(3);

    // Earth's surface
    auto surface = line_mesh({-100, 0.0}, {100, 0.0}).refine_repeatedly(10);

    auto constraint_matrix = surf_fault_constraints(surface.begin(), fault.begin());

    double shear_mod = 30e9;
    double poisson = 0.25;
    ElasticHypersingular<2> hyp(shear_mod, poisson);

    double slip = 1;
    std::vector<double> duxy(fault.n_dofs(), slip);
    std::vector<std::vector<double>> du{duxy, duxy};

    TIC
    auto p_rhs = make_problem<2>(fault, surface, hyp);
    auto rhs_op = mesh_to_mesh_operator(p_rhs, qs);
    auto res = apply_operator(rhs_op, du);
    auto all_dofs_rhs = apply_operator(rhs_op, du);
    auto rhs = concatenate({
        condense_vector(constraint_matrix, all_dofs_rhs[0]),
        condense_vector(constraint_matrix, all_dofs_rhs[1])
    });
    TOC("Building RHS");

    TIC2
    auto p_lhs = make_problem<2>(surface, surface, hyp);
    auto lhs = mesh_to_mesh_operator(p_lhs, qs);
    TOC("Building LHS matrices");
    TIC2
    auto condensed_lhs = condense_block_operator(
        {constraint_matrix, constraint_matrix},
        {constraint_matrix, constraint_matrix},
        lhs
    );
    auto combined_lhs = combine_components(condensed_lhs);
    TOC("Condensing LHS matrices");

    std::cout << "Condition number: " << arma_cond(combined_lhs.ops[0]) << std::endl;

    TIC2
    auto inv_lhs = arma_invert(combined_lhs.ops[0]);

    auto disp_reduced = apply_operator({1, 1, {inv_lhs}}, rhs.data);
    TOC("Solve");

    auto disp_reduced_vec = expand(rhs, disp_reduced);
    std::vector<std::vector<double>> soln{
        distribute_vector(constraint_matrix, disp_reduced_vec[0], surface.n_dofs()),
        distribute_vector(constraint_matrix, disp_reduced_vec[1], surface.n_dofs())
    };

    auto file = HDFOutputter("test_out/planestrain_u.hdf5");
    out_surface(file, surface, soln);
}
