#include "3bem.h"
#include "elastic_kernels.h"
#include "disloc_shared.h"

using namespace tbem;

int main() {
    double surf_width = 4;
    int refine_surf = 3;
    double far_threshold = 3.0;
    int near_steps = 5;
    int src_quad_pts = 2;
    int obs_quad_pts = 2;
    double near_tol = 1e-3;

    auto fault = rect_mesh(
        {-1, 0, -3.0}, {-1, 0, -0.0},
        {1, 0, -0.0}, {1, 0, -3.0}
    ).refine_repeatedly(refine_surf - 1);

    auto surface = rect_mesh(
        {-surf_width, -surf_width, 0}, {-surf_width, surf_width, 0},
        {surf_width, surf_width, 0}, {surf_width, -surf_width, 0}
    ).refine_repeatedly(refine_surf);

    auto constraint_matrix = surf_fault_constraints(surface.begin(), fault.begin());

    QuadStrategy<3> qs(obs_quad_pts, src_quad_pts,
                    near_steps, far_threshold, near_tol);

    ElasticHypersingular<3> hyp(30e9, 0.25);
    
    std::cout << "Number of surface DOFs: " << surface.n_dofs() << std::endl;

    VectorX dux(fault.n_dofs(), 1.0);
    VectorX duyz(fault.n_dofs(), 0.0);
    BlockVectorX du{dux, duyz, duyz};

    TIC
    auto p_rhs = make_problem<3>(surface, fault, hyp);
    auto rhs_op = mesh_to_mesh_operator(p_rhs, qs);
    auto all_dofs_rhs = apply_operator(rhs_op, du);
    BlockVectorX condensed{
        condense_vector(constraint_matrix, all_dofs_rhs[0]),
        condense_vector(constraint_matrix, all_dofs_rhs[1]),
        condense_vector(constraint_matrix, all_dofs_rhs[2])
    };
    auto dof_map = block_dof_map_from_functions(condensed);
    auto rhs = concatenate(dof_map, condensed);
    TOC("Building RHS");

    TIC2
    auto p_lhs = make_problem<3>(surface, surface, hyp);
    auto lhs = mesh_to_mesh_operator(p_lhs, qs);
    TOC("Building LHS matrices");

    auto disp_reduced = solve(lhs, rhs, dof_map, surface, constraint_matrix);

    auto disp_reduced_vec = expand(dof_map, disp_reduced);
    BlockVectorX soln{
        distribute_vector(constraint_matrix, disp_reduced_vec[0], surface.n_dofs()),
        distribute_vector(constraint_matrix, disp_reduced_vec[1], surface.n_dofs()),
        distribute_vector(constraint_matrix, disp_reduced_vec[2], surface.n_dofs())
    };

    auto file = HDFOutputter("test_out/rect_dislocation_u.hdf5");
    out_surface(file, surface, soln);
}
