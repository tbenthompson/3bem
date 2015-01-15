#include "3bem.h"
#include "elastic_kernels.h"
#include "disloc_shared.h"

using namespace tbem;

int main() {
    double surf_width = 4;
    int refine_surf = 4;
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

    std::vector<double> dux(fault.n_dofs(), 1.0);
    std::vector<double> duyz(fault.n_dofs(), 0.0);
    std::vector<std::vector<double>> du{dux, duyz, duyz};

    TIC
    auto p_rhs = make_problem<3>(fault, surface, hyp);
    auto rhs_op = mesh_to_mesh_operator(p_rhs, qs);
    auto all_dofs_rhs = apply_operator(rhs_op, du);
    auto rhs = concatenate({
        constraint_matrix.get_reduced(all_dofs_rhs[0]),
        constraint_matrix.get_reduced(all_dofs_rhs[1]),
        constraint_matrix.get_reduced(all_dofs_rhs[2])
    });
    TOC("Building RHS");

    TIC2
    auto p_lhs = make_problem<3>(surface, surface, hyp);
    auto lhs = mesh_to_mesh_operator(p_lhs, qs);
    TOC("Building LHS matrices");

    int count = 0;
    auto disp_reduced = solve_system(rhs.data, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            auto x_vec_reduced = expand(rhs, x);
            std::vector<std::vector<double>> x_vec{
                constraint_matrix.get_all(x_vec_reduced[0], surface.n_dofs()),
                constraint_matrix.get_all(x_vec_reduced[1], surface.n_dofs()),
                constraint_matrix.get_all(x_vec_reduced[2], surface.n_dofs())
            };
            auto y_vec = apply_operator(lhs, x_vec);
            auto out = concatenate({
                constraint_matrix.get_reduced(y_vec[0]),
                constraint_matrix.get_reduced(y_vec[1]),
                constraint_matrix.get_reduced(y_vec[2])
            });
            for (std::size_t i = 0; i < out.data.size(); i++) {
                y[i] = -out.data[i];
            }
        }
    );

    auto disp_reduced_vec = expand(rhs, disp_reduced);
    std::vector<std::vector<double>> soln{
        constraint_matrix.get_all(disp_reduced_vec[0], surface.n_dofs()),
        constraint_matrix.get_all(disp_reduced_vec[1], surface.n_dofs()),
        constraint_matrix.get_all(disp_reduced_vec[2], surface.n_dofs())
    };

    auto filex = HDFOutputter("test_out/rect_dislocation_ux.hdf5");
    out_surface(filex, surface, soln[0], 1);
    auto filey = HDFOutputter("test_out/rect_dislocation_uy.hdf5");
    out_surface(filey, surface, soln[1], 1);
    auto filez = HDFOutputter("test_out/rect_dislocation_uz.hdf5");
    out_surface(filez, surface, soln[2], 1);
}
