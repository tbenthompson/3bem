#include "3bem.h"
#include "elastic_kernels.h"

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

    auto continuity = mesh_continuity(surface.begin());
    auto cut_cont = cut_at_intersection(continuity, surface.begin(), fault.begin());
    auto constraints = convert_to_constraints(cut_cont);
    auto constraint_matrix = ConstraintMatrix::from_constraints(constraints);

    QuadStrategy<3> qs(obs_quad_pts, src_quad_pts,
                    near_steps, far_threshold, near_tol);

    ElasticHypersingular<3> hyp(30e9, 0.25);
    
    std::cout << "Number of surface DOFs: " << surface.n_dofs() << std::endl;

    std::vector<Vec3<double>> du(fault.n_dofs(), {1.0, 0.0, 0.0});
    std::vector<Vec3<double>> all_dofs_rhs(surface.n_dofs(), 
                                           zeros<Vec3<double>>::make());

    TIC
    auto p_rhs = make_problem<3>(fault, surface, hyp);
    auto rhs_op = interact_matrix(p_rhs, qs);
    auto res = bem_mat_mult(rhs_op, du);
    for (unsigned int i = 0; i < res.size(); i++) {
        all_dofs_rhs[i] += res[i];
    }
    auto rhs = constraint_matrix.get_reduced(all_dofs_rhs);
    TOC("Building RHS");

    TIC2
    auto p_lhs = make_problem<3>(surface, surface, hyp);
    auto lhs = interact_matrix(p_lhs, qs);
    TOC("Building LHS matrices");

    int count = 0;
    auto disp_reduced = solve_system(reinterpret_vector<double>(rhs), 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            auto x_vec_reduced = reinterpret_vector<Vec3<double>>(x);
            auto x_vec = constraint_matrix.get_all(x_vec_reduced, surface.n_dofs());
            auto y_vec = bem_mat_mult(lhs, x_vec);
            auto y_vec_reduced = constraint_matrix.get_reduced(y_vec);
            for (std::size_t i = 0; i < y_vec_reduced.size(); i++) {
                y[3 * i] = -y_vec_reduced[i][0];
                y[3 * i + 1] = -y_vec_reduced[i][1];
                y[3 * i + 2] = -y_vec_reduced[i][2];
            }
        }
    );

    auto disp_reduced_vec = reinterpret_vector<Vec3<double>>(disp_reduced);
    auto disp_vec = constraint_matrix.get_all(disp_reduced_vec, surface.n_dofs());

    auto file = HDFOutputter("test_out/rect_dislocation.hdf5");
    out_surface<3>(file, surface, disp_vec, 3);
}
