#include "disloc_shared.h"

using namespace tbem;

int main() {
    // HALF SPACE THRUST FAULT IN PLANE STRAIN

    // Fault mesh.
    auto fault = line_mesh({-1, -1}, {0, 0}).refine_repeatedly(0);

    // Simple low accuracy quadrature strategy.
    QuadStrategy<2> qs(3);

    // Earth's surface
    auto surface = line_mesh({-100, 0.0}, {100, 0.0}).refine_repeatedly(11);

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
        constraint_matrix.get_reduced(all_dofs_rhs[0]),
        constraint_matrix.get_reduced(all_dofs_rhs[1])
    });
    TOC("Building RHS");

    TIC2
    auto p_lhs = make_problem<2>(surface, surface, hyp);
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
                constraint_matrix.get_all(x_vec_reduced[1], surface.n_dofs())
            };
            auto y_vec = apply_operator(lhs, x_vec);
            auto out = concatenate({
                constraint_matrix.get_reduced(y_vec[0]),
                constraint_matrix.get_reduced(y_vec[1])
            });
            for (std::size_t i = 0; i < out.data.size(); i++) {
                y[i] = out.data[i];
            }
        }
    );

    auto disp_reduced_vec = expand(rhs, disp_reduced);
    std::vector<std::vector<double>> soln{
        constraint_matrix.get_all(disp_reduced_vec[0], surface.n_dofs()),
        constraint_matrix.get_all(disp_reduced_vec[1], surface.n_dofs())
    };

    auto filex = HDFOutputter("test_out/planestrain_ux.hdf5");
    out_surface(filex, surface, soln[0], 1);
    auto filey = HDFOutputter("test_out/planestrain_uy.hdf5");
    out_surface(filey, surface, soln[1], 1);
}
