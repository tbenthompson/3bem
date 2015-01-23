#ifndef __HJKHSDFLJSLLHA_LAPLACE_H
#define __HJKHSDFLJSLLHA_LAPLACE_H
#include "3bem.h"
#include "laplace_kernels.h"
#include "armadillo_interface.h"

using namespace tbem;

double error_inf(const std::vector<double>& a, const std::vector<double>& b);

double random_val(); 

Vec3<double> spherify_pt(Vec3<double> pt, Vec3<double> c, double r); 

Vec3<double> random_pt_sphere(Vec3<double> c, double r);

Vec2<double> random_pt_sphere(Vec2<double> c, double r);

template <size_t dim, typename Fnc, typename Deriv>
void dirichlet_laplace_test(const Mesh<dim>& mesh,
                      std::vector<Vec<double,dim>> test_interior_pts,
                      const Fnc& fnc,
                      const Deriv& deriv) {
    double far_threshold = 3.0;
    int near_steps = 8;
    int src_quad_pts = 2;
    int obs_quad_pts = 3;
    double tol = 1e-4;
    QuadStrategy<dim> qs(obs_quad_pts, src_quad_pts,
                         near_steps, far_threshold, tol);

    auto continuity = mesh_continuity(mesh.begin());
    auto constraints = convert_to_constraints(continuity);
    auto constraint_matrix = from_constraints(constraints);

    auto u = constrained_interpolate<dim>(mesh, fnc, constraint_matrix);
    auto dudn = constrained_interpolate<dim>(mesh, deriv, constraint_matrix);

    // Construct and evaluate the RHS for a Dirichlet Laplace problem:
    // The integral equation is: DoubleLayer(u) + u = SingleLayer(dudn)
    // First, the double layer potential is evaluated (DoubleLayer(u))
    // This is added to the mass matrix term (u)
    TIC
    LaplaceDouble<dim> double_kernel;
    auto p_double = make_problem<dim>(mesh, mesh, double_kernel);
    auto rhs_double_op = mesh_to_mesh_operator(p_double, qs);
    auto rhs_double = apply_operator(rhs_double_op, u);

    IdentityScalar<dim> identity;
    auto p_mass = make_problem<dim>(mesh, mesh, identity);
    auto mass_op = mass_operator(p_mass, qs);
    auto rhs_mass = apply_operator(mass_op, u);
    
    std::vector<double> rhs_full(mesh.n_dofs());
    double mass_factor[2] = {1.0, 1.0};
    for (unsigned int i = 0; i < rhs_full.size(); i++){
        rhs_full[i] = rhs_double[i] + mass_factor[dim - 2] * rhs_mass[i];
    }
    TOC("RHS Eval")

    // The LHS matrix for a Dirichlet Laplace problem.
    TIC2
    LaplaceSingle<dim> single_kernel;
    auto p_single = make_problem<dim>(mesh, mesh, single_kernel);
    auto matrix = mesh_to_mesh_operator(p_single, qs);
    TOC("Matrix construct on " + std::to_string(mesh.n_facets()) + " facets");
    TIC2
    auto condensed_matrix = condense_matrix(constraint_matrix, 
        {matrix.n_rows, matrix.n_cols, matrix.data[0]});
    MatrixOperator condensed_op{
        condensed_matrix.n_rows, condensed_matrix.n_cols, 1, 1, {condensed_matrix.data}
    };
    TOC("Matrix condense");

    TIC2
    auto inv_condensed_matrix = armadillo_invert(condensed_op.data[0]);
    MatrixOperator inv_condensed_op{
        condensed_matrix.n_rows, condensed_matrix.n_cols, 1, 1, {inv_condensed_matrix}
    };
    std::cout << "DOFs: " << condensed_matrix.n_rows << std::endl;
    TOC("Invert");


    // Apply the constraints to the RHS vector to get a condensed vector.
    auto rhs = condense_vector(constraint_matrix, rhs_full);
    
    // Actually solve the linear system by providing a function to evaluate
    // matrix vector products.
    // int count = 0;
    // double linear_solve_tol = 1e-5;
    // auto dudn_solved_subset = solve_system(rhs, linear_solve_tol,
    //     [&] (std::vector<double>& x, std::vector<double>& y) {
    //         std::cout << "iteration " << count << std::endl;
    //         count++;
    //         TIC
    //         auto y_mult = apply_operator(condensed_op, x); 
    //         std::copy(y_mult.begin(), y_mult.end(), y.begin());
    //         TOC("Matrix multiply on " + std::to_string(mesh.n_facets()) + " faces");
    //     });
    auto dudn_solved_subset = apply_operator(inv_condensed_op, rhs);
    // Get all the constrained DOFs from the reduced DOF vector.
    auto dudn_solved = distribute_vector(constraint_matrix, dudn_solved_subset, mesh.n_dofs());

    // Output the error and the solution 
    std::cout << error_inf(dudn_solved, dudn) << std::endl;

    auto file = HDFOutputter("test_out/laplace" + std::to_string(dim) + "d.hdf5");
    out_surface<dim>(file, mesh, {dudn_solved});

    for(int i = 0; i < test_interior_pts.size(); i++) {
        auto obs_pt = test_interior_pts[i]; 
        ObsPt<dim> obs = {0.001, obs_pt, zeros<Vec<double,dim>>::make(),
                          zeros<Vec<double,dim>>::make()};
       
        auto double_layer_op = mesh_to_point_operator(p_double, qs, obs);
        double double_layer = apply_operator(double_layer_op, u)[0];
        auto single_layer_op = mesh_to_point_operator(p_single, qs, obs);
        double single_layer = apply_operator(single_layer_op, dudn)[0];
        double result = single_layer - double_layer;
        double exact = fnc(obs_pt);
        double error = std::fabs(exact - result);
        if (error > 1e-2) {
            std::cout << "Failed with point: " << obs_pt << std::endl;
            std::cout << result << " " << exact << std::endl;
        }
    }
}

#endif
