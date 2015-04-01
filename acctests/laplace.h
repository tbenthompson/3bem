#ifndef __HJKHSDFLJSLLHA_LAPLACE_H
#define __HJKHSDFLJSLLHA_LAPLACE_H
#include "3bem.h"
#include "laplace_kernels.h"
#include "armadillo_facade.h"

using namespace tbem;

double error_inf(const VectorX& a, const VectorX& b);

double random_val(); 

Vec3<double> spherify_pt(Vec3<double> pt, Vec3<double> c, double r); 

Vec3<double> random_pt_sphere(Vec3<double> c, double r);

Vec2<double> random_pt_sphere(Vec2<double> c, double r);

template <size_t dim>
struct LaplaceSoln {
    Mesh<dim> mesh;
    VectorX u;
    VectorX dudn;
};

template <size_t dim>
QuadStrategy<dim> get_qs() 
{
    double far_threshold = 3.0;
    int near_steps = 8;
    int src_quad_pts = 2;
    int obs_quad_pts = 3;
    double tol = 1e-4;
    QuadStrategy<dim> qs(obs_quad_pts, src_quad_pts,
                         near_steps, far_threshold, tol);
    return qs;
}

template <size_t dim, typename Fnc, typename Deriv>
LaplaceSoln<dim> 
dirichlet_laplace_test(const Mesh<dim>& mesh, const Fnc& fnc, const Deriv& deriv) 
{
    auto qs = get_qs<dim>();
    auto continuity = mesh_continuity(mesh.begin());
    auto constraints = convert_to_constraints(continuity);
    auto constraint_matrix = from_constraints(constraints);

    auto u = interpolate<dim>(mesh, fnc);
    auto dudn = interpolate<dim>(mesh, deriv);

    // Construct and evaluate the RHS for a Dirichlet Laplace problem:
    // The integral equation is: DoubleLayer(u) + u = SingleLayer(dudn)
    // First, the double layer potential is evaluated (DoubleLayer(u))
    // This is added to the mass matrix term (u)
    TIC
    LaplaceDouble<dim> double_kernel;
    auto p_double = make_boundary_integral<dim>(mesh, mesh, double_kernel);
    auto rhs_double_op = mesh_to_mesh_operator(p_double, qs);
    auto rhs_double = rhs_double_op.apply({u})[0];

    IdentityScalar<dim> identity;
    auto p_mass = make_boundary_integral<dim>(mesh, mesh, identity);
    auto mass_op = mass_operator(p_mass, qs);
    auto rhs_mass = mass_op.apply({u})[0];
    
    VectorX rhs_full(mesh.n_dofs());
    for (unsigned int i = 0; i < rhs_full.size(); i++){
        rhs_full[i] = rhs_double[i] + rhs_mass[i];
    }

    // Apply the constraints to the RHS vector to get a condensed vector.
    auto rhs = condense_vector(constraint_matrix, rhs_full);

    TOC("RHS Eval")

    // The LHS matrix for a Dirichlet Laplace problem.
    TIC2
    LaplaceSingle<dim> single_kernel;
    auto p_single = make_boundary_integral<dim>(mesh, mesh, single_kernel);
    auto matrix = mesh_to_mesh_operator(p_single, qs);
    TOC("Matrix construct on " + std::to_string(mesh.n_facets()) + " facets");
    TIC2
    auto condensed_op = condense_matrix(
        constraint_matrix, constraint_matrix, matrix.ops[0]);
    TOC("Matrix condense");

    TIC2
    auto inv_condensed_matrix = arma_invert(condensed_op);
    auto dudn_solved_subset = inv_condensed_matrix.apply(rhs);
    TOC("Solve");

    // Get all the constrained DOFs from the reduced DOF vector.
    auto dudn_solved = distribute_vector(constraint_matrix, dudn_solved_subset, mesh.n_dofs());

    // Output the error and the solution 
    std::cout << error_inf(dudn_solved, dudn) << std::endl;

    auto file = HDFOutputter("test_out/laplace" + std::to_string(dim) + "d.hdf5");
    out_surface<dim>(file, mesh, {dudn_solved});

    return {mesh, u, dudn_solved};
}

template <size_t dim, typename Fnc>
void check_laplace_interior(const LaplaceSoln<dim>& soln, 
                            const std::vector<Vec<double,dim>>& pts,
                            Fnc exact_fnc) 
{
    LaplaceSingle<dim> single_kernel;
    LaplaceDouble<dim> double_kernel;

    auto p_single = make_boundary_integral<dim>(soln.mesh, soln.mesh, single_kernel);
    auto p_double = make_boundary_integral<dim>(soln.mesh, soln.mesh, double_kernel);

    auto qs = get_qs<dim>();

    for(int i = 0; i < pts.size(); i++) {
        auto obs_pt = pts[i]; 
        ObsPt<dim> obs = {0.001, obs_pt, zeros<Vec<double,dim>>::make(),
                          zeros<Vec<double,dim>>::make()};
       
        auto double_layer_op = mesh_to_point_operator(p_double, qs, obs);
        double double_layer = double_layer_op.apply({soln.u})[0][0];
        auto single_layer_op = mesh_to_point_operator(p_single, qs, obs);
        double single_layer = single_layer_op.apply({soln.dudn})[0][0];
        double result = single_layer - double_layer;
        double exact = exact_fnc(obs_pt);
        double error = std::fabs(exact - result);
        if (error > 1e-2) {
            std::cout << "Failed with point: " << obs_pt << std::endl;
            std::cout << result << " " << exact << std::endl;
        }
    }
}

#endif
