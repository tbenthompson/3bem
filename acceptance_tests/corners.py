from tbempy.TwoD import *
from laplace2d import theta_u, make_theta_dudn
import scipy.sparse.linalg as sp_la
import numpy as np

def test_corners():
    center = [5, 0]
    theta_dudn = make_theta_dudn(center)
    r = 3.0
    refine_level = 8;
    linear_solve_tol = 1e-7;
    switch_dof = 500
    qs = QuadStrategy(4, 10, 3.0, 1e-5);
    surface = circle_mesh(center, r, refine_level);
    n_dofs = surface.n_dofs();
    if switch_dof > n_dofs:
        switch_dof = n_dofs

    continuity = mesh_continuity(surface.begin());
    potential_constraints = convert_to_constraints(continuity)
    potential_bc_constraints = interpolate_bc_constraints(
        surface, range(0, switch_dof), theta_u
    )
    potential_constraints.extend(potential_bc_constraints)
    potential_cm = from_constraints(potential_constraints)

    flux_constraints = interpolate_bc_constraints(
        surface, range(switch_dof, n_dofs), theta_dudn
    )
    flux_cm = from_constraints(flux_constraints)

    single_kernel = LaplaceSingle()
    double_kernel = LaplaceDouble()

    single_mthd = make_adaptive_integration_mthd(qs, single_kernel)
    single_op = dense_integral_operator(surface, surface, single_mthd, surface)

    double_mthd = make_adaptive_integration_mthd(qs, double_kernel)
    double_op = dense_integral_operator(surface, surface, double_mthd, surface)

    mass_op = mass_operator_scalar(surface, 2)

    all_zeros = np.zeros(n_dofs)
    potential_rhs = condense_vector(potential_cm, all_zeros)
    flux_rhs = condense_vector(flux_cm, all_zeros)
    rhs = np.concatenate((potential_rhs, flux_rhs))

    def mv(v):
        potential_soln = v[:potential_rhs.shape[0]]
        flux_soln = v[potential_rhs.shape[0]:]
        potential = distribute_vector(potential_cm, potential_soln, surface.n_dofs())
        flux = distribute_vector(flux_cm, flux_soln, surface.n_dofs())
        single_eval = single_op.apply(flux)
        double_eval = double_op.apply(potential)
        mass_eval = mass_op.apply(potential)
        out = mass_eval - single_eval + double_eval
        potential_out = condense_vector(potential_cm, out)
        flux_out = condense_vector(flux_cm, out)
        res = np.concatenate((potential_out, flux_out))
        res -= mv.rhs
        # print("iteration " + str(mv.it))
        mv.it += 1
        return res
    mv.rhs = rhs
    mv.it = 0

    zeros_mv = mv(np.zeros_like(rhs))
    mv.rhs = zeros_mv
    np_rhs = -zeros_mv
    A = sp_la.LinearOperator((np_rhs.shape[0], np_rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, np_rhs, tol = linear_solve_tol)
    reduced_potential_soln = res[0][:potential_rhs.shape[0]]
    reduced_flux_soln = res[0][potential_rhs.shape[0]:]
    soln_potential = distribute_vector(
        potential_cm, reduced_potential_soln, surface.n_dofs()
    )
    soln_flux = distribute_vector(flux_cm, reduced_flux_soln, surface.n_dofs())

    vs = surface.facets.reshape((2 * surface.n_facets(), 2))
    xs = vs[:, 0]
    ys = vs[:, 1]
    exact_potential = np.array([theta_u([x, y]) for x, y in zip(xs, ys)])
    exact_flux = np.array([theta_dudn([x, y]) for x, y in zip(xs, ys)])
    error_flux = np.sqrt(np.mean((soln_flux - exact_flux) ** 2))
    error_potential = np.sqrt(np.mean((soln_potential - exact_potential) ** 2))
    assert(error_flux < 2e-4)
    assert(error_potential < 2e-6)

if __name__ == "__main__":
    test_corners()
