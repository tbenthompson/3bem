from tbempy import *
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
    pot_constraints = convert_to_constraints(continuity)
    pot_bc_constraints = interpolate_bc_constraints(
        surface, range(0, switch_dof), theta_u
    )
    pot_constraints.extend(pot_bc_constraints)
    pot_cm = from_constraints(pot_constraints)

    flux_constraints = interpolate_bc_constraints(
        surface, range(switch_dof, n_dofs), theta_dudn
    )
    flux_cm = from_constraints(flux_constraints)

    single_kernel = LaplaceSingle()
    double_kernel = LaplaceDouble()

    single_mthd = make_adaptive_integration_mthd(qs, single_kernel)
    single_op = integral_operator(surface, surface, single_mthd)

    double_mthd = make_adaptive_integration_mthd(qs, double_kernel)
    double_op = integral_operator(surface, surface, double_mthd)

    mass_op = mass_operator_scalar(surface, 2)

    all_zeros = VectorX([0.0] * n_dofs)
    condensed = BlockVectorX([
        condense_vector(pot_cm, all_zeros),
        condense_vector(flux_cm, all_zeros)
    ])
    dof_map = block_dof_map_from_functions(condensed)
    rhs = concatenate(dof_map, condensed)

    def mv(v):
        vec_v = VectorX(v)
        both = expand(dof_map, vec_v)
        pot = distribute_vector(pot_cm, both.storage[0], surface.n_dofs())
        flux = distribute_vector(flux_cm, both.storage[1], surface.n_dofs())
        single_eval = single_op.apply(flux)
        double_eval = double_op.apply(pot)
        mass_eval = mass_op.apply(pot)
        out = mass_eval - single_eval + double_eval
        pot_out = condense_vector(pot_cm, out)
        flux_out = condense_vector(flux_cm, out)
        res = np.array(concatenate(dof_map, BlockVectorX([pot_out, flux_out])).storage)
        res -= mv.rhs
        return res
    mv.rhs = np.array(rhs.storage)

    zeros_mv = mv([0.0] * len(rhs.storage))
    mv.rhs = zeros_mv
    np_rhs = -zeros_mv
    A = sp_la.LinearOperator((np_rhs.shape[0], np_rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, np_rhs, tol = linear_solve_tol)
    both_solns = expand(dof_map, VectorX(res[0]))
    soln_pot = distribute_vector(pot_cm, both_solns.storage[0], surface.n_dofs())
    soln_flux = distribute_vector(flux_cm, both_solns.storage[1], surface.n_dofs())

    vs = surface.facets.reshape((2 * surface.n_facets(), 2))
    xs = vs[:, 0]
    ys = vs[:, 1]
    exact_pot = np.array([theta_u([x, y]) for x, y in zip(xs, ys)])
    est_pot = np.array(soln_pot.storage)
    exact_flux = np.array([theta_dudn([x, y]) for x, y in zip(xs, ys)])
    est_flux = np.array(soln_flux.storage)
    error_flux = np.sqrt(np.mean((est_flux - exact_flux) ** 2))
    error_pot = np.sqrt(np.mean((est_pot - exact_pot) ** 2))
    assert(error_flux < 2e-4)
    assert(error_pot < 2e-6)

if __name__ == "__main__":
    test_corners()
