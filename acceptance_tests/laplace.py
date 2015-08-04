from tbempy import *
import tbempy.TwoD
import tbempy.ThreeD
import scipy.sparse.linalg as sp_la
import numpy as np

def np_matrix_from_tbem_matrix(matrix):
    np_matrix = matrix.data()
    data_shape = np_matrix.shape[0]
    rows = np.sqrt(data_shape)
    cols = rows
    np_matrix = np_matrix.reshape((rows,cols))
    return np_matrix

def solve_direct(tbem, constraint_matrix, matrix, rhs):
    matrix_condensed = tbem.condense_matrix(constraint_matrix, constraint_matrix,
                                       matrix);
    np_matrix = np_matrix_from_tbem_matrix(matrix_condensed)
    return np.linalg.solve(np_matrix, rhs)

def solve_iterative(tbem, constraint_matrix, matrix, rhs):
    def mv(v):
        distributed = tbem.distribute_vector(constraint_matrix, v, matrix.n_rows())
        applied = matrix.apply(distributed)
        condensed = tbem.condense_vector(constraint_matrix, applied)
        return condensed
    A = sp_la.LinearOperator((rhs.shape[0], rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, rhs, tol = 1e-8)
    assert(res[1] == 0) #Check that the iterative solver succeeded
    return res[0]

def solve(dim, mesh, linear_solver, operator_builder, obs_pts, u_fnc, dudn_fnc,
    far_threshold = 3.0):

    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD

    continuity = tbem.mesh_continuity(mesh.begin())
    constraints = tbem.convert_to_constraints(continuity)
    constraint_matrix = tbem.from_constraints(constraints)

    u = tbem.interpolate(mesh, u_fnc)
    dudn = tbem.interpolate(mesh, dudn_fnc)

    double_kernel = tbem.LaplaceDouble()
    double_mthd = tbem.make_adaptive_integrator(
        1e-5, 4, 4, 3, 8, far_threshold, double_kernel
    )
    rhs_double = operator_builder(mesh, mesh, double_mthd, mesh).apply(u)

    rhs_mass = tbem.mass_operator_scalar(mesh, 3).apply(u);

    rhs = rhs_double + rhs_mass
    rhs_condensed = tbem.condense_vector(constraint_matrix, rhs);

    single_kernel = tbem.LaplaceSingle()
    single_mthd = tbem.make_adaptive_integrator(
        1e-5, 4, 4, 3, 8, far_threshold, single_kernel
    )
    matrix = operator_builder(mesh, mesh, single_mthd, mesh);

    soln_condensed = linear_solver(tbem, constraint_matrix, matrix, rhs_condensed)
    soln = tbem.distribute_vector(constraint_matrix, soln_condensed, mesh.n_dofs())
    np_soln = soln
    np_dudn = dudn
    boundary_error = np.sqrt(np.mean((np_soln - np_dudn) ** 2))

    single_int = tbem.dense_interior_operator(
        obs_pts['locs'], obs_pts['normals'], mesh, single_mthd, mesh
    ).apply(soln)
    double_int = tbem.dense_interior_operator(
        obs_pts['locs'], obs_pts['normals'], mesh, double_mthd, mesh
    ).apply(u)

    result = single_int - double_int

    exact = [u_fnc(p) for p in obs_pts['locs']]
    max_interior_error = np.max(np.abs(exact - result))
    return boundary_error, max_interior_error

