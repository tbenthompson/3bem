from tbempy import *
import tbempy.TwoD
import tbempy.ThreeD
import scipy.sparse.linalg as sp_la
import numpy as np

def np_matrix_from_tbem_matrix(matrix):
    np_matrix = np.array(matrix.data())
    data_shape = np_matrix.shape[0]
    rows = np.sqrt(data_shape)
    cols = rows
    np_matrix = np_matrix.reshape((rows,cols))
    return np_matrix

def solve_direct(constraint_matrix, matrix, rhs):
    matrix_condensed = condense_matrix(constraint_matrix, constraint_matrix,
                                       matrix.get_block(0,0));
    np_matrix = np_matrix_from_tbem_matrix(matrix_condensed)
    np_rhs = np.array(rhs.storage)
    return np.linalg.solve(np_matrix, np_rhs)

def solve_iterative(constraint_matrix, matrix, rhs):
    np_rhs = np.array(rhs.storage)
    def mv(v):
        vec_v = VectorX(v)
        distributed = distribute_vector(constraint_matrix, vec_v, matrix.n_total_rows())
        applied = matrix.apply(distributed)
        condensed = condense_vector(constraint_matrix, applied)
        res = np.array(condensed.storage)
        return res
    A = sp_la.LinearOperator((np_rhs.shape[0], np_rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, np_rhs, tol = 1e-8)
    assert(res[1] == 0) #Check that the iterative solver succeeded
    return res[0]

def solve(dim, mesh, linear_solver, operator_builder, obs_pts, u_fnc, dudn_fnc):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    qs = tbem.QuadStrategy(3, 8, 3.0, 1e-5)

    continuity = tbem.mesh_continuity(mesh.begin())
    constraints = tbem.convert_to_constraints(continuity)
    constraint_matrix = from_constraints(constraints)

    u = tbem.interpolate(mesh, u_fnc)
    dudn = tbem.interpolate(mesh, dudn_fnc)

    double_kernel = tbem.LaplaceDouble()
    double_mthd = tbem.make_adaptive_integration_mthd(qs, double_kernel)
    rhs_double = operator_builder(mesh, mesh, double_mthd).apply(u)

    rhs_mass = tbem.mass_operator_scalar(mesh, 3).apply(u);

    rhs = rhs_double + rhs_mass
    rhs_condensed = condense_vector(constraint_matrix, rhs);

    single_kernel = tbem.LaplaceSingle()
    single_mthd = tbem.make_adaptive_integration_mthd(qs, single_kernel)
    matrix = operator_builder(mesh, mesh, single_mthd);

    soln_condensed = linear_solver(constraint_matrix, matrix, rhs_condensed)
    soln = distribute_vector(constraint_matrix, VectorX(soln_condensed), mesh.n_dofs())
    np_soln = np.array(soln.storage)
    np_dudn = np.array(dudn.storage)
    boundary_error = np.sqrt(np.mean((np_soln - np_dudn) ** 2))

    single_int = tbem.mesh_to_points_operator(obs_pts, mesh, single_mthd)\
                    .apply(soln)
    double_int = tbem.mesh_to_points_operator(obs_pts, mesh, double_mthd)\
                    .apply(u)
    result = np.array(single_int.storage) - np.array(double_int.storage)

    exact = [u_fnc(p.loc) for p in obs_pts]
    max_interior_error = np.max(np.abs(exact - result))
    return boundary_error, max_interior_error

