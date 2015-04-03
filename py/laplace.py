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

def solve_direct(matrix, rhs):
    np_matrix = np_matrix_from_tbem_matrix(matrix)
    np_rhs = np.array(rhs.storage)
    return np.linalg.solve(np_matrix, np_rhs)

def solve_iterative(matrix, rhs):
    np_rhs = np.array(rhs.storage)
    def mv(v):
        return np.array(matrix.apply(VectorX(v)).storage)
    A = sp_la.LinearOperator((np_rhs.shape[0], np_rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, np_rhs, tol = 1e-8)
    assert(res[1] == 0) #Check that the iterative solver succeeded
    return res[0]

def solve(dim, mesh, linear_solver, obs_pts, u_fnc, dudn_fnc):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    qs = tbem.QuadStrategy(3, 2, 8, 3.0, 1e-5)

    continuity = tbem.mesh_continuity(mesh.begin())
    constraints = tbem.convert_to_constraints(continuity)
    constraint_matrix = from_constraints(constraints)

    u = tbem.interpolate(mesh, u_fnc)
    dudn = tbem.interpolate(mesh, dudn_fnc)

    double_kernel = tbem.LaplaceDouble()
    double_layer = tbem.make_boundary_integral(mesh, mesh, double_kernel)
    rhs_double = tbem.mesh_to_mesh_operator(double_layer, qs).apply_scalar(u)

    identity = tbem.IdentityScalar()
    mass_integral = tbem.make_boundary_integral(mesh, mesh, identity)
    rhs_mass = tbem.mass_operator(mass_integral, qs).apply_scalar(u);

    rhs = rhs_double + rhs_mass
    rhs_condensed = condense_vector(constraint_matrix, rhs);

    single_kernel = tbem.LaplaceSingle()
    single_layer = tbem.make_boundary_integral(mesh, mesh, single_kernel)
    matrix = tbem.mesh_to_mesh_operator(single_layer, qs);
    matrix_condensed = condense_matrix(constraint_matrix, constraint_matrix,
                                       matrix.get_block(0,0));

    soln_condensed = linear_solver(matrix_condensed, rhs_condensed)
    soln = distribute_vector(constraint_matrix, VectorX(soln_condensed), mesh.n_dofs())
    np_soln = np.array(soln.storage)
    np_dudn = np.array(dudn.storage)
    boundary_error = np.sqrt(np.mean((np_soln - np_dudn) ** 2))

    single_int = tbem.mesh_to_points_operator(single_layer, qs, obs_pts)\
                    .apply_scalar(soln)
    double_int = tbem.mesh_to_points_operator(double_layer, qs, obs_pts)\
                    .apply_scalar(u)
    result = np.array(single_int.storage) - np.array(double_int.storage)

    exact = [u_fnc(p.loc) for p in obs_pts]
    max_interior_error = np.max(np.abs(exact - result))
    return boundary_error, max_interior_error

