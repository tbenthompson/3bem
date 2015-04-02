from tbempy import *
import numpy as np
import time
import random
import matplotlib.pyplot as plt
import scipy.sparse.linalg as sp_la

def log_u(x, y):
    return np.log(np.sqrt(x * x + y * y))

def theta_u(x, y):
    return np.arctan2(y, x)

def make_log_dudn(center):
    def log_dudn(x, y):
        loc = np.array([x, y])
        dist_to_origin2 = x * x + y * y
        v = center - loc
        v /= np.linalg.norm(v)
        return v.dot(loc) / dist_to_origin2
    return log_dudn

def make_theta_dudn(center):
    def theta_dudn(x, y):
        loc = np.array([x, y])
        n = center - loc
        n /= np.linalg.norm(n)
        dy = 1.0 / (x * (1 + (y * y / (x * x))))
        dx = (-y / x) * dy
        return n.dot([dx, dy])
    return theta_dudn

def random_pt_circle(center, max_r):
    theta = random.uniform(0, 2 * np.pi)
    r = random.uniform(0, max_r)
    x = center[0] + r * np.cos(theta)
    y = center[1] + r * np.sin(theta)
    return [x, y]

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
    A = sp_la.LinearOperator((np_rhs.shape[0], np_rhs.shape[0]), matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, np_rhs, tol = 1e-8)
    assert(res[1] == 0) #Check that the iterative solver succeeded
    return res[0]

def run(solver, refine, u_fnc, dudn_fnc_builder):
    center = [5.0, 0.0]
    r = 3
    obs_radius = 2.9
    n_test_pts = 100
    qs = QuadStrategy(3, 2, 8, 3.0, 1e-5)

    circle = circle_mesh(center, r, refine)
    continuity = mesh_continuity(circle.begin())
    constraints = convert_to_constraints(continuity)
    constraint_matrix = from_constraints(constraints)

    u = interpolate(circle, u_fnc)
    dudn = interpolate(circle, dudn_fnc_builder(center))

    double_kernel = LaplaceDouble()
    double_layer = make_boundary_integral(circle, circle, double_kernel)
    rhs_double = mesh_to_mesh_operator(double_layer, qs).apply_scalar(u)

    identity = IdentityScalar()
    mass_integral = make_boundary_integral(circle, circle, identity)
    rhs_mass = mass_operator(mass_integral, qs).apply_scalar(u);

    rhs = rhs_double + rhs_mass
    rhs_condensed = condense_vector(constraint_matrix, rhs);

    single_kernel = LaplaceSingle()
    single_layer = make_boundary_integral(circle, circle, single_kernel)
    matrix = mesh_to_mesh_operator(single_layer, qs);
    matrix_condensed = condense_matrix(constraint_matrix, constraint_matrix,
                                       matrix.get_block(0,0));

    soln_condensed = solver(matrix_condensed, rhs_condensed)
    soln = distribute_vector(constraint_matrix, VectorX(soln_condensed), circle.n_dofs())
    np_soln = np.array(soln.storage)
    np_dudn = np.array(dudn.storage)
    max_error = np.max(np.abs(np_soln - np_dudn))

    n = 100
    pts = [random_pt_circle(center, obs_radius) for i in range(n)]
    # plt.plot(np_soln)
    # plt.plot(np_dudn)
    # plt.show()
    return max_error

def test_convergence():
    es = []
    rs = np.arange(2,8)
    hs = 1.0 / (2 ** rs)
    for refine in rs:
        es.append(run(solve_iterative, refine, log_u, make_log_dudn))
    # plt.loglog(hs, es)
    # plt.show()
    loges = np.log(es)
    loghs = np.log(hs)
    dlogh = loghs[1:] - loghs[:-1]
    dlogerr = loges[1:] - loges[:-1]
    rate = dlogerr / dlogh
    for r in rate:
        assert(r > 1.5)

def test_log_u_direct():
    assert(run(solve_direct, 7, log_u, make_log_dudn) < 1e-4)

def test_log_u_iterative():
    assert(run(solve_iterative, 7, log_u, make_log_dudn) < 1e-4)

def test_theta_u_direct():
    assert(run(solve_direct, 7, theta_u, make_theta_dudn) < 1e-4)

def test_theta_u_iterative():
    assert(run(solve_iterative, 7, theta_u, make_theta_dudn) < 1e-4)

if __name__ == "__main__":
    test_log_u_direct()
