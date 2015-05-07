import numpy as np
import random
from tbempy.TwoD import *
from laplace import *

def log_u(x):
    return np.log(np.sqrt(x[0] * x[0] + x[1] * x[1]))

def theta_u(x):
    return np.arctan2(x[1], x[0])

def make_log_dudn(center):
    def log_dudn(x):
        loc = np.array([x[0], x[1]])
        dist_to_origin2 = x[0] * x[0] + x[1] * x[1]
        n = center - loc
        n /= np.linalg.norm(n)
        return n.dot(loc) / dist_to_origin2
    return log_dudn

def make_theta_dudn(center):
    def theta_dudn(x):
        loc = np.array([x[0], x[1]])
        n = center - loc
        n /= np.linalg.norm(n)
        dy = 1.0 / (x[0] * (1 + (x[1] * x[1] / (x[0] * x[0]))))
        dx = (-x[1] / x[0]) * dy
        return n.dot([dx, dy])
    return theta_dudn

def random_pt_circle(center, max_r):
    theta = random.uniform(0, 2 * np.pi)
    r_integral = random.uniform(0, 0.5 * max_r * max_r)
    r = np.sqrt(2.0 * r_integral)
    x = center[0] + r * np.cos(theta)
    y = center[1] + r * np.sin(theta)
    return [x, y]

def make_interior_pts(n, center, max_r):
    pts = [random_pt_circle(center, max_r) for i in range(n)]
    inward_normal = center - np.array(pts)
    inward_normal /= np.reshape(np.linalg.norm(inward_normal, axis = 1), (n, 1))
    obs_pts = [ObsPt(0.001, p, n, n) for (p, n) in zip(pts, inward_normal)]
    return pts, obs_pts

def run(linear_solver, operator_builder, refine, u_fnc, dudn_fnc_builder):
    center = [5.0, 0.0]
    r = 3
    obs_radius = 2.9
    n_test_pts = 100

    pts, obs_pts = make_interior_pts(100, center, obs_radius)
    circle = circle_mesh(center, r, refine)
    return solve(2, circle, linear_solver, operator_builder,
                 obs_pts, u_fnc, dudn_fnc_builder(center))

def test_convergence():
    es = []
    rs = np.arange(2,8)
    hs = 1.0 / (2 ** rs)
    for refine in rs:
        es.append(run(solve_iterative, integral_operator,
                      refine, log_u, make_log_dudn)[0])
    # plt.loglog(hs, es)
    # plt.show()
    loges = np.log(es)
    loghs = np.log(hs)
    dlogh = loghs[1:] - loghs[:-1]
    dlogerr = loges[1:] - loges[:-1]
    rate = dlogerr / dlogh
    for r in rate:
        assert(r > 1.85)

def test_log_u():
    for mthd in [(solve_direct, dense_integral_operator),
                 (solve_iterative, integral_operator)]:
        bdry_error, int_error = run(mthd[0], mthd[1], 7, log_u, make_log_dudn)
        assert(bdry_error < 1e-5)
        assert(int_error < 2e-5)

def test_theta_u():
    for mthd in [(solve_direct, dense_integral_operator),
                 (solve_iterative, integral_operator)]:
        bdry_error, int_error = run(mthd[0], mthd[1], 7, theta_u, make_theta_dudn)
        assert(bdry_error < 1e-5)
        assert(int_error < 2e-5)

if __name__ == "__main__":
    test_log_u()
