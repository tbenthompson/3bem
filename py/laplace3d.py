from tbempy import *
from tbempy.ThreeD import *
import numpy as np
from laplace import*
import random
import matplotlib.pyplot as plt


def harmonic_u(x):
    return 1.0 / np.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])

def make_harmonic_dudn(center):
    def harmonic_dudn(x):
        loc = np.array([x[0], x[1], x[2]])
        dist_to_origin = np.linalg.norm(loc)
        n = center - loc
        n /= np.linalg.norm(n)
        return n.dot(loc) / (dist_to_origin ** 3)
    return harmonic_dudn

def random_pt_sphere(center, max_r):
    D = np.array([random.gauss(0, 1) for i in range(3)])
    D /= np.sqrt(sum(D ** 2))
    u = random.uniform(0.0, 1.0)
    r = max_r * (u ** (1.0 / 3.0))
    pt = np.array(center) + r * D
    return pt

def make_interior_pts(n, center, max_r):
    pts = [random_pt_sphere(center, max_r) for i in range(n)]
    inward_normal = center - np.array(pts)
    inward_normal /= np.reshape(np.linalg.norm(inward_normal, axis = 1), (n, 1))
    obs_pts = [ObsPt(0.001, p, n, n) for (p, n) in zip(pts, inward_normal)]
    return pts, obs_pts

def run():
    refine = 3
    center = np.array([5.0, 0.0, 0.0])
    r = 3.0
    obs_radius = 2.9
    qs = QuadStrategy(3, 2, 8, 3.0, 1e-5)

    pts, obs_pts = make_interior_pts(2000, center, obs_radius)
    pts = np.array(pts)
    sphere = sphere_mesh(center, r, refine)
    return solve(3, sphere, solve_direct, obs_pts, harmonic_u, make_harmonic_dudn(center))

def test_laplace3d():
    bdry_error, int_error = run()
    assert(bdry_error < 1e-3)
    assert(int_error < 4e-3)


if __name__ == "__main__":
    run()
