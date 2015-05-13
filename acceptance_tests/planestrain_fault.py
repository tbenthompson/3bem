from tbempy import *
from tbempy.TwoD import *
from dislocation import *
import numpy as np

def exact_displacements(x):
    s = np.sqrt(2)
    delta = 3 * np.pi / 4
    d = 1
    xd = d / np.tan(delta)
    xsi = (x - xd) / d
    exact_ux = (-s / np.pi) * (
            np.cos(delta) * (np.arctan(xsi) - (np.pi / 2) * np.sign(x)) +
            (np.sin(delta) - xsi * np.cos(delta)) / (1 + xsi ** 2))
    exact_uy = (s / np.pi) * (
            np.sin(delta) * (np.arctan(xsi) - (np.pi / 2) * np.sign(x)) +
            (np.cos(delta) + xsi * np.sin(delta)) / (1 + xsi ** 2))
    return exact_ux, exact_uy

def check_planestrain_error(surface, soln):
    xs = surface.facets[:, :, 0].reshape((surface.n_facets() * 2))
    indices = [i for i in range(xs.shape[0])
               if 0 < np.abs(xs[i]) < 10]
    exact_ux, exact_uy = exact_displacements(xs[indices])
    ux = soln[0][indices]
    uy = soln[1][indices]
    ux_error = np.sqrt(np.mean((ux - exact_ux) ** 2))
    uy_error = np.sqrt(np.mean((uy - exact_uy) ** 2))
    assert(ux_error < 2e-3)
    assert(uy_error < 6e-3)

def test_planestrain():
    fault = line_mesh([-1, -1], [0, 0]).refine_repeatedly(2)
    surface = line_mesh([100, 0], [-100, 0]).refine_repeatedly(9)
    slip = np.ones(2 * fault.n_dofs())

    qs = QuadStrategy(3, 8, 3.0, 1e-4)

    hyp = ElasticHypersingular(30e9, 0.25)
    soln = solve(2, surface, fault, hyp, qs, slip)
    check_planestrain_error(surface, soln)



if __name__ == "__main__":
    test_planestrain()
