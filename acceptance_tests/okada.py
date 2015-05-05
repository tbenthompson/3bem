from tbempy import *
from tbempy.ThreeD import *
from dislocation import *
import numpy as np
from okada_wrapper import dc3dwrapper

def compute_okada(mu, pr, vertices):
    n_pts = vertices.shape[0]
    disp = np.empty((n_pts, 3))
    for i in range(n_pts):
        m = 30e9
        pr = 0.25
        l = (2 * m * pr) / (1 - 2 * pr);
        alpha = (l+m) / (l + 2 * m)
        v = vertices[i,:]
        success, u, grad_u = dc3dwrapper(alpha, v, 2.0, 90, [-1.0, 1.0],
                                         [-1.0, 2.0], [-1.0, 0.0, 0.0])
        if success != 0:
            pass
        disp[i, :] = u
    return disp

def test_okada():
    surf_width = 4
    refine_surf = 3
    fault = rect_mesh(
        [-1, 0, -3.0], [-1, 0, -0.0],
        [1, 0, -0.0], [1, 0, -3.0]
    ).refine_repeatedly(refine_surf - 1)
    slip = np.concatenate((
        np.ones(fault.n_dofs()),
        np.zeros(fault.n_dofs()),
        np.zeros(fault.n_dofs())
    ))
    qs = QuadStrategy(2, 5, 3.0, 1e-3)

    surface = rect_mesh(
        [-surf_width, -surf_width, 0], [-surf_width, surf_width, 0],
        [surf_width, surf_width, 0], [surf_width, -surf_width, 0]
    ).refine_repeatedly(refine_surf);

    mu = 30e9
    pr = 0.25
    hyp = ElasticHypersingular(mu, pr)
    soln = solve(3, surface, fault, hyp, qs, slip)
    np.save('okada.dat', soln)
    # soln = np.load('okada.dat')

    vs = surface.facets.reshape((surface.n_facets() * 3, 3))
    facets = np.arange(vs.shape[0]).reshape((surface.n_facets(), 3))

    exact = compute_okada(mu, pr, vs)

    vmax = 0.2
    opts = dict(shading = 'gouraud', vmin = -vmax, vmax = vmax)
    # plt.figure()
    # plt.tripcolor(vs[:, 0], vs[:, 1], facets, soln[0], **opts)
    # plt.title('est')
    # plt.figure()
    # plt.tripcolor(vs[:, 0], vs[:, 1], facets, exact[:, 0], **opts)
    # plt.title('exact')
    # plt.show()

if __name__ == "__main__":
    test_okada()
