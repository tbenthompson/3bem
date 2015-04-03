from tbempy import *
from tbempy.ThreeD import *
from dislocation import *
from antiplane_fault import get_vertices
import numpy as np

def test_okada():
    surf_width = 4
    refine_surf = 3
    fault = rect_mesh(
        [-1, 0, 3], [-1, 0, 0], [1, 0, 0], [1, 0, -3]
    ).refine_repeatedly(refine_surf - 1)

    surface = rect_mesh(
        [-surf_width, -surf_width, 0], [-surf_width, surf_width, 0],
        [surf_width, surf_width, 0], [surf_width, -surf_width, 0]
    ).refine_repeatedly(refine_surf);
    print(surface.n_dofs())

    hyp = ElasticHypersingular(30e9, 0.25)
    soln = solve(3, surface, fault, hyp)

    x = get_vertices(3, surface, 0)
    import matplotlib.pyplot as plt
    plt.plot(x, soln[0])
    plt.show()

if __name__ == "__main__":
    test_okada()
