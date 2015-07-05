from antiplane_fault import solve_half_space
from tbempy.TwoD import line_mesh, MeshPreprocessor, Mesh
import numpy as np

def test_split():
    refine = 6
    intersect = 1.0
    fault = line_mesh([intersect, -1], [intersect, 0])
    slip = np.ones(fault.n_dofs())
    surface = line_mesh([50, 0.0], [-50, 0.0]).refine_repeatedly(refine)

    # Split the surface
    mp = MeshPreprocessor()
    intersections = mp.find_intersections(surface.facets, fault.facets)
    split_surface = mp.split_facets_at_intersections(surface.facets, intersections)
    assert(split_surface.shape[0] > surface.facets.shape[0])

    surface = Mesh(split_surface)
    soln = solve_half_space(slip, fault, surface)
    xs = surface.facets[:, :, 0].reshape((surface.n_facets() * 2))
    xs -= intersect
    indices = [i for i in range(len(xs)) if 0 < np.abs(xs[i]) < 10]
    xs = xs[indices]
    exact = np.arctan(1.0 / xs) / np.pi
    error = np.sqrt(np.mean((exact - soln[indices]) ** 2))
    assert(error < 0.01)
