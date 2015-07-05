from tbempy.interior_meshing import *
import numpy as np

def box_mesh():
    return np.array([
        [[0, 0], [1, 0]],
        [[1, 0], [1, 1]],
        [[1, 1], [0, 1]],
        [[0, 1], [0, 0]]
    ])

def test_determine_extents():
    min_corner, max_corner = determine_extents(box_mesh())
    np.testing.assert_almost_equal(min_corner, [0, 0])
    np.testing.assert_almost_equal(max_corner, [1, 1])

def test_expand_extents():
    min_corner, max_corner = expand_extents(*determine_extents(box_mesh()))
    np.testing.assert_almost_equal(min_corner, [-0.5, -0.5])
    np.testing.assert_almost_equal(max_corner, [1.5, 1.5])

def test_extents_to_box():
    result = extents_to_box_2d(*determine_extents(box_mesh()))
    np.testing.assert_almost_equal(result, box_mesh())

def test_interior_meshing():
    facets = box_mesh()
    create_mesh(facets)
    create_mesh(add_extent_surface(facets))
