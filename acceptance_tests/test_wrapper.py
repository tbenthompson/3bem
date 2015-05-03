from tbempy import *
from tbempy.TwoD import *
import numpy as np

def test_vector_to_array():
    circle = circle_mesh([0, 0], 1.0, 2)
    np.testing.assert_almost_equal(
        circle.facets[3,:,:],
        [[0.31622777, 0.9486833], [0.0, 1.0]],
        6
    )

def test_mesh_constructor():
    facets = np.random.rand(4,2,2)
    m = Mesh(facets)
    np.testing.assert_almost_equal(facets, m.facets)

def test_mesh_union():
    facets1 = np.random.rand(4,2,2)
    m1 = Mesh(facets1)
    facets2 = np.random.rand(4,2,2)
    m2 = Mesh(facets2)
    m_union = Mesh.create_union([m1, m2])
    assert(m_union.facets.shape[0] == 8)

def test_dim():
    import tbempy.TwoD
    import tbempy.ThreeD
    assert(tbempy.TwoD.dim == 2)
    assert(tbempy.ThreeD.dim == 3)

if __name__ == "__main__":
    test_mesh_union()
