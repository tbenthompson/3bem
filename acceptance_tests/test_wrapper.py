from tbempy import *
import numpy as np
import pytest

def test_vector_to_array():
    circle = TwoD.circle_mesh([0, 0], 1.0, 2)
    np.testing.assert_almost_equal(
        circle.facets[3,:,:],
        [[0.31622777, 0.9486833], [0.0, 1.0]],
        6
    )

def test_mesh_constructor():
    facets = np.random.rand(4,2,2)
    m = TwoD.Mesh(facets)
    np.testing.assert_almost_equal(facets, m.facets)

def test_mesh_union():
    facets1 = np.random.rand(4,2,2)
    m1 = TwoD.Mesh(facets1)
    facets2 = np.random.rand(4,2,2)
    m2 = TwoD.Mesh(facets2)
    m_union = TwoD.Mesh.create_union([m1, m2])
    assert(m_union.facets.shape[0] == 8)

def test_mix_dims_failure():
    facets1 = np.random.rand(4,2,2)
    m1 = TwoD.Mesh(facets1)
    with pytest.raises(TypeError) as excinfo:
        m_union = ThreeD.Mesh.create_union([m1])

def test_dim():
    assert(TwoD.dim == 2)
    assert(ThreeD.dim == 3)

def test_fmm_config():
    config = TwoD.FMMConfig(0.3, 20, 100, 0.05, False)

if __name__ == "__main__":
    test_mesh_union()
