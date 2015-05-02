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

if __name__ == "__main__":
    test_vector_to_array()
