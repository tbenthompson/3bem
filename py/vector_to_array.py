from tbempy import *
from tbempy.TwoD import *
import matplotlib.pyplot as plt

def test_vector_to_array():
    circle = circle_mesh([0, 0], 1.0, 2)
    for i in range(circle.facets.shape[0]):
        f = circle.facets[i, :, :]
        plt.plot([f[0][0], f[1][0]], [f[0][1], f[1][1]])
    plt.show()

if __name__ == "__main__":
    test_vector_to_array()
