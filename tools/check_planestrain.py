import numpy as np
from matplotlib import pyplot as plt
from hdf_load import load_surface, load_volume
from testing import check_error_norm

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

def test_planestrain():
    filename_u = "test_out/planestrain_u.hdf5"
    vertices, data_u = load_surface(filename_u)

    indices = [i for i in range(vertices.shape[0])
               if 0 < np.abs(vertices[i, 0]) < 10]
    x = vertices[indices, 0]
    ux = data_u[indices, 0]
    uy = data_u[indices, 1]

    exact_ux, exact_uy = exact_displacements(x)

    np.testing.assert_almost_equal(ux, exact_ux, 1)
    np.testing.assert_almost_equal(uy, exact_uy, 1)

    check_error_norm(ux, exact_ux, 0.06)
    check_error_norm(uy, exact_uy, 0.06)

