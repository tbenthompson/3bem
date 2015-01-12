import numpy as np
from matplotlib import pyplot as plt
from hdf_load import load_surface, load_volume
from testing import check_error_norm

def plot_volumetric(locs, data):
    n_pts = int(np.sqrt(locs.shape[0]))
    x = locs[:, 0].reshape(n_pts, n_pts)
    y = locs[:, 1].reshape(n_pts, n_pts)
    d = data.reshape(n_pts, n_pts)
    opts = dict()
    # opts['levels'] = np.linspace(-0.5, 0.5, 20)
    opts['levels'] = np.linspace(15, 28, 20)
    plt.figure()
    plt.contour(x, y, d,
                      colors = '#333333',
                      linestyles = 'solid',
                      **opts)
    plt.contourf(x, y, d, **opts)
    plt.colorbar()

shear_modulus = 30e9
d = 1.0
s = 1.0

def test_fullspace():
    filename = 'test_out/antiplane_full_space.hdf5'
    vertices, data = load_surface(filename)
    indices = [i for i in range(vertices.shape[0]) if np.abs(vertices[i, 0]) > 0]
    x = vertices[indices, 0]
    uz = data[indices, 0]
    exact_uz = s * np.arctan((d / 2) / x) / np.pi
    np.testing.assert_almost_equal(uz, exact_uz, 5)

def test_halfspace():
    filename = 'test_out/antiplane_half_space.hdf5'
    vertices, data = load_surface(filename)
    indices = [i for i in range(vertices.shape[0])
               if 0 < np.abs(vertices[i, 0]) < 10]
    x = vertices[indices, 0]
    uz = data[indices, 0]
    exact_uz = s * np.arctan(d / x) / np.pi
    np.testing.assert_almost_equal(uz, exact_uz, 1)

def test_halfspace_vol_disp():
    filename = 'test_out/antiplane_half_space_volu.hdf5'
    locs, data = load_volume(filename)
    x = locs[:, 0]
    y = locs[:, 1]
    exact_uz = (-s / (2 * np.pi)) * (
        np.arctan((y - d) / x) -
        np.arctan((y + d) / x))
    np.testing.assert_almost_equal(data[:, 0], exact_uz, 2)

def test_halfspace_vol_trac():
    tracx_filename = 'test_out/antiplane_half_space_voltx.hdf5'
    tracy_filename = 'test_out/antiplane_half_space_volty.hdf5'
    locs, tracx = load_volume(tracx_filename)
    locs, tracy = load_volume(tracy_filename)
    tracx = tracx[:, 0]
    tracy = tracy[:, 0]
    x = locs[:, 0]
    y = locs[:, 1]
    exact_tracx = -(s * shear_modulus) / (2 * np.pi) * (
        ((y + d) / (x ** 2 + (y + d) ** 2)) -
        ((y - d) / (x ** 2 + (y - d) ** 2)))
    exact_tracy = (s * shear_modulus) / (2 * np.pi) * (
        (x / (x ** 2 + (y + d) ** 2)) -
        (x / (x ** 2 + (y - d) ** 2)))

    check_error_norm(tracx, exact_tracx, 0.01)
    check_error_norm(tracy, exact_tracy, 0.01)
