import h5py
import numpy as np

def get_2d_vertices(facets):
    x_index = [0, 2]
    y_index = [1, 3]
    vertices = np.array([
        facets[:, x_index].flatten(),
        facets[:, y_index].flatten()
    ]).T
    return vertices

def get_3d_vertices(facets):
    x_index = [0, 3, 6]
    y_index = [1, 4, 7]
    z_index = [2, 5, 8]
    vertices = np.array([
        facets[:, x_index].flatten(),
        facets[:, y_index].flatten(),
        facets[:, z_index].flatten()
    ]).T
    return vertices

def dataset_name(i):
    return "values" + str(i)

def load_data(f):
    datasets = []
    i = 0
    while dataset_name(i) in f:
        datasets.append(f[dataset_name(i)][:])
        i += 1
    assert(len(datasets) > 0)
    data = np.hstack(datasets)
    return data

def load_surface(filename):
    f = h5py.File(filename)
    facets = f['locations']
    if (facets.shape[1] == 9):
        vertices = get_3d_vertices(facets)
    elif (facets.shape[1] == 4):
        vertices = get_2d_vertices(facets)
    data = load_data(f)
    return vertices, data

def load_volume(filename):
    f = h5py.File(filename)
    locs = f['locations']
    data = load_data(f)
    return locs, data
