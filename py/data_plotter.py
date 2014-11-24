import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import mayavi.mlab as mlab

x_index = [0, 3, 6]
y_index = [1, 4, 7]
z_index = [2, 5, 8]

def main(filename, values_dim):
    f = h5py.File(filename)
    facets = f['facets']
    data = f['values'][:, values_dim]
    vertices = np.array([
        facets[:, x_index].flatten(),
        facets[:, y_index].flatten(),
        facets[:, z_index].flatten()
    ]).T

    n_v = vertices.shape[0]
    faces = np.arange(n_v).reshape((n_v / 3, 3))

    mlab.triangular_mesh(vertices[:, 0], vertices[:, 1], vertices[:, 2],
                       faces[:,:], scalars = data)
    mlab.show()

if __name__ == "__main__":
    advice = "Usage is: python py/data_plotter.py filename column" +\
             "\n where the column specifies the column of the values" +\
             " dataspace to display"
    if len(sys.argv) < 3:
        print(advice)
        sys.exit()
    filename = sys.argv[1]
    try:
        values_dim = int(sys.argv[2])
    except:
        print(advice)
    main(sys.argv[1], int(sys.argv[2]))
