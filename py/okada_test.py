import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import mayavi.mlab as mlab
from matplotlib.tri import Triangulation
from okada_wrapper import dc3dwrapper

x_index = [0, 3, 6]
y_index = [1, 4, 7]
z_index = [2, 5, 8]

def compute_okada(vertices):
    n_pts = vertices.shape[0]
    disp = np.empty((n_pts, 3))
    for i in range(n_pts):
        m = 30e9
        pr = 0.25
        l = (2 * m * pr) / (1 - 2 * pr);
        alpha = (l+m) / (l + 2 * m)
        v = vertices[i,:]
        success, u, grad_u = dc3dwrapper(alpha, v, 2.0, 90, [-1.0, 1.0],
                                         [-1.0, 2.0], [-1.0, 0.0, 0.0])
        if success != 0:
            print "WHOA"
        disp[i, :] = u
    return disp

def main(filename, values_dim):
    f = h5py.File(filename)
    facets = f['locations']
    data = f['values' + str(values_dim)][:,0]
    vertices = np.array([
        facets[:, x_index].flatten(),
        facets[:, y_index].flatten(),
        facets[:, z_index].flatten()
    ]).T

    n_v = vertices.shape[0]
    faces = np.arange(n_v).reshape((n_v / 3, 3))

    # mlab.triangular_mesh(vertices[:, 0], vertices[:, 1], vertices[:, 2],
    #                    faces[:,:], scalars = data)
    # mlab.show()
    okada_exact = compute_okada(vertices)

    vmax = 0.2

    opts = dict(shading = 'gouraud', vmin = -vmax, vmax = vmax)
    plt.figure()
    plt.tripcolor(vertices[:,0], vertices[:,1], faces, data, **opts)
    plt.colorbar()
    plt.title('Output')
    plt.figure()
    plt.tripcolor(vertices[:,0], vertices[:,1], faces, okada_exact[:,values_dim],
                          **opts)
    plt.title('Exact')
    plt.colorbar()
    plt.show()


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
