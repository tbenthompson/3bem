import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
from okada_wrapper import dc3dwrapper

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

    # mlab.triangular_mesh(vertices[:, 0], vertices[:, 1], vertices[:, 2],
    #                    faces[:,:], scalars = data)
    # mlab.show()
    n_pts = vertices.shape[0]
    disp = np.empty((n_pts, 3))
    for i in range(n_pts):
        m = 30e9
        pr = 0.25
        l = (2 * m * pr) / (1 - 2 * pr);
        alpha = (l+m) / (l + 2 * m)
        v = vertices[i,:]
        success, u, grad_u = dc3dwrapper(alpha, v, 2.0, 90, [-1.0, 1.0],
                                         [-1.0, 1.0], [-1.0, 0.0, 0.0])
        disp[i,:] = u
    print disp[0,:]
    plt.figure()
    trip1 = plt.tripcolor(vertices[:,0], vertices[:,1], data, shading = 'gouraud', vmin = -0.04, vmax = 0.04)
    # plt.tricontourf(vertices[:,0], vertices[:,1], data, shading = 'gouraud')
    # plt.tricontour(vertices[:,0], vertices[:,1], data,
    #                colors = ['k'], linestyles = 'solid', shading = 'gouraud')
    plt.colorbar()
    plt.figure()
    trip2 = plt.tripcolor(vertices[:,0], vertices[:,1], disp[:,0], shading = 'gouraud', vmin = -0.04, vmax = 0.04)
    # plt.tricontourf(vertices[:,0], vertices[:,1], disp[:,0], shading = 'gouraud')
    # plt.tricontour(vertices[:,0], vertices[:,1], disp[:,0],
    #                colors = ['k'], linestyles = 'solid', shading = 'gouraud')
    plt.colorbar()

    diff = disp[:,0] - data#np.log(np.abs((disp[:,0] - data) / 1)) / np.log(10)
    plt.figure()
    trip3 = plt.tripcolor(vertices[:,0], vertices[:,1], diff, shading = 'gouraud', vmin = -0.04, vmax = 0.04)
    plt.colorbar()
    plt.figure()
    trip3 = plt.tripcolor(vertices[:,0], vertices[:,1], diff, shading = 'gouraud')
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
