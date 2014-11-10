import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
import mayavi.mlab as mlab
from okada_wrapper import dc3dwrapper

def main(filename, values_dim):
    f = h5py.File(filename)
    faces = f['faces']
    vertices = f['vertices']
    data = f['values'][:,values_dim]

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
    plt.tricontourf(vertices[:,0], vertices[:,1], data, shading = 'gouraud')
    plt.tricontour(vertices[:,0], vertices[:,1], data,
                   colors = ['k'], linestyles = 'solid', shading = 'gouraud')
    plt.figure()
    plt.tricontourf(vertices[:,0], vertices[:,1], disp[:,0], shading = 'gouraud')
    plt.tricontour(vertices[:,0], vertices[:,1], disp[:,0],
                   colors = ['k'], linestyles = 'solid', shading = 'gouraud')
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
