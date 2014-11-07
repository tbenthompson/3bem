import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
import mayavi.mlab as mlab

def main(filename, values_dim):
    f = h5py.File(filename)
    faces = f['faces']
    vertices = f['vertices']
    data = f['values'][:,values_dim]
    print len(data)
    print vertices.shape

    plt.tripcolor(vertices[:,0], vertices[:,1], data, shading = 'gouraud')
    plt.show()
    # mlab.triangular_mesh(vertices[:, 0], vertices[:, 1], vertices[:, 2],
    #                    faces[:,:], scalars = data)
    # mlab.show()

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
