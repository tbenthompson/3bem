import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import sys

def main(filename, values_dim):
    f = h5py.File(filename)
    locs = f['locations']
    data = f['values'][:, values_dim]

    x = locs[:,0].reshape(200,200)
    y = locs[:,1].reshape(200,200)
    d = data.reshape(200,200)
    opts = dict(
        levels = np.linspace(-0.5, 0.5, 20)
        )
    plt.contourf(x, y, d, **opts)
    plt.contour(x, y, d,
                      colors = '#333333',
                      linestyles = 'solid',
                      **opts)
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
