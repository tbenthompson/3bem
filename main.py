import pyopencl as cl
import pyopencl.clrandom
import pyopencl.array
import pyopencl.characterize
import time
import numpy as np
import warnings
import matplotlib.pyplot as plt

def cheb_polys(x, n_max, a = -1, b = 1):
    x_hat = 2 * ((x - a) / (b - a)) - 1
    try:
        length = x_hat.shape[0]
    except AttributeError:
        length = 1
    except IndexError:
        length = 1
    res = np.empty((length, n_max + 1))
    res[:, 0] = 1.0
    if n_max == 0:
        return res
    res[:, 1] = x_hat
    if n_max == 1:
        return res
    for i in range(2, n_max + 1):
        res[:, i] = 2 * x_hat * res[:, i - 1] - res[:, i - 2]
    return res

def s_n(x, y, n, a = -1, b = 1):
    x_cheb = cheb_polys(x, n - 1, a, b)[0, 1:]
    y_cheb = cheb_polys(y, n - 1, a, b)[0, 1:]
    return (1.0 / n) + (2.0 / n) * np.sum(x_cheb * y_cheb)

def cheb_pts(n, a = -1, b = 1):
    """Chebyshev points (of the second kind). Includes -1 and 1."""
    m = np.arange(1, n + 1)
    return a + (b - a) * 0.5 * (np.cos((2 * m - 1) * np.pi / (2 * n)) + 1)


def tests():
    theta = np.linspace(0.0, 1.0, 4)
    cos_theta = np.cos(theta)
    cheb_test = cheb_polys(cos_theta, 10)
    np.testing.assert_almost_equal(cheb_polys(0.5, 0), [[1]])
    exact = np.cos(np.outer(np.arange(11), theta).T)
    np.testing.assert_almost_equal(exact, cheb_test)

    fnc = lambda x: np.sin(2 * x)
    n_pts_interp = 10

    n_locs = 100
    a = 0.0
    b = 0.5
    eval_locs = np.linspace(-1, 1, n_locs)
    exact = fnc(eval_locs)

    cheb_nodes = cheb_pts(n_pts_interp, a = a, b = b)

    est = np.zeros_like(exact)
    for i in range(n_locs):
        for m in range(n_pts_interp):
            est[i] += fnc(cheb_nodes[m]) * \
                      s_n(cheb_nodes[m], eval_locs[i], n_pts_interp, a, b)

    plt.plot(eval_locs, exact, 'k-')
    plt.plot(eval_locs, est, 'b-')
    plt.plot(cheb_nodes, fnc(cheb_nodes), 'o')
    plt.show()
    plt.plot(eval_locs, np.log(np.abs(exact - est) + 1e-16) / np.log(10))
    plt.show()
tests()
ctx = cl.create_some_context(interactive = False)
queue = cl.CommandQueue(ctx)

n = np.int(4e3)
depth = 1
leaves_1d = 2 ** depth

x = np.random.rand(n, 3).astype(np.float32)
x[:, 1] -= 1
x[:, 2] += 1
values = np.random.rand(n).astype(np.float32)
# values = np.ones(n).astype(np.float32)
kernel = lambda x, y: 1.0 / np.sqrt(np.sum((x - y) ** 2, axis = 1))
kernel2 = lambda x, y: 1.0 / np.sqrt(np.sum((x - y) ** 2))
# kernel = lambda x, y: np.ones(y.shape[0])
# kernel2 = lambda x, y: 1.0
# x_g = cl.array.to_device(queue, np.empty((n, 3)).astype(np.float32))
# cl.clrandom.fill_rand(x_g)
# x = np.empty((n, 3)).astype(np.float32)
# cl.enqueue_copy(queue, x, x_g.data)


min_x = np.min(x, axis = 0)
max_x = np.max(x, axis = 0)
center = (max_x + min_x) / 2
half_width = (max_x - min_x) / 2
eps = 1e-6
half_width *= 1 + eps

min_octree = center - half_width
max_octree = center + half_width
index = np.floor(((x - center) / (2 * half_width) + 0.5) * leaves_1d).astype(np.uint32)

# start = time.time()
#
# mf = cl.mem_flags
# # x_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=x)
# c_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=center)
# hw_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=half_width)
# index_g = cl.Buffer(ctx, mf.WRITE_ONLY, x.nbytes)
#
# prg = cl.Program(ctx, """
# __kernel void octree_loc(__global const float *x_g,
#                          __global const float *c_g,
#                          __global const float *hw_g,
#                          const uint leaves_1d,
#                          __global int *index_g) {
#     int g_addr = get_global_id(0) * 3;
#     for (int i = 0; i < 3; ++i) {
#         index_g[g_addr + i] = floor(((x_g[g_addr + i] - c_g[i])
#                                     / (2 * hw_g[i]) + 0.5) * leaves_1d);
#     }
# }
# """).build()
#
# prg.octree_loc(queue, (n, 1), None, x_g.data, c_g,
#                hw_g, np.uint32(leaves_1d), index_g)
#
# index = np.empty_like(x).astype(np.uint32)
# cl.enqueue_copy(queue, index, index_g)

# end = time.time()
# print("ELAPSED: " + str(end - start))

start = time.time()
vec = np.empty(n)
for i in range(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        K = kernel(x[i, :], x)
    K[K == np.inf] = 0
    vec[i] = np.sum(values * K)
print("Direct sum: " + str(time.time() - start))

start = time.time()
expansion_order = 4

cell_nodes = []
for i in range(0, depth + 1):
    n_cells_1d = 2 ** i
    cell_width = 2 * half_width / n_cells_1d
    cell_nodes.append(np.empty((n_cells_1d, 3, expansion_order)))
    for d in range(3):
        min_cell = np.linspace(min_octree[d], max_octree[d], n_cells_1d + 1)
        max_cell = min_cell + cell_width[d]
        for c in range(n_cells_1d):
            cell_nodes[i][c, d, :] =\
                    cheb_pts(expansion_order, min_cell[c], max_cell[c])
cell_nodes = np.array(cell_nodes)
import ipdb; ipdb.set_trace()

leaf_weights = np.zeros((leaves_1d, leaves_1d, leaves_1d, expansion_order, expansion_order, expansion_order))
for i in range(index.shape[0]):
    idx = index[i]
    leaf_weights[idx[0], idx[1], idx[2], 0, 0, 0] += values[i]

vec_treecode = np.empty(n)
for i in range(n):
    K = np.empty((leaves_1d, leaves_1d, leaves_1d))
    for l1 in range(leaves_1d):
        for l2 in range(leaves_1d):
            for l3 in range(leaves_1d):
                center = np.array([cell_centers[-1][l1, 0],
                                   cell_centers[-1][l2, 1],
                                   cell_centers[-1][l3, 2]])
                K[l1, l2, l3] = kernel2(x[i, :], center)
    vec_treecode[i] = np.sum(leaf_weights * K)
print("Treecode: " + str(time.time() - start))

# np.testing.assert_almost_equal(vec_treecode, vec)
