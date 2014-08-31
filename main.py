import pyopencl as cl
import pyopencl.clrandom
import pyopencl.array
import pyopencl.characterize
import time
import numpy as np
ctx = cl.create_some_context(interactive = False)
queue = cl.CommandQueue(ctx)

n = np.int(1e4)
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
    K = kernel(x[i, :], x)
    K[K == np.inf] = 0
    vec[i] = np.sum(values * K)
print("Direct sum: " + str(time.time() - start))

start = time.time()

cell_centers = []
for i in range(0, depth + 1):
    n_cells_1d = 2 ** i
    half_cell_width = half_width / n_cells_1d
    cell_centers.append(np.empty((n_cells_1d, 3)))
    for d in range(3):
        cell_centers[i][:, d] = np.linspace(min_octree[d], max_octree[d],
                                n_cells_1d + 1)[:-1] + half_cell_width[d]
cell_centers = np.array(cell_centers)

leaf_weights = np.zeros((leaves_1d, leaves_1d, leaves_1d, ))
for i in range(index.shape[0]):
    idx = index[i]
    leaf_weights[idx[0], idx[1], idx[2]] += values[i]

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
