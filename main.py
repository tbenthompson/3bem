import pyopencl as cl
import pyopencl.clrandom
import pyopencl.array
import pyopencl.characterize
import time, sys, warnings
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.sparse as spsp

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

    if len(x_hat.shape) > 1:
        res[:, 1:2] = x_hat
    else:
        res[:, 1] = x_hat
    if n_max == 1:
        return res
    for i in range(2, n_max + 1):
        res[:, i] = 2 * x_hat * res[:, i - 1] - res[:, i - 2]
    return res

def s_n(x, y, n, a = -1, b = 1):
    x_cheb = cheb_polys(x, n - 1, a, b)
    y_cheb = cheb_polys(y, n - 1, a, b)
    # chop the first value because it is always one
    x_cheb = x_cheb[:, 1:]
    y_cheb = y_cheb[:, 1:]
    return (1.0 / n) + (2.0 / n) * np.sum(x_cheb * y_cheb, axis = 1)

def cheb_pts(n, a = -1, b = 1):
    m = np.arange(1, n + 1)
    return a + (b - a) * 0.5 * (np.cos((2 * m - 1) * np.pi / (2 * n)) + 1)

def tests():
    theta = np.linspace(0.0, 1.0, 4)
    cos_theta = np.cos(theta)
    cheb_test = cheb_polys(cos_theta, 10)
    np.testing.assert_almost_equal(cheb_polys(0.5, 0), [[1]])
    exact = np.cos(np.outer(np.arange(11), theta).T)
    np.testing.assert_almost_equal(exact, cheb_test)

    fnc = lambda x: np.sin(16 * x)
    n_pts_interp = 32

    n_locs = 300
    a = -1.0
    b = 1.0
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
# tests()

def sep_idxs_to_full(sep, depth):
    leaves_1d = 2 ** depth
    leaf_index = sep[0, :] + sep[1, :] * leaves_1d + sep[2, :] * (leaves_1d ** 2)
    return leaf_index

def get_index(x, center, half_width, depth):
    leaves_1d = 2 ** depth
    index = np.floor(((x - center) / (2 * half_width) + 0.5) * leaves_1d).astype(np.uint32)
    leaf_index = sep_idxs_to_full(index, depth)
    return leaf_index

def direct(n, x, v, kernel):
    start = time.time()
    start = time.time()
    pos = np.swapaxes(np.tile(x, (x.shape[1], 1, 1)), 0, 1)
    posT = np.swapaxes(pos, 1, 2)
    K = kernel(pos, posT)
    np.fill_diagonal(K, 0)
    vec = K.dot(v)
    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")
    #     for i in range(n):
    #         K = kernel(x[:, i:(i + 1)], x)
    #         K[K == np.inf] = 0
    #         vec[i] = np.sum(v * K)
    print("Direct sum: " + str(time.time() - start))
    return vec

def level_nodes(level, order, min_corner, max_corner):
    n_cells_1d = 2 ** level
    n_cells = 8 ** level
    n_nodes_1d = order + 1
    n_nodes_per_cell = (order + 1) ** 3
    width = (max_corner - min_corner).reshape(-1, 1)
    cell_width = width / n_cells_1d
    min_ref_1d = [np.linspace(0.0, 1.0, n_cells_1d + 1)[:-1] for d in [0,1,2]]
    min_ref = [min_dim.reshape(n_cells) for min_dim in np.meshgrid(*min_ref_1d)]
    min_ref[0], min_ref[1] = min_ref[1], min_ref[0]
    min_ref[0], min_ref[2] = min_ref[2], min_ref[0]
    mins = min_corner.reshape(-1, 1) + width * min_ref
    maxs = mins + cell_width
    C = cheb_pts(order + 1, a = 0.0, b = 1.0)
    node_ref = np.array([Cdim.reshape(n_nodes_per_cell) for Cdim in np.meshgrid(C, C, C)])
    node_ref_all = np.tile(node_ref, (1, n_cells))
    node_min = np.repeat(mins, n_nodes_per_cell, axis = 1)
    node_max = np.repeat(maxs, n_nodes_per_cell, axis = 1)
    nodes = node_min + (node_max - node_min) * node_ref_all
    return mins, maxs, nodes

def particle_to_multipole(x, values, index, min_corner, max_corner, nodes, order):
    # print x.shape, values.shape, index.shape, min_corner.shape, max_corner.shape, nodes.shape, order
    op_shape = (nodes.shape[1], values.shape[0])


    n_nodes_1d = order + 1
    n_nodes_per_cell = n_nodes_1d ** 3
    start_node_idx = index * n_nodes_per_cell
    node_ref_range = np.arange(n_nodes_per_cell, dtype = np.uint32).reshape(-1, 1)
    row_2d = start_node_idx + node_ref_range
    row = row_2d.T.reshape(-1, 1)[:, 0]
    column = np.repeat(np.arange(x.shape[1], dtype = np.uint32), n_nodes_per_cell)
    cell = np.repeat(index, n_nodes_per_cell)
    data = []
    for d in [0, 1, 2]:
        data.append(s_n(nodes[d, row], x[d, column], order + 1,
                        a = min_corner[d, cell], b = max_corner[d, cell]))
    data = np.prod(data, axis = 0)
    matrix = spsp.csr_matrix((data, (row, column)), shape = op_shape)
    # plt.spy(matrix)
    # plt.show()
    result = matrix.dot(values)
    return result

def multipole_to_particle(x, weights, min_corner, max_corner, nodes, order, kernel):
    n = x.shape[1]
    vec = np.empty(n)
    for i in range(n):
        K = kernel(x[:, i:(i + 1)], nodes)
        vec[i] = np.sum(K * weights)
    return vec

def multipole_to_particle2(x, values, index, weights, centers, widths, nodes, children, order, kernel, leaf_buckets):
    n = x.shape[1]
    vec = np.zeros(n)
    depth = len(weights) - 1
    n_nodes_1d = order + 1
    n_nodes_per_cell = n_nodes_1d ** 3
    all_weights = np.hstack(weights)
    all_nodes = np.hstack(nodes)
    chunk_start = [0]
    for w_chunk in weights:
        chunk_start.append(len(w_chunk) + chunk_start[-1])

    for i in range(n):
        interact_with_level = [0]
        interact_with_cell = [0]
        pt = x[:, i]
        neighbors = []
        m2p = []
        while interact_with_level:
            level = interact_with_level.pop()
            cell = interact_with_cell.pop()
            center_a = centers[level][:, cell]
            width_a = widths[level][:, cell]
            dist_ratio = (pt - center_a) ** 2 / width_a ** 2
            dist_ratio = sum(dist_ratio)
            if dist_ratio < 3:
                if level == depth:
                    neighbors.extend(leaf_buckets[cell])
                else:
                    interact_with_level.extend([level + 1] * 8)
                    interact_with_cell.extend(children[level][(8 * cell):(8 * (cell + 1))])
            else:
                # M2P
                start_node = chunk_start[level] + cell * n_nodes_per_cell
                end_node = start_node + n_nodes_per_cell
                m2p.extend(range(start_node, end_node))

        if neighbors:
            res = kernel(x[:, i:(i + 1)], x[:, neighbors])
            res[res == np.inf] = 0
            vec[i] += np.sum(res * values[neighbors])

        if m2p:
            K = kernel(x[:, i:(i + 1)], all_nodes[:, m2p])
            vec[i] += np.sum(K * all_weights[m2p])
    return vec



def main():
    # kernel = lambda x, y: np.ones(y.shape[1])
    kernel = lambda x, y: 1.0 / np.sqrt(np.sum((x - y) ** 2, axis = 0))
    # kernel = lambda x, y: np.sqrt(np.sum((x - y) ** 2, axis = 0))


    n = np.int(1e4)
    # xx, yy, zz = np.mgrid[0:1:10j, 0:1:20j, 0:1:5j]
    # locs = np.zeros((3, n))
    # locs[0, :] = xx.reshape((n))
    # locs[1, :] = yy.reshape((n))
    # locs[2, :] = zz.reshape((n))
    locs = np.random.rand(3, n)
    locs[1, :] -= 1
    locs[2, :] += 1
    # values = np.ones(n)
    values = np.random.rand(n)

    depth = 2
    order = 2

    start = time.time()
    # find the extents of the octree
    min_x = np.min(locs, axis = 1)
    max_x = np.max(locs, axis = 1)
    octree_center = ((max_x + min_x) / 2).reshape(-1, 1)
    half_width = ((max_x - min_x) / 2).reshape(-1, 1)

    # ensure that all points in the *interior* and
    # none are on the boundary of the octree
    width_fudge = 1e-2
    half_width *= 1 + width_fudge

    min_octree = octree_center - half_width
    max_octree = octree_center + half_width
    all_index = [get_index(locs, octree_center, half_width, cur_depth) for cur_depth in range(0, depth + 1)]
    leaf_index = all_index[-1]

    sorted_order = np.argsort(leaf_index)
    leaf_index = leaf_index[sorted_order]
    locs = locs[:, sorted_order]
    values = values[sorted_order]
    for i in range(0, depth + 1):
        all_index[i] = all_index[i][sorted_order]

    leaf_buckets = [list([]) for _ in xrange(8 ** depth)]
    for i, li in enumerate(leaf_index):
        leaf_buckets[li].append(i)

    maxs = []
    mins = []
    nodes = []
    centers = []
    widths = []
    index_up = []
    index_down = []
    child_nodes = []
    eps = 1e-13
    for cur_depth in range(depth + 1):
        min_locs, max_locs, nodes_locs = level_nodes(cur_depth, order, min_octree, max_octree)
        mins.append(min_locs)
        maxs.append(max_locs)
        nodes.append(nodes_locs)
        centers.append((min_locs + max_locs) / 2.0)
        widths.append((max_locs - min_locs) / 2.0)
        index_up.append(get_index(nodes[-1], octree_center, half_width, cur_depth - 1))
        offsets = [g.reshape(8) for g in np.mgrid[-eps:eps:2j, -eps:eps:2j, -eps:eps:2j]]
        rep_offsets = np.tile(np.array(offsets), (1, centers[-1].shape[1]))
        rep_centers = np.repeat(centers[-1], 8, axis = 1)
        child_nodes.append(get_index(rep_centers + rep_offsets, octree_center, half_width, cur_depth + 1))
    print "Treecode1: " + str(time.time() - start)

    # print child_nodes

    # for i in range(n):
    #     for cur_depth in range(depth + 1):
    #         n_nodes_per_cell = (order + 1) ** 3
    #         my_i = all_index[cur_depth][i]
    #         up_i = index_up[cur_depth][n_nodes_per_cell * my_i]
    #         up_depth = cur_depth - 1
    #         for d in range(3):
    #             min_my = mins[cur_depth][d, my_i]
    #             max_my = maxs[cur_depth][d, my_i]
    #             # print min_my, locs[d, i], max_my
    #             assert(min_my < locs[d, i] < max_my)
    #             if cur_depth >= 1:
    #                 min_up = mins[up_depth][d, up_i]
    #                 max_up = maxs[up_depth][d, up_i]
    #                 # print min_up, locs[d, i], max_up
    #                 assert(min_up < locs[d, i] < max_up)
    # print "Treecode2: " + str(time.time() - start)


    weights = [0] * (depth + 1)
    weights2 = [0] * (depth + 1)
    weights[-1] = particle_to_multipole(locs, values, leaf_index, mins[-1],
                                               maxs[-1], nodes[-1], order)
    print "Treecode3: " + str(time.time() - start)

    for up_pass in range(1, depth + 1):
        cur_level = depth - up_pass
        down_level = depth - up_pass + 1

        weights[cur_level] = particle_to_multipole(nodes[down_level], weights[down_level], index_up[down_level], mins[cur_level], maxs[cur_level], nodes[cur_level], order)
        # weights2[cur_level] = particle_to_multipole(locs, values, all_index[cur_level], mins[cur_level], maxs[cur_level], nodes[cur_level], order)
        # error = np.mean(np.abs(weights[cur_level] - weights2[cur_level])) / np.mean(np.abs(weights2[cur_level]))
        # assert(error < 1e-14)

    # from_level = depth - 1
    # vec_treecode = multipole_to_particle(locs, weights[from_level],
    #                                      mins[from_level],
    #                                      maxs[from_level],
    #                                      nodes[from_level],
    #                                      order, kernel)

    # print "Treecode4: " + str(time.time() - start)
    # vec_treecode = multipole_to_particle2(locs, values, all_index, weights, centers, widths, nodes, child_nodes, order, kernel, leaf_buckets)

    print "Treecode5: " + str(time.time() - start)
    vec_exact = direct(n, locs, values, kernel)
    # np.testing.assert_almost_equal(vec_treecode, vec_exact, 2)

if __name__ == "__main__":
    main()
