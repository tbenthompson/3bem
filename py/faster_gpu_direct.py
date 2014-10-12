import time
import pyopencl as cl
import numpy as np
import numpy.linalg as la

def calc_speed(total_time, n):
    interactions = n ** 2
    ops_per_interact = 15
    if dev == 0:
        proc = 1.0e9
    elif dev == 1:
        proc = 3.4e9
    cycles = total_time * proc
    instr_per_cycle = (interactions * ops_per_interact) / cycles;
    gigaflops = instr_per_cycle * proc / 1e9
    print("# of instructions: " + str(interactions))
    print("Instructions/cycle: " + str(instr_per_cycle))
    print("GFlop/s: " + str(gigaflops))
    print("Time: " + str(total_time))


items_per_tile = 256 * 4
tiles_per_row = 128
n = items_per_tile * tiles_per_row
check = n < 10001
src = np.random.rand(n, 4).astype(np.float32)
obs = np.random.rand(n, 4).astype(np.float32)
strength = np.random.rand(n).astype(np.float32)

dev = 0
ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

mf = cl.mem_flags
src_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=src)
obs_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=obs)
intermediate_buf = cl.Buffer(ctx, mf.READ_WRITE, strength.nbytes * tiles_per_row)
lmem = cl.LocalMemory(items_per_tile * 4 * np.dtype('float32').itemsize)
dest_buf = cl.Buffer(ctx, mf.WRITE_ONLY, strength.nbytes)

# Going from 3 float arrays to one float4 array boosts from 800 to 1200 gflops
# Tiling the computation properly boosts from 1200 -> 2600 on "loose" with one GTX 780 Ti
prg = cl.Program(ctx, """
constant float factor = 1.0 / (4.0 * M_PI);
__kernel void direct_n_body1(global const float4 *src,
                            global const float4 *obs,
                            const unsigned n_src,
                            const unsigned n_obs,
                            global float *result)
{
    int gid = get_global_id(0);
    float sum = 0.0f;
    float4 cur_obs = obs[gid];
    for (int i = 0; i < n_src; i++) {
        float4 cur_src = src[i];
        float dx = cur_obs.x - cur_src.x;
        float dy = cur_obs.y - cur_src.y;
        float dz = cur_obs.z - cur_src.z;
        float r2 = dx * dx;
        r2 += dy * dy;
        r2 += dz * dz;
        float kernel_val = factor * native_rsqrt(r2);
        sum += cur_src.w * kernel_val;
    }
    result[gid] = sum;
}

inline float interact(float4 src, float4 obs) {
    float dx = src.x - obs.x;
    float dy = src.y - obs.y;
    float dz = src.z - obs.z;
    float r2 = dx * dx;
    r2 += dy * dy;
    r2 += dz * dz;
    return factor * rsqrt(r2) * src.w;
}

#define ROW_DIM 0
#define COL_DIM 1
// Tiled direct n body calculation.
// src.w is the strength of the source
__kernel void direct_n_body2(global const float4 *src,
                             global const float4 *obs,
                             const unsigned n_src,
                             const unsigned n_obs,
                             global float *intermediate_result,
                             local float4 *workspace)
{
    // Determine which tile we are working on.
    int gobs_id = get_global_id(ROW_DIM);
    int gtile_id = get_global_id(COL_DIM);
    int lsrc_id = get_local_id(ROW_DIM);
    int nsrcs_per_tile = get_local_size(ROW_DIM);
    int tile_first_src = nsrcs_per_tile * gtile_id;
    int idx = tile_first_src + lsrc_id;

    // add if statement if the number of sources is not a multiple of the tile size.
    workspace[lsrc_id] = src[idx];

    // barrier to ensure local memory is all loaded
    barrier(CLK_LOCAL_MEM_FENCE);

    // Compute effect of this block on the observation point.
    float sum = 0.0;
    float4 cur_obs = obs[gobs_id];
    for (int k = 0; k < nsrcs_per_tile; k++)
    {
        sum += interact(workspace[k], cur_obs);
    }

    // Store in Y (P columns per row)
    intermediate_result[gobs_id + n_obs * gtile_id] = sum;
}


// Reduce M = get_global_size(0) rows of P values in matrix Y.
// Stores the result in first column of Y.
__kernel void reduce_rows(__global float* intermediate_result,
                          __global float* result,
                          int n_obs,
                          int n_tiles_per_row)
{
    int gobs_id = get_global_id(0);
    float sum = 0.0f;
    for (int tile = 0; tile < n_tiles_per_row; tile++) {
        sum += intermediate_result[gobs_id + n_obs * tile];
    }
    result[gobs_id] = sum;
}
""").build()

print "START GPU"
start = time.time()

# args = [src_buf, obs_buf, np.uint32(n), np.uint32(n), dest_buf]
# global_size = (n,)
# local_size = None#(items_per_tile,)
# print global_size, local_size
# prg.direct_n_body1(queue, global_size, local_size, *args)

direct_args = [src_buf, obs_buf, np.uint32(n), np.uint32(n), intermediate_buf, lmem]
global_size = (n, tiles_per_row)
local_size = (items_per_tile, 1)
print global_size, local_size
prg.direct_n_body2(queue, global_size, local_size, *direct_args)

reduce_size = (n,)
reduce_args = [intermediate_buf, dest_buf, np.uint32(n), np.uint32(tiles_per_row)]
prg.reduce_rows(queue, reduce_size, None, *reduce_args)

result = np.empty(n).astype(np.float32)
cl.enqueue_copy(queue, result, dest_buf)

end = time.time()
total_time = end - start
calc_speed(total_time, n)


if check:
    start = time.time()
    print "START Python"
    exact = np.zeros_like(strength)
    for i in range(strength.shape[0]):
        r2 = (obs[i][0] - src[:, 0]) ** 2 + (obs[i][1] - src[:, 1]) ** 2 + (obs[i][2] - src[:, 2]) ** 2
        kernel = 1.0 / (4 * np.pi * np.sqrt(r2))
        exact[i] = np.sum(kernel * strength)
    end = time.time()
    print "TIME Python: " + str(end - start)
    np.testing.assert_almost_equal(exact / n, result / n, 4)

# print(la.norm(a_plus_b - (a+b)), la.norm(a_plus_b))
