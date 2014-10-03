import time
import pyopencl as cl
import numpy as np
import numpy.linalg as la

n = 300000
check = n < 10001
srcx = np.random.rand(n).astype(np.float32)
srcy = np.random.rand(n).astype(np.float32)
srcz = np.random.rand(n).astype(np.float32)
obsx = np.random.rand(n).astype(np.float32)
obsy = np.random.rand(n).astype(np.float32)
obsz = np.random.rand(n).astype(np.float32)
strength = np.random.rand(n).astype(np.float32)

ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

mf = cl.mem_flags
srcx_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=srcx)
srcy_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=srcy)
srcz_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=srcz)
obsx_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=obsx)
obsy_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=obsy)
obsz_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=obsz)
str_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=strength)
dest_buf = cl.Buffer(ctx, mf.WRITE_ONLY, srcx.nbytes)

prg = cl.Program(ctx, """
constant float factor = 1.0 / (4.0 * M_PI);
__kernel void direct_n_body(global const float *srcx,
                            global const float *srcy,
                            global const float *srcz,
                            global const float *obsx,
                            global const float *obsy,
                            global const float *obsz,
                            global const float *str,
                            const unsigned n_src,
                            const unsigned n_obs,
                            global float *result)
{
    int gid = get_global_id(0);
    float sum = 0.0f;
    float x = obsx[gid];
    float y = obsy[gid];
    float z = obsz[gid];
    for (int i = 0; i < n_src; i++) {
        float dx = x - srcx[i];
        float dy = y - srcy[i];
        float dz = z - srcz[i];
        float r2 = dx * dx;
        r2 += dy * dy;
        r2 += dz * dz;
        float kernel_val = factor * native_rsqrt(r2);
        sum += str[i] * kernel_val;
    }
    result[gid] = sum;
}
""").build()

start = time.time()
print "START GPU"
prg.direct_n_body(queue, obsx.shape, None, srcx_buf, srcy_buf, srcz_buf, obsx_buf, obsy_buf, obsz_buf, str_buf, np.uint32(n), np.uint32(n), dest_buf)

result = np.empty_like(srcx)
cl.enqueue_copy(queue, result, dest_buf)
end = time.time()
total_time = end - start
interactions = n ** 2
ops_per_interact = 15
proc = 3.4e9
cycles = total_time * proc
instr_per_cycle = (interactions * ops_per_interact) / cycles;
gigaflops = instr_per_cycle * proc / 1e9
print("Instructions/cycle: " + str(instr_per_cycle))
print("GFlop/s: " + str(gigaflops))
print("Time: " + str(total_time))

if check:
    start = time.time()
    print "START Python"
    exact = np.zeros_like(obsx)
    for i in range(obsx.shape[0]):
        r2 = (obsx[i] - srcx) ** 2 + (obsy[i] - srcy) ** 2 + (obsz[i] - srcz) ** 2
        kernel = 1.0 / (4 * np.pi * np.sqrt(r2))
        exact[i] = np.sum(kernel * strength)
    end = time.time()
    print "TIME Python: " + str(end - start)
    np.testing.assert_almost_equal(exact / n, result / n, 4)

# print(la.norm(a_plus_b - (a+b)), la.norm(a_plus_b))
