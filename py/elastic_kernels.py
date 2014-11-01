import numpy as np
from math import pi
import sympy as sp

# I was unable to find a good source on the exact form of the
# fundamental solutions for the elastic boundary integral equations
# Some are slightly wrong, etc.
#
# So, I symbolically derive an ugly form here and then test it against
# the hand derived version.

x1, x2, y1, y2, z1, z2 = sp.symbols('x1, x2, y1, y2, z1, z2')
sm, pr = sp.symbols('sm, pr')
nx, ny, nz = sp.symbols('nx, ny, nz')
mx, my, mz = sp.symbols('mx, my, mz')
sm_val = 30e9
pr_val = 0.25

def disp_creator(k, j):
    C1 = 1 / (16 * sp.pi * sm * (1 - pr))
    C2 = (3 - 4 * pr)
    delta = [x1 - x2, y1 - y2, z1 - z2]
    r2 = sum([d**2 for d in delta])

    disp = (C1 / sp.sqrt(r2)) *\
           (C2 * (1 if k == j else 0) + delta[k] * delta[j] / r2)
    return disp

def traction_operator(disp_vec, pos_vec, index, normal):
    disp_grad = [[sp.diff(disp_vec[l], pos_vec[m])
                    for m in [0,1,2]] for l in [0,1,2]]
    strain = [[(disp_grad[l][m] + disp_grad[m][l]) / 2
                    for m in [0,1,2]] for l in [0,1,2]]
    strain_trace = (strain[0][0] + strain[1][1] + strain[2][2])

    lame_lambda = (2 * sm * pr) / (1 - 2 * pr)
    def stress_from_strain(l, m):
        return lame_lambda * strain_trace * (1 if l == m else 0) +\
               2 * sm * strain[l][m]
    stress = [[stress_from_strain(l, m) for m in [0,1,2]] for l in [0,1,2]]
    trac = stress[index][0] * normal[0] +\
           stress[index][1] * normal[1] +\
           stress[index][2] * normal[2]
    return trac

def derive_kernels(k, j):
    disp = disp_creator(k, j)

    src_disp_vec = [disp_creator(k, dir) for dir in [0,1,2]]
    src_pos_vec = [x1, y1, z1]
    src_n = [nx, ny, nz]

    obs_disp_vec = [disp_creator(dir, j) for dir in [0,1,2]]
    obs_pos_vec = [x2, y2, z2]
    obs_n = [mx, my, mz]

    trac = traction_operator(src_disp_vec, src_pos_vec, j, src_n)

    adj_trac = traction_operator(obs_disp_vec, obs_pos_vec, k, obs_n)

    tracx = traction_operator([disp_creator(0, d) for d in [0,1,2]],
                              src_pos_vec, j, src_n)
    tracy = traction_operator([disp_creator(1, d) for d in [0,1,2]],
                              src_pos_vec, j, src_n)
    tracz = traction_operator([disp_creator(2, d) for d in [0,1,2]],
                              src_pos_vec, j, src_n)
    hyp_disp_vec = [tracx, tracy, tracz]
    hyp = traction_operator(hyp_disp_vec, obs_pos_vec, k, obs_n)

    return disp, trac, adj_trac, hyp

def get_arg_sets():

    n_arg_sets = 100

    arg_sets = []
    for i in range(n_arg_sets):
        vec1 = np.random.rand(3)
        vec2 = np.random.rand(3)
        n_vec = np.random.rand(3)
        n_vec /= np.linalg.norm(n_vec)
        m_vec = np.random.rand(3)
        m_vec /= np.linalg.norm(m_vec)
        arg_sets.append([sm_val, pr_val] + vec1.tolist() + vec2.tolist() +
                        n_vec.tolist() + m_vec.tolist())
    return arg_sets

def main(pair, data_file):
    k, j = pair
    arg_sets = get_arg_sets()

    disp, trac, adj_trac, hyp = derive_kernels(k, j)
    print k,j

    args = (sm, pr, x1, y1, z1, x2, y2, z2, nx, ny, nz, mx, my, mz)
    disp_fnc = sp.utilities.lambdify(args, disp)
    trac_fnc = sp.utilities.lambdify(args, trac)
    adj_trac_fnc = sp.utilities.lambdify(args, adj_trac)
    hyp_fnc = sp.utilities.lambdify(args, hyp)

    for fnc, name in [(disp_fnc, "disp"), (trac_fnc, "trac"),
                      (adj_trac_fnc, "adj_trac"), (hyp_fnc, "hyp")]:
        for args in arg_sets:
            value = fnc(*args)
            data_list = [value, name, k, j, sm_val, pr_val] + args[2:]
            str_data = [str(thing) for thing in data_list]
            data = ','.join(str_data)
#output value, "trac", i, j, sm, pr, x1, x2, y1, y2, z1, z2, nx, ny, nz, mx, my, mz
            data_file.write(data + '\n')


if __name__ == "__main__":
    test_data = open('test/test_data', 'w')
    for k in range(3):
        for j in range(3):
            main((k,j), test_data)
    test_data.close()
