from tbempy import *
import tbempy.TwoD
import tbempy.ThreeD
import scipy.sparse.linalg as sp_la
import numpy as np

import time
class Timer:
    def __init__(self):
        self.start = time.time()
    def print_reset(self, comment):
        print(comment + ": " + str(time.time() - self.start))
        self.reset()
    def reset(self):
        self.start = time.time()

def get_tbem(dim):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    return tbem

def faulted_surface_constraints(dim, surface, fault):
    tbem = get_tbem(dim)
    continuity = tbem.mesh_continuity(surface.begin())
    cut_continuity = tbem.cut_at_intersection(
        continuity, surface.begin(), fault.begin()
    )
    constraints = tbem.convert_to_constraints(cut_continuity)
    constraint_matrix = from_constraints(constraints)
    return constraint_matrix

def solve(dim, surface, fault, hyp, qs, slip):
    linear_solve_tol = 1e-8
    n_dofs = surface.n_dofs()
    print("Number of dofs: " + str(dim * n_dofs))
    tbem = get_tbem(dim)
    t = Timer()

    mthd = tbem.make_adaptive_integration_mthd(qs, hyp)
    cm = faulted_surface_constraints(dim, surface, fault)
    t.print_reset("build constraints")
    rhs_op = tbem.integral_operator(surface, fault, mthd)
    t.print_reset("build rhs op")
    all_dofs_rhs = rhs_op.apply(slip)
    t.print_reset("apply rhs")
    rhs = BlockVectorX([
        condense_vector(cm, all_dofs_rhs.storage[i]) for i in range(dim)
    ])

    dof_map = block_dof_map_from_functions(rhs)
    rhs = concatenate(dof_map, rhs)

    lhs = tbem.integral_operator(surface, surface, mthd)
    t.print_reset("Construct lhs")

    def mv(v):
        t.reset()
        vec_v = VectorX(v)
        t.print_reset('convert to VectorX')
        both = expand(dof_map, vec_v)
        t.print_reset('expand')

        test = BlockVectorX([
            distribute_vector(cm, both.storage[i], n_dofs)
            for i in range(dim)
        ])
        t.print_reset('distribute')
        applied = lhs.apply(test)
        t.print_reset('apply')
        condensed = BlockVectorX([
            condense_vector(cm, applied.storage[i])
            for i in range(dim)
        ])
        t.print_reset('condense')
        out = concatenate(dof_map, condensed)
        t.print_reset('concatenate')
        print("IT:" + str(mv.it))
        mv.it+=1
        return np.array(out.storage)
    mv.it = 0

    np_rhs = np.array(rhs.storage)
    A = sp_la.LinearOperator((np_rhs.shape[0], np_rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, np_rhs, tol = linear_solve_tol)

    res_vec = VectorX(res[0])
    soln_expanded = expand(dof_map, res_vec)

    soln = [
        np.array(distribute_vector(cm, soln_expanded.storage[i], n_dofs).storage)
        for i in range(dim)
    ]
    return soln
