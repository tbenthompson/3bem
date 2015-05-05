import tbempy.TwoD
import tbempy.ThreeD
import scipy.sparse.linalg as sp_la
import numpy as np

import time

def get_tbem(dim):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    return tbem

def faulted_surface_constraints(tbem, surface, fault):
    continuity = tbem.mesh_continuity(surface.begin())
    cut_continuity = tbem.cut_at_intersection(
        continuity, surface.begin(), fault.begin()
    )
    constraints = tbem.convert_to_constraints(cut_continuity)
    constraint_matrix = tbem.from_constraints(constraints)
    return constraint_matrix

def condense(tbem, input, cm, n_dofs):
    out_components = []
    for d in range(tbem.dim):
        start_dof = d * n_dofs
        end_dof = (d + 1) * n_dofs
        component = input[start_dof:end_dof]
        out_components.append(tbem.condense_vector(cm, component))
    return out_components

def distribute(tbem, input, cm, n_dofs, n_reduced_dofs):
    out_components = []
    for d in range(tbem.dim):
        start_dof = d * n_reduced_dofs
        end_dof = (d + 1) * n_reduced_dofs
        component = input[start_dof:end_dof]
        out_components.append(tbem.distribute_vector(cm, component, n_dofs))
    return out_components

def solve(dim, surface, fault, hyp, qs, slip):
    linear_solve_tol = 1e-8
    n_dofs = surface.n_dofs()
    print("Number of dofs: " + str(dim * n_dofs))
    tbem = get_tbem(dim)

    mthd = tbem.make_adaptive_integration_mthd(qs, hyp)
    cm = faulted_surface_constraints(tbem, surface, fault)

    rhs_op = tbem.integral_operator(surface, fault, mthd)
    all_dofs_rhs = rhs_op.apply(slip)
    rhs = np.concatenate(condense(tbem, all_dofs_rhs, cm, n_dofs))
    n_reduced_dofs = rhs.shape[0] / dim

    lhs = tbem.integral_operator(surface, surface, mthd)

    def mv(v):
        print("IT: " + str(mv.it))
        mv.it+=1
        soln = np.concatenate(distribute(tbem, v, cm, n_dofs, n_reduced_dofs))
        applied = lhs.apply(soln)
        out = np.concatenate(condense(tbem, applied, cm, n_dofs))
        return out
    mv.it = 0

    A = sp_la.LinearOperator((rhs.shape[0], rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, rhs, tol = linear_solve_tol)
    soln = distribute(tbem, res[0], cm, n_dofs, n_reduced_dofs)
    return soln
