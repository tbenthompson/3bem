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

def faulted_surface_constraints(tbem, surface, fault, disp_dim):
    continuity = tbem.mesh_continuity(surface.begin())
    cut_continuity = tbem.cut_at_intersection(
        continuity, surface.begin(), fault.begin()
    )
    one_component = tbem.convert_to_constraints(cut_continuity)
    all_components = []
    for d in range(disp_dim):
        all_components += tbem.shift_constraints(one_component, d * surface.n_dofs())
    constraint_matrix = tbem.from_constraints(all_components)
    return constraint_matrix

def default_linear_solver(tbem, cm, lhs, rhs):
    def mv(v):
        mv.it+=1
        soln = tbem.distribute_vector(cm, v, lhs.n_rows())
        applied = lhs.apply(soln)
        out = tbem.condense_vector(cm, applied)
        print(mv.it)
        return out
    mv.it = 0

    A = sp_la.LinearOperator((rhs.shape[0], rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, rhs, tol = 1e-7)
    return res[0]

def solve(dim, surface, fault, mthd, slip, **kwargs):
    tbem = get_tbem(dim)

    linear_solver = kwargs.get('linear_solver', default_linear_solver)

    all_mesh = tbem.Mesh.create_union([surface, fault])

    n_dofs = tbem.dim * surface.n_dofs()

    cm = faulted_surface_constraints(tbem, surface, fault, dim)

    rhs_op = tbem.integral_operator(surface, fault, mthd, all_mesh)
    all_dofs_rhs = rhs_op.apply(slip)
    rhs = -tbem.condense_vector(cm, all_dofs_rhs)

    lhs = tbem.integral_operator(surface, surface, mthd, all_mesh)
    soln = linear_solver(tbem, cm, lhs, rhs)
    full_soln = tbem.distribute_vector(cm, soln, n_dofs)
    soln_components = []
    for d in range(tbem.dim):
        start_idx = surface.n_dofs() * d
        end_idx = surface.n_dofs() * (d + 1)
        soln_components.append(full_soln[start_idx:end_idx])
    return soln_components
