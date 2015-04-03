from tbempy import *
import tbempy.TwoD
import tbempy.ThreeD
import scipy.sparse.linalg as sp_la
import numpy as np

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

def solve(dim, surface, fault, hyp):
    linear_solve_tol = 1e-8
    n_dofs = surface.n_dofs()
    tbem = get_tbem(dim)

    qs = tbem.QuadStrategy(2, 2, 5, 3.0, 1e-3)
    cm = faulted_surface_constraints(dim, surface, fault)

    slip = BlockVectorX(dim, VectorX(fault.n_dofs(), 1.0))

    p_rhs = tbem.make_boundary_integral(surface, fault, hyp);
    print("Building RHS operator")
    rhs_op = tbem.mesh_to_mesh_operator(p_rhs, qs);
    print("Done building RHS operator")
    all_dofs_rhs = rhs_op.apply(slip);
    rhs = BlockVectorX([
        condense_vector(cm, all_dofs_rhs.storage[i]) for i in range(dim)
    ])

    dof_map = block_dof_map_from_functions(rhs)
    rhs = concatenate(dof_map, rhs)

    p_lhs = tbem.make_boundary_integral(surface, surface, hyp);
    print("Building LHS operator")
    lhs = tbem.mesh_to_mesh_operator(p_lhs, qs);
    print("Done building LHS operator")

    def mv(v):
        vec_v = VectorX(v)
        both = expand(dof_map, vec_v)
        test = BlockVectorX([
            distribute_vector(cm, both.storage[i], n_dofs)
            for i in range(dim)
        ])
        applied = lhs.apply(test)
        condensed = BlockVectorX([
            condense_vector(cm, applied.storage[i])
            for i in range(dim)
        ])
        out = concatenate(dof_map, condensed)
        print("ITERATION")
        return np.array(out.storage)

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
