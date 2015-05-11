import tbempy.TwoD
from dislocation import solve
from planestrain_fault import exact_displacements, check_planestrain_error
import scipy.sparse.linalg
import scipy.sparse
import numpy as np

# Note: need to scale rows so that the residual is properly evaluated. Simply
# dividing by the shear modulus should work.
shear_modulus = 30e9

def form_preconditioner(tbem, constraint_matrix, op):
    nearfield_condensed = tbem.condense_matrix(
        constraint_matrix, constraint_matrix, op.get_nearfield_matrix()
    )

    scaled_values = nearfield_condensed.values / shear_modulus

    np_matrix = scipy.sparse.csr_matrix((
        scaled_values,
        nearfield_condensed.column_indices,
        nearfield_condensed.row_ptrs
    ))

    np_matrix = np_matrix.tocsc()

    P = scipy.sparse.linalg.spilu(np_matrix)
    def precondition(x):
        # return x
        # print("PREC")
        return P.solve(x)
    return precondition

class Solver(object):
    def __init__(self):
        self.iterations = 0

    def solve(self, tbem, constraint_matrix, op, rhs):
        rhs /= shear_modulus
        def mv(v):
            distributed = tbem.distribute_vector(
                constraint_matrix, v, op.n_total_rows()
            )
            applied = op.apply(distributed)
            applied /= shear_modulus
            condensed = tbem.condense_vector(constraint_matrix, applied)
            self.iterations += 1
            return condensed
        n = rhs.shape[0]

        prec_fnc = form_preconditioner(tbem, constraint_matrix, op)
        M = scipy.sparse.linalg.LinearOperator((n, n), matvec = prec_fnc)
        A = scipy.sparse.linalg.LinearOperator((n, n), matvec = mv)

        res = scipy.sparse.linalg.lgmres(A, rhs, tol = 1e-8, M = M)
        assert(res[1] == 0) #Check that the iterative solver succeeded
        return res[0]

def test_ILU():
    fault = tbempy.TwoD.line_mesh([-1, -1], [0, 0])
    surface = tbempy.TwoD.line_mesh([-100, 0], [100, 0]).refine_repeatedly(9)
    slip = np.ones(2 * fault.n_dofs())
    qs = tbempy.TwoD.QuadStrategy(3, 8, 5.0, 1e-4)
    hyp = tbempy.TwoD.ElasticHypersingular(shear_modulus, 0.25)

    solver = Solver()
    soln = solve(2, surface, fault, hyp, qs, slip, linear_solver = solver.solve)
    check_planestrain_error(surface, soln)
    assert(solver.iterations < 30)

if __name__ == "__main__":
    test_ILU()
