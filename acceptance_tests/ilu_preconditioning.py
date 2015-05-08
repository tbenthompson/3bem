from tbempy.TwoD import integral_operator
from laplace import solve_iterative, np_matrix_from_tbem_matrix
from laplace2d import run, log_u, make_log_dudn
import scipy.sparse.linalg
import scipy.sparse
import numpy as np

def test_ILU():
    class Solver(object):
        def __init__(self):
            self.iterations = 0

        def solve(self, tbem, constraint_matrix, op, rhs):
            nearfield = op.get_nearfield_matrix()
            np_matrix = scipy.sparse.csr_matrix((
                nearfield.values,
                nearfield.column_indices,
                nearfield.row_ptrs
            )).tocsc()
            P = scipy.sparse.linalg.spilu(np_matrix, drop_tol = 1e-5)
            # M_x = lambda x: P.solve(x)
            def mv(v):
                distributed = tbem.distribute_vector(
                    constraint_matrix, v, op.n_total_rows()
                )
                applied = op.apply(distributed)
                condensed = tbem.condense_vector(constraint_matrix, applied)
                self.iterations += 1
                return condensed
            # M = scipy.sparse.linalg.LinearOperator((rhs.shape[0], rhs.shape[0]), M_x)
            A = scipy.sparse.linalg.LinearOperator((rhs.shape[0], rhs.shape[0]),
                                     matvec = mv, dtype = np.float64)
            res = scipy.sparse.linalg.gmres(A, rhs, tol = 1e-8)
            assert(res[1] == 0) #Check that the iterative solver succeeded
            return res[0]

    solver = Solver()

    bdry_error, int_error = run(
        solver.solve, integral_operator, 7, log_u, make_log_dudn
    )
    assert(solver.iterations < 20)
    assert(bdry_error < 1e-5)
    assert(int_error < 2e-5)
