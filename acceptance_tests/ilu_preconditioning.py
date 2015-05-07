from tbempy.TwoD import dense_integral_operator
from laplace import solve_iterative, np_matrix_from_tbem_matrix
from laplace2d import run, log_u, make_log_dudn
import scipy.sparse.linalg as sp_la
import numpy as np

def test_ILU():
    class Solver(object):
        def __init__(self):
            self.iterations = 0

        def solve(self, tbem, constraint_matrix, matrix, rhs):
            #TODO: Use the sparse nearfield matrix instead
            #TODO: Do my own preconditioning so that I can apply the matrix
            # without condensing it at all... Or should I condense first?
            matrix_condensed = tbem.condense_matrix(
                constraint_matrix, constraint_matrix, matrix.get_block(0,0)
            );
            np_matrix = np_matrix_from_tbem_matrix(matrix_condensed)
            P = sp_la.spilu(np_matrix, drop_tol = 1e-5)
            M_x = lambda x: P.solve(x)
            def mv(v):
                distributed = tbem.distribute_vector(
                    constraint_matrix, v, matrix.n_total_rows()
                )
                applied = matrix.apply(distributed)
                condensed = tbem.condense_vector(constraint_matrix, applied)
                self.iterations += 1
                return condensed
            M = sp_la.LinearOperator((rhs.shape[0], rhs.shape[0]), M_x)
            A = sp_la.LinearOperator((rhs.shape[0], rhs.shape[0]),
                                     matvec = mv, dtype = np.float64)
            res = sp_la.gmres(A, rhs, tol = 1e-8, M = M)
            assert(res[1] == 0) #Check that the iterative solver succeeded
            return res[0]

    solver = Solver()

    bdry_error, int_error = run(
        solver.solve, dense_integral_operator, 7, log_u, make_log_dudn
    )
    assert(solver.iterations < 20)
    assert(bdry_error < 1e-5)
    assert(int_error < 2e-5)
