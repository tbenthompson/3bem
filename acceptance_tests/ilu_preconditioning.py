from tbempy.TwoD import dense_integral_operator, integral_operator
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
            # nearfield_condensed = tbem.condense_matrix(
            #     constraint_matrix, constraint_matrix, op.get_nearfield_matrix()
            # )

            # np_matrix = scipy.sparse.csr_matrix((
            #     nearfield_condensed.values,
            #     nearfield_condensed.column_indices,
            #     nearfield_condensed.row_ptrs
            # )).tocsc()

            # print(np.linalg.cond(np_matrix.todense()))
            # nnz = np_matrix.nnz
            # max_nnz = np_matrix.shape[0] * np_matrix.shape[1]
            # density = nnz / float(max_nnz)
            # print('matrix density ' + str(density))

            # P = scipy.sparse.linalg.splu(np_matrix)
            def precondition(x):
                return x
                print("PREC")
                return P.solve(x)

            def mv(v):
                distributed = tbem.distribute_vector(
                    constraint_matrix, v, op.n_total_rows()
                )
                applied = op.apply(distributed)
                condensed = tbem.condense_vector(constraint_matrix, applied)
                self.iterations += 1
                return condensed
            n = rhs.shape[0]
            M = scipy.sparse.linalg.LinearOperator((n, n), precondition)
            A = scipy.sparse.linalg.LinearOperator((n, n), matvec = mv,
                dtype = np.float64)

            def callback(residual):
                print(residual)

            res = scipy.sparse.linalg.gmres(
                A, rhs, tol = 1e-7,
                callback = callback, M = M
            )
            assert(res[1] == 0) #Check that the iterative solver succeeded
            return res[0]

    solver = Solver()

    bdry_error, int_error = run(
        solver.solve, dense_integral_operator, 10, log_u, make_log_dudn, 3.0
    )
    assert(bdry_error < 2e-5)
    assert(int_error < 2e-5)
    assert(solver.iterations < 20)

if __name__ == "__main__":
    test_ILU()
