#include "petsc_facade.h"
#include "petscksp.h"
#include "matrix_entry.h"
#include "vectorx.h"
#include <algorithm>
#include <iostream>

namespace tbem {

std::vector<int> build_sparsity_pattern(size_t n_rows, 
    const std::vector<MatrixEntry>& entries) 
{
    std::vector<int> row_nnz(n_rows, 0.0);
    for (const auto& e: entries) {
        row_nnz[e.loc[0]] += 1;
    }
    return row_nnz;
}

void init_petsc() {
    int argc = 0;
    char** args = new char*[0];
    //TODO: Remove PetscInitialize from in here. It will be needed outside I think.
    PetscInitialize(&argc, &args, (char*)0, "bem code"); 
}

PETScSparseMatWrapper::PETScSparseMatWrapper(size_t n_rows, size_t n_cols,
    const std::vector<MatrixEntry>& entries)
{
    init_petsc();

    MPI_Comm comm = PETSC_COMM_WORLD;
    MatCreate(comm, &internal_mat);
    MatSetType(internal_mat, MATSEQAIJ); 
    MatSetSizes(internal_mat, PETSC_DECIDE, PETSC_DECIDE, n_rows, n_cols);

    auto sparsity_pattern = build_sparsity_pattern(n_rows, entries);
    MatSeqAIJSetPreallocation(internal_mat, 0, sparsity_pattern.data());

    for (const auto& e: entries) {
        int row = e.loc[0];
        int col = e.loc[1];
        MatSetValues(internal_mat, 1, &row, 1, &col, &e.value, ADD_VALUES);
    }

    MatAssemblyBegin(internal_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(internal_mat, MAT_FINAL_ASSEMBLY);
}

std::pair<int,int> get_matrix_shape(Mat m) {
    std::pair<int,int> shape;
    MatGetSize(m, &shape.first, &shape.second);
    return shape;
}

size_t PETScSparseMatWrapper::n_rows() const {
    return static_cast<size_t>(get_matrix_shape(internal_mat).first);
}

size_t PETScSparseMatWrapper::n_cols() const {
    return static_cast<size_t>(get_matrix_shape(internal_mat).second);
}

VectorX PETScSparseMatWrapper::mat_vec_prod(const VectorX& v) const {
    MPI_Comm comm = PETSC_COMM_WORLD;

    Vec x;
    VecCreateMPIWithArray(comm, v.size(), v.size(), PETSC_DECIDE, v.data(), &x);

    VectorX out(n_cols());
    Vec out_petsc;
    VecCreateMPIWithArray(comm, out.size(), out.size(), PETSC_DECIDE, out.data(), &out_petsc);
    MatMult(internal_mat, x, out_petsc);

    return out;
}

/* Useful for debugging PETSc code.
 */
void print_vec(Vec& x) {
    int size;
    VecGetSize(x, &size);
    double* vals;
    VecGetArray(x, &vals);
    for (int i = 0; i < size; i++) {
        std::cout << vals[i] << std::endl;
    }
    VecRestoreArray(x, &vals);
}

template <class T>
void wrapArrayInVector(T *sourceArray, size_t arraySize, std::vector<T> &targetVector) 
{
    typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *vectorPtr =
        (typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *)(
            (void *) &targetVector);
    vectorPtr->_M_start = sourceArray;
    auto finish = vectorPtr->_M_start + arraySize;
    vectorPtr->_M_finish = finish;
    vectorPtr->_M_end_of_storage = finish;
}

template <class T>
void releaseVectorWrapper(std::vector<T> &targetVector) 
{
    typename std::_Vector_base<T, std::allocator<T>>::_Vector_impl* vectorPtr =
        (typename std::_Vector_base<T, std::allocator<T>>::_Vector_impl*)(
            (void*) &targetVector);
    vectorPtr->_M_start = NULL;
    vectorPtr->_M_finish = NULL;
    vectorPtr->_M_end_of_storage = NULL;
}

PetscErrorCode mult_wrapper(Mat A, Vec x, Vec y) {
    MatVecFnc* fnc;
    MatShellGetContext(A, (void**) &fnc);

    // Extract the data from the Petsc Vec objects
    double* x_vals;
    double* y_vals;
    int x_size;
    int y_size;
    VecGetSize(x, &x_size);
    VecGetSize(y, &y_size);
    VecGetArray(x, &x_vals);
    VecGetArray(y, &y_vals);

    // Place the data in an std::vector
    std::vector<double> x_wrap;
    std::vector<double> y_wrap;
    wrapArrayInVector(x_vals, x_size, x_wrap);
    wrapArrayInVector(y_vals, y_size, y_wrap);
    (*fnc)(x_wrap, y_wrap);

    // Remove the data from the std::vector
    releaseVectorWrapper(x_wrap);
    releaseVectorWrapper(y_wrap);

    // Replace the data in the Petsc Vec objects
    VecRestoreArray(x, &x_vals);
    VecRestoreArray(y, &y_vals);

    return 0;
}

void setup_ksp(MPI_Comm comm, KSP& ksp, Mat& mat, double tolerance) {
    PC pc;
    PetscErrorCode ierr; 
    int maxiter = 1000;

    ierr = KSPCreate(comm, &ksp);
    CHKERRABORT(comm,ierr);
    ierr = KSPSetOperators(ksp, mat, mat, DIFFERENT_NONZERO_PATTERN);
    CHKERRABORT(comm,ierr);
    ierr = KSPGetPC(ksp, &pc);
    CHKERRABORT(comm,ierr);
    //ierr = PCSetType(pc, ??????);
    CHKERRABORT(comm,ierr);
    ierr = KSPSetTolerances(ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, maxiter);
    CHKERRABORT(comm,ierr);
}

std::vector<double>
solve_system(const double* rhs, int n_dofs, double tolerance, MatVecFnc fnc) 
{
    init_petsc();

    PetscErrorCode ierr;
    MPI_Comm comm = PETSC_COMM_WORLD;

    Vec b;
    VecCreateMPIWithArray(comm, n_dofs, n_dofs, PETSC_DECIDE, rhs, &b);

    // USE VecCreateMPIWithArray to set the data location in the PETSc vector
    // so that I don't need to copy.
    std::vector<double> retval(n_dofs, 0.0);
    Vec xo;
    VecCreateMPIWithArray(comm, 1, n_dofs, PETSC_DECIDE, retval.data(), &xo);

    Mat mat;
    MatCreateShell(comm, n_dofs, n_dofs, PETSC_DETERMINE, PETSC_DETERMINE, &fnc, &mat);
    MatShellSetOperation(mat, MATOP_MULT, (void(*)(void))mult_wrapper);

    KSP ksp;
    setup_ksp(comm, ksp, mat, tolerance);

    ierr = KSPSolve(ksp, b, xo);CHKERRABORT(comm,ierr);
    return retval;
}

std::vector<double> solve_system(const VectorX& rhs, double tolerance, MatVecFnc fnc) 
{
    return solve_system(rhs.data(), rhs.size(), tolerance, fnc);
}

} // END namespace tbem
