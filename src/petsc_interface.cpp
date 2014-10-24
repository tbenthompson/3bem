#include "petsc_interface.h"
#include "petscksp.h"
#include <algorithm>
#include <iostream>

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

int setup_vec(Vec& x, int n, MPI_Comm& comm) {
    PetscErrorCode ierr;
    ierr = VecCreate(comm, &x);CHKERRQ(ierr);
    ierr = VecSetSizes(x, n, PETSC_DECIDE);CHKERRQ(ierr);
    ierr = VecSetType(x, "mpi"); CHKERRQ(ierr);
    return ierr;
}

template <class T>
void wrapArrayInVector( T *sourceArray, size_t arraySize, std::vector<T, std::allocator<T> > &targetVector ) {
    typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *vectorPtr =
    (typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *)((void *) &targetVector);
    vectorPtr->_M_start = sourceArray;
    vectorPtr->_M_finish = vectorPtr->_M_end_of_storage = vectorPtr->_M_start + arraySize;
}

template <class T>
void releaseVectorWrapper( std::vector<T, std::allocator<T> > &targetVector ) {
    typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *vectorPtr =
        (typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *)((void *) &targetVector);
    vectorPtr->_M_start = vectorPtr->_M_finish = vectorPtr->_M_end_of_storage = NULL;
}

template <typename T>
PetscErrorCode mult_wrapper(Mat A, Vec x, Vec y) {
    T* fnc;
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

    ierr = KSPCreate(comm, &ksp);
    CHKERRABORT(comm,ierr);
    ierr = KSPSetOperators(ksp, mat, mat);
    CHKERRABORT(comm,ierr);
    ierr = KSPGetPC(ksp, &pc);
    CHKERRABORT(comm,ierr);
    //ierr = PCSetType(pc, ??????);
    CHKERRABORT(comm,ierr);
    ierr = KSPSetTolerances(ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRABORT(comm,ierr);
}

std::vector<double> solve_system(std::vector<double> rhs,
                                 double tolerance,
                                 MatVecFnc fnc) {
    int argc = 0;
    char** args = new char*[1];
    args[1] = (char*)"./test/test_petsc";
    PetscInitialize(&argc, &args, (char*)0, "bem code"); 

    unsigned int n = rhs.size();
    std::cout << "Vector length: " << n << std::endl;
    PetscErrorCode ierr;
    MPI_Comm comm = PETSC_COMM_WORLD;

    Vec b;
    VecCreateMPIWithArray(comm, n, n, PETSC_DECIDE, rhs.data(), &b);

    Mat mat;
    MatCreateShell(comm, n, n, PETSC_DETERMINE, PETSC_DETERMINE, &fnc, &mat);
    MatShellSetOperation(mat, MATOP_MULT, (void(*)(void))mult_wrapper<MatVecFnc>);

    KSP ksp;
    setup_ksp(comm, ksp, mat, tolerance);

    // USE VecCreateMPIWithArray to set the data location in the PETSc vector
    // so that I don't need to copy.
    std::vector<double> retval(n, 0.0);
    Vec xo;
    VecCreateMPIWithArray(comm, 1, n, PETSC_DECIDE, retval.data(), &xo);

    ierr = KSPSolve(ksp, b, xo);CHKERRABORT(comm,ierr);
    return retval;
}
