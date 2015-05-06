#include "blas_wrapper.h"
#include <cmath>
#include <cassert>

extern "C" void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv,
    int* info);
extern "C" void dgetrs_(char* TRANS, int* N, int* NRHS, double* A, int* LDA,
    int* IPIV, double* B, int* LDB, int* INFO);
extern "C" void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A,
    int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK,
    int* LWORK, int* INFO);
extern "C" void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, 
    double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA,
    double* C, int* LDC);

namespace tbem {

struct LU {
    std::vector<double> LU;
    std::vector<int> pivots;
};

void LUDeleter::operator()(LU* thing) {
    delete thing;
}


LUPtr lu_decompose(const std::vector<double>& matrix) 
{
    int n = std::sqrt(matrix.size()); 
    std::vector<double> A = matrix;
    std::vector<int> pivots(n);
    int info;
    dgetrf_(&n, &n, A.data(), &n, pivots.data(), &info);
    assert(info == 0);
    return LUPtr(new LU{A, pivots});
}

std::vector<double> lu_solve(const LUPtr& lu, const std::vector<double>& b) 
{
    std::vector<double> x = b;
    int n = b.size(); 
    int n_rhs = 1;
    int info;
    // type = 'T' means that we solve A^T x = b rather than Ax = b. This is good
    // because blas operates on column major data while this code sticks to 
    // row major data.
    char type = 'T';
    dgetrs_(&type, &n, &n_rhs, lu->LU.data(), &n, lu->pivots.data(), x.data(), &n, &info);
    assert(info == 0);
    return x;
}

struct SVD {
    std::vector<double> singular_values;
    std::vector<double> left_singular_vectors;
    std::vector<double> right_singular_vectors;
    double threshold;
};

void SVDDeleter::operator()(SVD* thing) {
    delete thing;
}

SVDPtr svd_decompose(const std::vector<double>& matrix) 
{
    std::vector<double> A = matrix;
    int n = std::sqrt(matrix.size()); 
    char jobu = 'A';
    char jobvt = 'A';
    auto svd = SVDPtr(new SVD{
        std::vector<double>(n),
        std::vector<double>(n * n),
        std::vector<double>(n * n),
        1e-16
    });
    int lwork = 5 * n * n;
    std::vector<double> work_space(lwork);
    int info;
    dgesvd_(
        &jobu, &jobvt, &n, &n, A.data(), &n, svd->singular_values.data(),
        svd->left_singular_vectors.data(), &n, svd->right_singular_vectors.data(),
        &n, work_space.data(), &lwork, &info
    );
    assert(info == 0);
    return std::move(svd);
}

void set_threshold(const SVDPtr& svd, double threshold) 
{
    svd->threshold = threshold;
}

std::vector<double> svd_pseudoinverse(const SVDPtr& svd)
{
    auto n = svd->singular_values.size();
    std::vector<double> left_singular_vectors = svd->left_singular_vectors;

    double cutoff = svd->threshold * svd->singular_values[0];
    for (size_t i = 0; i < n; i++) {
        if (svd->singular_values[i] > cutoff) {
            for (size_t j = 0; j < n; j++) {
                left_singular_vectors[i * n + j] /= svd->singular_values[i];
            }
        } else {
            for (size_t j = 0; j < n; j++) {
                left_singular_vectors[i * n + i] = 0.0;
            }
        }
    }
    
    std::vector<double> out(n * n);
    char transa = 'T';
    char transb = 'T';
    int mn = static_cast<int>(n);
    double alpha = 1;
    double beta = 0;
    dgemm_(
        &transa, &transb, &mn, &mn, &mn, &alpha,
        svd->right_singular_vectors.data(), &mn,
        left_singular_vectors.data(), &mn,
        &beta, out.data(), &mn
    );
    return out;
}

std::vector<double> svd_solve(const SVDPtr& svd, const std::vector<double>& b)
{
    // With SVD = U(S)V^T
    // SVD inverse = V(S^{-1})U^T
    // But BLAS input is column-major while my input is row-major, so I want
    // (A^T)^-1 = (V(S^{-1})U^T)^T 
    //          = U(S^{-1})V^T
    auto n = b.size();
    std::vector<double> mult_u_transpose(n, 0.0);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            mult_u_transpose[i] += svd->right_singular_vectors[i * n + j] * b[j];
        }
    }

    double cutoff = svd->threshold * svd->singular_values[0];
    for (size_t i = 0; i < n; i++) {
        if (svd->singular_values[i] > cutoff) {
            mult_u_transpose[i] /= svd->singular_values[i];
        } else {
            mult_u_transpose[i] = 0.0;
        }
    }

    std::vector<double> out(n, 0.0);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            out[i] += svd->left_singular_vectors[i * n + j] * mult_u_transpose[j];
        }
    }

    return out;
}

double condition_number(const SVDPtr& svd)
{
    auto first = svd->singular_values.front();
    auto last = svd->singular_values.back();
    return first / last;
}

}// end namespace tbem
