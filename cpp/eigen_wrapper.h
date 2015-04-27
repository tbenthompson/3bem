#ifndef __ZBZBBZBZBZBZBZBBZBZ_BLAS_WRAPPER_H
#define __ZBZBBZBZBZBZBZBBZBZ_BLAS_WRAPPER_H

#include <vector>
#include <memory>

namespace tbem {

struct LU;
struct LUDeleter {
    void operator()(LU* thing);
};
typedef std::unique_ptr<LU,LUDeleter> LUPtr;

LUPtr LU_decompose(const std::vector<double>& matrix);
std::vector<double> LU_solve(const LUPtr& lu, const std::vector<double>& b);
double condition_number(const std::vector<double>& matrix); 

struct SVD {
    
};

/* Form the SVD. Assumes a square matrix. */
SVD svd_decompose(const std::vector<double>& matrix);

} // end namespace tbem

#endif
