#ifndef TBEMZBZBBZBZBZBZBZBBZBZ_BLAS_WRAPPER_H
#define TBEMZBZBBZBZBZBZBZBBZBZ_BLAS_WRAPPER_H

#include <vector>
#include <memory>

namespace tbem {

/* By forward declaring the LU and SVD structs, I am effectively using
 * a PImpl pattern with a very thin std::unique_ptr wrapper.
 * This allows the implementation details to be hidden in the .cpp file
 * with no leakage.
 */
struct LU;
struct LUDeleter {
    void operator()(LU* thing);
};
typedef std::unique_ptr<LU,LUDeleter> LUPtr;

/* Construct an LU decomposition assuming a square matrix. */
LUPtr lu_decompose(const std::vector<double>& matrix);
std::vector<double> lu_solve(const LUPtr& lu, const std::vector<double>& b);

struct SVD;
struct SVDDeleter {
    void operator()(SVD* thing);
};
typedef std::unique_ptr<SVD,SVDDeleter> SVDPtr;

/* Construct an SVD decomposition assuming a square matrix */
SVDPtr svd_decompose(const std::vector<double>& matrix);
void set_threshold(const SVDPtr& svd, double threshold);
std::vector<double> svd_pseudoinverse(const SVDPtr& svd);
std::vector<double> svd_solve(const SVDPtr& svd, const std::vector<double>& b);
double condition_number(const SVDPtr& svd); 

} // end namespace tbem

#endif
