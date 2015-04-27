#include "eigen_wrapper.h"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace tbem {

typedef Eigen::Matrix<double,-1,-1,Eigen::RowMajor> EigenMatrix;

struct LU {
    Eigen::PartialPivLU<EigenMatrix> lu;
};

void LUDeleter::operator()(LU* thing) {
    delete thing;
}


LUPtr lu_decompose(const std::vector<double>& matrix) 
{
    int n = std::sqrt(matrix.size()); 
    Eigen::Map<const EigenMatrix> eigen_mat(matrix.data(), n, n);
    auto decomp = eigen_mat.lu();
    return LUPtr(new LU{decomp});
}

std::vector<double> lu_solve(const LUPtr& lu, const std::vector<double>& b) 
{
    Eigen::Map<const Eigen::VectorXd> eigen_b(b.data(), lu->lu.cols());
    Eigen::VectorXd soln = lu->lu.solve(eigen_b);
    return std::vector<double>(soln.data(), soln.data() + lu->lu.rows());
}

struct SVD {
    Eigen::JacobiSVD<EigenMatrix> svd;
};

void SVDDeleter::operator()(SVD* thing) {
    delete thing;
}

SVDPtr svd_decompose(const std::vector<double>& matrix) 
{
    int n = std::sqrt(matrix.size()); 
    Eigen::Map<const EigenMatrix> eigen_mat(matrix.data(), n, n);
    Eigen::JacobiSVD<EigenMatrix> svd(
        eigen_mat, Eigen::ComputeThinU | Eigen::ComputeThinV
    );
    return SVDPtr(new SVD{svd});
}

void set_threshold(const SVDPtr& svd, double threshold) 
{
    svd->svd.setThreshold(threshold);
}

std::vector<double> svd_solve(const SVDPtr& svd, const std::vector<double>& b)
{
    Eigen::Map<const Eigen::VectorXd> eigen_b(b.data(), svd->svd.cols());
    Eigen::VectorXd soln = svd->svd.solve(eigen_b);
    return std::vector<double>(soln.data(), soln.data() + svd->svd.rows());
}

double condition_number(const SVDPtr& svd)
{
    auto first = svd->svd.singularValues()[0];
    auto last = svd->svd.singularValues()[svd->svd.rows() - 1];
    return first / last;
}

}// end namespace tbem
