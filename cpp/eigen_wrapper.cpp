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


LUPtr LU_decompose(const std::vector<double>& matrix) 
{
    int n = std::sqrt(matrix.size()); 
    Eigen::Map<const EigenMatrix> eigen_mat(matrix.data(), n, n);
    auto decomp = eigen_mat.lu();
    return LUPtr(new LU{decomp});
}

std::vector<double> LU_solve(const LUPtr& lu, const std::vector<double>& b) 
{
    Eigen::Map<const Eigen::VectorXd> eigen_b(b.data(), lu->lu.cols());
    Eigen::VectorXd soln = lu->lu.solve(eigen_b);
    return std::vector<double>(soln.data(), soln.data() + lu->lu.rows());
}

double condition_number(const std::vector<double>& matrix) 
{
    int n = std::sqrt(matrix.size()); 
    Eigen::Map<const EigenMatrix> eigen_mat(matrix.data(), n, n);
    Eigen::JacobiSVD<EigenMatrix> svd(eigen_mat);
    auto first = svd.singularValues()[0];
    auto last = svd.singularValues()[eigen_mat.rows() - 1];
    return first / last;
}

SVD svd_decompose(const std::vector<double>& matrix) 
{
    char* jobz; 
}

}// end namespace tbem
