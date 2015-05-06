#include "catch.hpp"
#include "blas_wrapper.h"

using namespace tbem;

TEST_CASE("LU solve", "[blas_wrapper]") 
{
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto lu = lu_decompose(matrix);
    auto soln = lu_solve(lu, {1,1});
    std::vector<double> correct{
        -0.25, 1.5
    };
    REQUIRE_ARRAY_CLOSE(soln, correct, 2, 1e-14);
}

TEST_CASE("SVD solve", "[blas_wrapper]") 
{
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto svd = svd_decompose(matrix);
    auto soln = svd_solve(svd, {1,1});
    std::vector<double> correct{
        -0.25, 1.5
    };
    REQUIRE_ARRAY_CLOSE(soln, correct, 2, 1e-14);
}

TEST_CASE("Pseudoinverse", "[blas_wrapper]") 
{
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto svd = svd_decompose(matrix);
    auto pseudoinv = svd_pseudoinverse(svd);
    std::vector<double> inv{
        0.25, -0.5, 0.5, 1.0
    };
    REQUIRE_ARRAY_CLOSE(pseudoinv, inv, 4, 1e-14);
}

TEST_CASE("Condition number", "[blas_wrapper]") 
{
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto svd = svd_decompose(matrix);
    double cond = condition_number(svd);
    REQUIRE(cond == Approx(2.7630857945186595).epsilon(1e-12));
}
