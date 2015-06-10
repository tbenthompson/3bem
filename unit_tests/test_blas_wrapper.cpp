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

TEST_CASE("Thresholded pseudoinverse", "[blas_wrapper]") 
{
    // Matrix has two singular values: 1.0 and 1e-5
    std::vector<double> matrix{
        0.0238032718573239, 0.1524037864980028,
        0.1524037864980028, 0.9762067281426762
    };
    auto svd = svd_decompose(matrix);
    auto no_threshold_pseudoinv = svd_pseudoinverse(svd);
    std::vector<double> correct_no_threshold{
        97620.6728142285282956, -15240.3786497941800917,
        -15240.3786497941782727, 2380.3271857314393856
    };
    REQUIRE_ARRAY_CLOSE(no_threshold_pseudoinv, correct_no_threshold, 4, 1e-4);
    set_threshold(svd, 1e-4);
    auto thresholded_pseudoinv = svd_pseudoinverse(svd);
    std::vector<double> correct_thresholded{
        0.0237935097924219, 0.1524053105511083,
        0.1524053105511085, 0.9762064902075779
    };
    REQUIRE_ARRAY_CLOSE(thresholded_pseudoinv, correct_thresholded, 4, 1e-12);
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

TEST_CASE("matrix vector product", "[blas_wrapper]")
{
    
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto result = matrix_vector_product({2, 1, -1, 0.5}, {4, -2});
    REQUIRE_ARRAY_CLOSE(result, std::vector<double>{6, -5}, 2, 1e-15);
}

TEST_CASE("matrix vector non-square", "[blas_wrapper]")
{
    std::vector<double> matrix{
        2, 1, 1, -1, 0.5, 10
    };
    auto result = matrix_vector_product(matrix, {4, -2, 0.5});
    REQUIRE_ARRAY_CLOSE(result, std::vector<double>{6.5, 0}, 2, 1e-15);
}
