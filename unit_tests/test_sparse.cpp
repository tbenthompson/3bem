#include "catch.hpp"
#include "sparse_operator.h"
#include "dense_operator.h"

using namespace tbem;

TEST_CASE("Empty sparse apply", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(2, 2, {});
    auto out = op.apply({0.5, 4.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.0, 0.0}, 2);
}

TEST_CASE("Sparse apply", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(2, 2, {
        {0, 0, 1.0}, {1, 0, 2.0}
    });
    auto out = op.apply({0.5, 0.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.5, 1.0}, 2);
}

TEST_CASE("Sparse apply different", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(3, 3, {
        {0, 0, 1.0}, {1, 0, 2.0}, {1, 1, 4.0}
    });
    auto out = op.apply({0.5, 7.0, 1.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.5, 29.0, 0.0}, 2);
}

TEST_CASE("Last CSR row_ptr is nnz", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(3, 3, {
        {0, 0, 1.0}, {1, 0, 2.0}, {1, 1, 4.0}
    });
    REQUIRE(op.row_ptrs[op.n_rows()] == 3);
}

TEST_CASE("nnz", "[sparse]") 
{

    auto op = SparseOperator::csr_from_coo(3, 3, {
        {0, 0, 1.0}, {1, 0, 2.0}
    });
    REQUIRE(op.nnz() == 2);
}

TEST_CASE("Fewer rows than entries", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(3, 3, {
        {0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 4.0},
        {1, 0, 1.0}, {1, 1, 2.0}, {1, 2, 4.0},
        {2, 0, 1.0}, {2, 1, 2.0}, {2, 2, 4.0}
    });
}

TEST_CASE("To dense", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(3, 3, {
        {0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 4.0}, {1, 2, 4.0}, {2, 0, 1.0}, 
    });
    auto d = op.to_dense();
    std::vector<double> correct = {
        1.0, 2.0, 4.0, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0
    };
    REQUIRE_ARRAY_EQUAL(d.data(), correct, 9);
}

TEST_CASE("To dense -- multiple values for one entry", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(3, 3, {
        {0, 0, 1.0}, {0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 4.0}, {1, 2, 4.0}, {2, 0, 1.0}, 
    });
    auto d = op.to_dense();
    std::vector<double> correct = {
        2.0, 2.0, 4.0, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0
    };
    REQUIRE_ARRAY_EQUAL(d, correct, 9);
}

TEST_CASE("left multiply with dense", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(2, 2, {
        {0, 0, 1.0}, {1, 1, 2.0}, {0, 1, 1.0}
    });
    auto out = op.left_multiply_with_dense(DenseOperator(3, 2, {
        1, 2, 3, 4, 5, 6
    }));
    REQUIRE(out.n_rows() == 3);
    REQUIRE(out.n_cols() == 2);
    std::vector<double> correct{
        1, 5, 3, 11, 5, 17
    };
    REQUIRE_ARRAY_EQUAL(out.data(), correct, 4);
}

TEST_CASE("right multiply with dense", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(2, 2, {
        {0, 0, 1.0}, {1, 1, 2.0}, {0, 1, 1.0}
    });
    auto out = op.right_multiply_with_dense(DenseOperator(2, 3, {
        1, 2, 3, 4, 5, 6
    }));
    REQUIRE(out.n_rows() == 2);
    REQUIRE(out.n_cols() == 3);
    std::vector<double> correct{
        5, 7, 9, 8, 10, 12
    };
    REQUIRE_ARRAY_EQUAL(out.data(), correct, 4);
}

TEST_CASE("add with dense", "[sparse]") 
{
    auto op = SparseOperator::csr_from_coo(2, 2, {
        {0, 0, 1.0}, {1, 1, 2.0}, {0, 1, 1.0}
    });
    auto out = op.add_with_dense(DenseOperator(2, 2, {
        1, 2, 3, 4
    }));
    REQUIRE(out.n_rows() == 2);
    REQUIRE(out.n_cols() == 2);
    std::vector<double> correct{
        2, 3, 3, 6
    };
    REQUIRE_ARRAY_EQUAL(out.data(), correct, 4);
}

TEST_CASE("sparse times sparse", "[sparse]") 
{
    auto op1 = SparseOperator::csr_from_coo(2, 2, {
        {0, 0, 1.0}, {1, 1, 2.0}, {0, 1, 1.0}
    });
    auto op2 = SparseOperator::csr_from_coo(2, 4, {
        {0, 1, 1}, {0, 3, 1}, {1, 0, 7}, {1, 1, -1}
    });
    auto out = op1.right_multiply(op2);
    REQUIRE(out.n_rows() == 2);
    REQUIRE(out.n_cols() == 4);
    std::vector<double> correct{
        7, 0, 0, 1, 14, -2, 0, 0
    };
    REQUIRE_ARRAY_EQUAL(out.to_dense().data(), correct, 8);
}
