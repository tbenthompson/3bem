#include "catch.hpp"
#include "sparse_operator.h"

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
    CHECK_ARRAY_EQUAL(d, correct, 9);
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
    CHECK_ARRAY_EQUAL(d, correct, 9);
}
