#include "catch.hpp"
#include "sparse_operator.h"

using namespace tbem;

TEST_CASE("Empty sparse apply", "[sparse]") 
{
    auto op = BlockSparseOperator::csr_from_coo(1, 1, 2, 2, {});
    auto out = op.apply({0.5, 4.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.0, 0.0}, 2);
}

TEST_CASE("Sparse apply", "[sparse]") 
{
    auto op = BlockSparseOperator::csr_from_coo(1, 1, 2, 2, {
        {0, 0, 1.0}, {1, 0, 2.0}
    });
    auto out = op.apply({0.5, 0.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.5, 1.0}, 2);
}

TEST_CASE("Block sparse apply", "[sparse]") 
{
    auto op = BlockSparseOperator::csr_from_coo(2, 2, 1, 1, {
        {0, 0, 1.0}, {1, 0, 2.0}, {1, 1, 4.0}
    });
    auto out = op.apply({0.5, 7.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.5, 29.0}, 2);
}

TEST_CASE("Block sparse apply different", "[sparse]") 
{
    auto op = BlockSparseOperator::csr_from_coo(3, 3, 1, 1, {
        {0, 0, 1.0}, {1, 0, 2.0}, {1, 1, 4.0}
    });
    auto out = op.apply({0.5, 7.0, 1.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.5, 29.0, 0.0}, 2);
}

TEST_CASE("Last CSR row_ptr is nnz", "[sparse]") 
{
    auto op = BlockSparseOperator::csr_from_coo(3, 3, 1, 1, {
        {0, 0, 1.0}, {1, 0, 2.0}, {1, 1, 4.0}
    });
    REQUIRE(op.row_ptrs[op.n_total_rows()] == 3);
}

TEST_CASE("nnz", "[sparse]") 
{

    auto op = BlockSparseOperator::csr_from_coo(3, 3, 1, 1, {
        {0, 0, 1.0}, {1, 0, 2.0}
    });
    REQUIRE(op.nnz() == 2);
}

TEST_CASE("Fewer rows than entries", "[sparse]") 
{
    auto op = BlockSparseOperator::csr_from_coo(3, 3, 1, 1, {
        {0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 4.0},
        {1, 0, 1.0}, {1, 0, 2.0}, {1, 2, 4.0},
        {2, 0, 1.0}, {2, 0, 2.0}, {2, 2, 4.0}
    });
}
