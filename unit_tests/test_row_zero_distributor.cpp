#include "catch.hpp"
#include "row_zero_distributor.h"
#include "constraint_matrix.h"
#include "dense_operator.h"

using namespace tbem;

TEST_CASE("identify ignored dofs", "[row_zero_distributor]")
{
    auto cm = from_constraints({
        boundary_condition(0, 5),
        boundary_condition(2, 10),
        continuity_constraint(1, 4)
    });
    auto constrained = identify_ignored_dofs(cm);
    REQUIRE(constrained.size() == 2);
    REQUIRE_ARRAY_EQUAL(constrained, std::vector<double>{0, 2}, 2);
}

TEST_CASE("identify ignored dofs none", "[row_zero_distributor]")
{
    auto cm = from_constraints({
        continuity_constraint(2, 3),
        continuity_constraint(4, 3),
        continuity_constraint(1, 4)
    });
    auto constrained = identify_ignored_dofs(cm);
    REQUIRE(constrained.size() == 0);
}

TEST_CASE("identify ignored dofs with dependencies", "[row_zero_distributor]")
{
    auto cm = from_constraints({
        boundary_condition(0, 5),
        boundary_condition(2, 10),
        continuity_constraint(2, 3),
        continuity_constraint(3, 5)
    });
    auto constrained = identify_ignored_dofs(cm);
    REQUIRE(constrained.size() == 4);
    REQUIRE_ARRAY_EQUAL(constrained, std::vector<double>{0, 2, 3, 5}, 4);
}

TEST_CASE("distribute row zeros", "[row_zero_distributor]")
{
    std::vector<double> matrix{
        {1,0,0  ,  0,1,0  ,  0,0,1}
    };
    auto result = distribute_row_zeros(
        DenseOperator(3, 3, matrix),
        from_constraints({
            boundary_condition(0, 5),
            boundary_condition(1, 10),
        })
    );
    REQUIRE(result.n_rows() == 5);
    REQUIRE(result.n_cols() == 3);
    std::vector<double> correct{
        {0,0,0, 0,0,0, 1,0,0, 0,1,0, 0,0,1}
    };
    REQUIRE_ARRAY_EQUAL(result, correct, 15);
}

TEST_CASE("row zero distributor", "[row_zero_distributor]")
{
    DenseOperator op(3, 3, {{1,0,0  ,  0,1,0  ,  0,0,1}});
    auto cm = from_constraints({boundary_condition(0, 5), boundary_condition(1, 10)});
    RowZeroDistributor rzd(cm, op);
    auto result = rzd.apply({1, 1, 1});
    std::vector<double> correct{{0, 0, 1, 1, 1}};
    REQUIRE_ARRAY_EQUAL(result, correct, 5);
}
