#include "catch.hpp"
#include "block_dof_map.h"
#include "util.h"

using namespace tbem;

TEST_CASE("BuildBlockDOFMap", "[block_dof_map]") {
    auto dof_map = build_block_dof_map({1,2,3,4});
    REQUIRE_ARRAY_EQUAL(dof_map.start_positions, (std::vector<size_t>{0,1,3,6}), 4);
    REQUIRE(dof_map.n_dofs == 10);
    REQUIRE(dof_map.n_components == 4);
}
