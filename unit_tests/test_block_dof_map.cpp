#include "catch.hpp"
#include "block_dof_map.h"
#include "vectorx.h"
#include "util.h"

using namespace tbem;

TEST_CASE("BuildBlockDOFMap", "[block_dof_map]") {
    auto dof_map = build_block_dof_map({1,2,3,4});
    REQUIRE_ARRAY_EQUAL(dof_map.start_positions, (std::vector<size_t>{0,1,3,6}), 4);
    REQUIRE(dof_map.n_dofs == 10);
    REQUIRE(dof_map.n_components == 4);
}

TEST_CASE("FromVectorXs", "[block_dof_map]") {
    BlockVectorX input{{1,2}, {3,4,5}, {6}};
    auto dof_map = block_dof_map_from_functions(input);
    REQUIRE(dof_map.n_components == 3);
    REQUIRE(dof_map.n_dofs == 6);
    REQUIRE_ARRAY_EQUAL(dof_map.start_positions, (std::vector<size_t>{0,2,5}), 3);
}

TEST_CASE("Concatenate", "[block_dof_map]") {
    BlockVectorX input{{1,2}, {3,4,5}, {6}};
    auto dof_map = block_dof_map_from_functions(input);
    auto concat_fnc = concatenate(dof_map, input);
    REQUIRE_ARRAY_EQUAL(concat_fnc, (VectorX{1,2,3,4,5,6}), 6);
}

TEST_CASE("Expand", "[block_dof_map]") {
    VectorX input{1,2,3,4,5};
    auto dof_map = build_block_dof_map({1,3,1});
    auto expanded = expand(dof_map, input);
    REQUIRE_ARRAY_EQUAL(expanded[0], (VectorX{1}), 1);
    REQUIRE_ARRAY_EQUAL(expanded[1], (VectorX{2,3,4}), 3);
    REQUIRE_ARRAY_EQUAL(expanded[2], (VectorX{5}), 1);
}

bool expand_concat_identity(const std::vector<std::vector<double>>& A) {
    BlockVectorX f(A.size());
    for (size_t i = 0; i < A.size(); i++) {
        f[i] = VectorX(A[i]);
    }
    auto dof_map = block_dof_map_from_functions(f);
    auto block_fnc = concatenate(dof_map, f); 
    auto expanded = expand(dof_map, block_fnc);
    if (A.size() != expanded.size()) {
        return false;
    }
    for (size_t i = 0; i < A.size(); i++) {
        if (A[i].size() != expanded[i].size()) {
            return false;
        }
        for (size_t j = 0; j < A[i].size(); j++) {
            if (std::fabs(A[i][j] - expanded[i][j]) > 1e-12) {
                return false;
            }
        }
    }
    return true;
}

TEST_CASE("ExpandConcatProperty", "[block_dof_map]") {
    for (size_t i = 0; i < 100; i++) {
        auto inner_size = random<size_t>(0, 10);
        std::vector<std::vector<double>> input;
        for (size_t j = 0; j < inner_size; j++) {
            input.push_back(random_list(random<size_t>(0, 10)));
        }
        expand_concat_identity(input);
    }
}
