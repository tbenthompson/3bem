#include "catch.hpp"
#include "interior_mesh.h"

using namespace tbem;

TEST_CASE("pt to tri map", "[interior_mesh]")
{
    std::vector<Vec<size_t,3>> tris{{0, 1, 2}, {2, 1, 3}, {1, 5, 4}, {1, 4, 3}};
    auto result = build_pt_to_tri_map(tris);
    for (auto t_idx: std::vector<size_t>{0, 1, 2, 3}) {
        REQUIRE(result[1].count(t_idx) == 1);
    }
    for (auto t_idx: std::vector<size_t>{0, 1}) {
        REQUIRE(result[2].count(t_idx) == 1);
    }
}

TEST_CASE("test adjacent tris", "[interior_mesh]")
{
    std::vector<Vec<size_t,3>> tris{{0, 1, 2}, {2, 1, 3}, {1, 5, 4}, {1, 4, 3}};
    auto pt_tri_map = build_pt_to_tri_map(tris);
    auto tri_indices1 = find_adjacent_tris(tris, {1, 3, 4}, pt_tri_map);
    REQUIRE(tri_indices1.size() == 2);
    REQUIRE(tri_indices1.count(1) == 1);
    REQUIRE(tri_indices1.count(2) == 1);

    auto tri_indices2 = find_adjacent_tris(tris, {0, 1, 2}, pt_tri_map);
    REQUIRE(tri_indices2.size() == 1);
    REQUIRE(tri_indices2.count(1) == 1);
}

TEST_CASE("identify regions", "[interior_mesh]")
{
    std::vector<Vec<size_t,3>> tris{
        {0, 1, 2}, {2, 1, 3},
        {0, 4, 5}, {0, 5, 1},
        {1, 5, 6}, {1, 6, 3},
        {3, 6, 7}, {2, 7, 3},
        {2, 7, 4}, {2, 4, 0}
    };
    std::vector<Vec<size_t,2>> bdry{
        {0, 1}, {1, 3}, {3, 2}, {2, 0}
    };
    auto regions = identify_regions(tris, bdry);
    REQUIRE_ARRAY_EQUAL(regions, std::vector<size_t>{
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0
    }, regions.size());
}
