#include "catch.hpp"
#include "mesh_converter.h"

using namespace tbem;

TEST_CASE("facet_to_pt_index simple", "[mesh_converter]")
{
    Mesh<2> in{
        {{{{0, 0}, {1, 0}}}, {{{1, 0}, {0, 1}}}, {{{0, 1}, {0, 0}}}}
    };
    auto out = convert_facet_to_pt_index<2>(in);

    REQUIRE(out.points.size() == 3);
    REQUIRE_ARRAY_EQUAL(out.points[0], std::vector<double>{0, 0}, 2);
    REQUIRE_ARRAY_EQUAL(out.points[1], std::vector<double>{1, 0}, 2);
    REQUIRE_ARRAY_EQUAL(out.points[2], std::vector<double>{0, 1}, 2);

    REQUIRE(out.facets.size() == 3);
    REQUIRE_ARRAY_EQUAL(out.facets[0], std::vector<size_t>{0, 1}, 2);
    REQUIRE_ARRAY_EQUAL(out.facets[1], std::vector<size_t>{1, 2}, 2);
    REQUIRE_ARRAY_EQUAL(out.facets[2], std::vector<size_t>{2, 0}, 2);
}

TEST_CASE("facet_to_pt_index harder", "[mesh_converter]")
{
    Mesh<2> in{{
        {{{0, 0}, {1, 0}}}, {{{1, 0}, {1, 1}}}, 
        {{{1, 1}, {0, 1}}}, {{{0, 1}, {0, 0}}}, {{{0, 0}, {1, 1}}}}
    };
    auto out = convert_facet_to_pt_index<2>(in);
    REQUIRE(out.points.size() == 4);
    REQUIRE(out.facets.size() == 5);
}
