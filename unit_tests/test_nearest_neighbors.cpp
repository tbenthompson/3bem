#include "catch.hpp"
#include "nearest_neighbors.h"
#include "util.h"
#include "mesh.h"

using namespace tbem;

TEST_CASE("nearest facet brute force", "[nearest_neighbors]")
{
    auto result = nearest_facet_brute_force<2>({0, 0}, {
        {{{0, 1}, {0, 2}}}, {{{1, 1}, {1, 2}}}
    });
    REQUIRE(result.idx == 0);
    REQUIRE(result.distance == 1.0);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0, 1}, 2, 1e-12);
}
// 
// TEST_CASE("test nearest facet", "[nearest_neighbors]") 
// {
//     //Test the fast version against the brute force version through many
//     //random examples
//     size_t n_facets = 10;
//     for (size_t i = 0; i < 100; i++) {
//         std::vector<Facet<2>> fs;
//         for (size_t j = 0; j < n_facets; j++) {
//             auto vs = random_pts<2>(2);
//             fs.push_back({vs[0], vs[1]});
//         }
//         // auto brute_force = nearest_facet_brute_force(
//     }
// }

