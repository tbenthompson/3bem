#include "catch.hpp"
#include "nearest_neighbors.h"
#include "util.h"
#include "mesh.h"

using namespace tbem;

TEST_CASE("nearest facet brute force", "[nearest_neighbors]")
{
    auto result = nearest_facet_brute_force<2>({0, 0}, NearestNeighborData<2>({
        {{{0, 1}, {0, 2}}}, {{{1, 1}, {1, 2}}}
    }, 20));
    REQUIRE(result.idx == 0);
    REQUIRE(result.distance == 1.0);
    REQUIRE_ARRAY_CLOSE(result.pt, Vec<double,2>{0, 1}, 2, 1e-12);
}

TEST_CASE("test nearest facet", "[nearest_neighbors]") 
{
    const size_t n = 100;
    std::vector<Facet<2>> fs(n);
    for (size_t i = 0; i < n; i++) {
        auto vs = random_pts<2>(2);
        auto center = (vs[0] + vs[1]) / 2.0;
        auto length = 1.0 / static_cast<double>(n);
        fs[i] = {center + (vs[0] - center) * length, center + (vs[1] - center) * length};
    }
    auto pts = random_pts<2>(n);
    NearestNeighborData<2> nn_data(fs, 1);
    for (size_t i = 0; i < n; i++) {
        auto brute = nearest_facet_brute_force(pts[i], nn_data);
        auto fast = nearest_facet(pts[i], nn_data);
        REQUIRE(brute.idx == fast.idx);
    }
}

TEST_CASE("nearest facet timing", "[nearest_neighbors]") 
{
    //TODO: Make a capacity test
    const size_t n = 5000;
    std::vector<Facet<2>> fs(n);
    for (size_t i = 0; i < n; i++) {
        auto vs = random_pts<2>(2);
        auto center = (vs[0] + vs[1]) / 2.0;
        auto length = 1.0 / static_cast<double>(n);
        fs[i] = {center + (vs[0] - center) * length, center + (vs[1] - center) * length};
    }
    auto pts = random_pts<2>(n);

    std::vector<size_t> fast(n);
    TIC
    NearestNeighborData<2> nn_data(fs, 30);
#pragma omp parallel for
    for (size_t i = 0; i < n; i++) {
        fast[i] = nearest_facet(pts[i], nn_data).idx;
    }
    TOC(std::to_string(n) + " nn-queries");
}
