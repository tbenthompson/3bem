#include "benchmark/benchmark.h"
#include "octree.h"
#include "geometry.h"
#include "util.h"

#include "fmm.h"
#include "nbody_operator.h"
#include "elastic_kernels.h"

using namespace tbem;

static void construct_octree(benchmark::State& state)
{
    int n = state.range_x();
    auto pts = random_pts<3>(n);
    while (state.KeepRunning())
    {
        auto oct = make_octree(pts, 20); 
    }
}
BENCHMARK(construct_octree)->Range(1, 250000);

template <size_t dim>
static void fmm_hypersingular_apply(benchmark::State& state)
{
    int n = state.range_x();
    size_t n_per_cell = state.range_y();
    size_t order = 30;
    if (dim == 3) {
        order = 100;
    }

    auto src_pts = random_pts<dim>(n);
    auto obs_pts = random_pts<dim>(n);
    auto normals = random_pts<dim>(n);
    std::vector<double> weights(n, 1.0);
    NBodyData<dim> data{obs_pts, normals, src_pts, normals, weights};
    std::vector<double> x(dim * data.src_locs.size(), 1.0);

    FMMOperator<dim,dim,dim> tree(
        ElasticHypersingular<dim>(30e9, 0.25),
        data,
        {0.35, order, n_per_cell, 0.05, true}
    );

    while (state.KeepRunning())
    {
        auto out = tree.apply(x);
    }
}
BENCHMARK_TEMPLATE(fmm_hypersingular_apply, 2)
    ->ArgPair(1000, 20)
    ->ArgPair(1000, 60)
    ->ArgPair(1000, 100)
    ->ArgPair(1000, 250)
    ->ArgPair(10000, 20)
    ->ArgPair(10000, 60)
    ->ArgPair(10000, 100)
    ->ArgPair(10000, 250)
    ->ArgPair(100000, 20)
    ->ArgPair(100000, 60)
    ->ArgPair(100000, 100)
    ->ArgPair(100000, 250)
    ->ArgPair(500000, 20)
    ->ArgPair(500000, 60)
    ->ArgPair(500000, 100)
    ->ArgPair(500000, 250);

BENCHMARK_TEMPLATE(fmm_hypersingular_apply, 3)
    // ->ArgPair(1000, 20)
    // ->ArgPair(1000, 60)
    // ->ArgPair(1000, 100)
    // ->ArgPair(1000, 250)
    // ->ArgPair(10000, 20)
    // ->ArgPair(10000, 60)
    ->ArgPair(10000, 100)
    ->ArgPair(10000, 250)
    ->ArgPair(100000, 20)
    ->ArgPair(100000, 60)
    ->ArgPair(100000, 100)
    ->ArgPair(100000, 250);
    // ->ArgPair(500000, 20)
    // ->ArgPair(500000, 60)
    // ->ArgPair(500000, 100)
    // ->ArgPair(500000, 250);

BENCHMARK_MAIN()
