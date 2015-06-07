#ifndef ASDJLKJKLQJJQQJQJQJQ_LIMIT_DIRECTION_H
#define ASDJLKJKLQJJQQJQJQJQ_LIMIT_DIRECTION_H

#include <cstdlib>
#include <vector>
#include "vec.h"
#include "octree.h"

namespace tbem {

template <size_t dim> struct Ball;

template <size_t dim>
struct NearfieldFacets {
    const std::vector<Vec<Vec<double,dim>,dim>> facets;
    const Vec<Vec<double,dim>,dim> nearest_facet;
    const Vec<double,dim> pt;
    const double distance;
};

template <size_t dim>
struct NearfieldFacetFinder {
    const static size_t n_facets_per_leaf = 20;
    const std::vector<Vec<Vec<double,dim>,dim>> facets;
    const std::vector<Ball<dim>> facet_balls;
    const Octree<dim> oct;

    NearfieldFacetFinder(const std::vector<Vec<Vec<double,dim>,dim>>& facets);

    NearfieldFacets<dim> find(const Vec<double,dim>& pt);
};


// A safety factor of 0.4 should work for all cases.
template <size_t dim>
Vec<double,dim> decide_limit_dir(const Vec<double,dim>& pt,
    const NearfieldFacets<dim>& nearest_facets, double safety_factor,
    double epsilon = 1e-12);

}

#endif
