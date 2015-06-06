#ifndef ASDJLKJKLQJJQQJQJQJQ_LIMIT_DIRECTION_H
#define ASDJLKJKLQJJQQJQJQJQ_LIMIT_DIRECTION_H

#include <cstdlib>
#include "vec.h"

namespace tbem {

//TODO: feature envy -- this structure and its respective builder should be
//in this header/compilation unit
template <size_t dim> struct NearestFacets;

// TODO: safety_factor = 0.5 should work, but doesn't. Why?
template <size_t dim>
Vec<double,dim> decide_limit_dir(const Vec<double,dim>& pt,
    const NearestFacets<dim>& nearest_facets, double epsilon = 1e-12,
    double safety_factor = 0.5);

}

#endif
