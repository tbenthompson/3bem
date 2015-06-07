#ifndef ZXCVBNMLKJHGFDSAQWE_GTE_WRAPPER_H
#define ZXCVBNMLKJHGFDSAQWE_GTE_WRAPPER_H

#include <vector>
#include "vec.h"

namespace tbem {

template <size_t dim> 
std::vector<Vec<double,dim>> seg_facet_intersection(const Vec<Vec<double,dim>,dim>& f,
    const Vec<Vec<double,dim>,2>& seg);

bool in_polygon(const std::vector<Vec<double,2>>& poly, const Vec<double,2>& pt);

template <size_t dim>
Vec<double,dim-1> closest_pt_facet(const Vec<double,dim>& pt,
    const Vec<Vec<double,dim>,dim> tri);

}

#endif
