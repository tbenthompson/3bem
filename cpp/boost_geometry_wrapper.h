#ifndef ZXCVBNMLKJHGFDSAQWE_BOOST_GEOMETRY_WRAPPER_H
#define ZXCVBNMLKJHGFDSAQWE_BOOST_GEOMETRY_WRAPPER_H

#include <vector>
#include "vec.h"

namespace tbem {

std::vector<Vec<double,2>> seg_seg_intersection(const Vec<Vec<double,2>,2>& A,
    const Vec<Vec<double,2>,2>& B);

std::vector<Vec<double,3>> seg_tri_intersection(const Vec<Vec<double,3>,3>& A,
    const Vec<Vec<double,3>,2>& B);

template <size_t dim> 
std::vector<Vec<double,dim>> seg_facet_intersection(const Vec<Vec<double,dim>,dim>& f,
    const Vec<Vec<double,dim>,2>& seg);

bool in_polygon(const std::vector<Vec<double,2>>& poly, const Vec<double,2>& pt);

}

#endif
