#ifndef TBEMASHKLHLA_CLOSEST_PT_H
#define TBEMASHKLHLA_CLOSEST_PT_H

#include "vec.h"

namespace tbem {

Vec<double,1> closest_pt_seg(const Vec<double,2>& pt, const Vec<Vec<double,2>,2> tri);
Vec<double,2> closest_pt_tri(const Vec<double,3>& pt, const Vec<Vec<double,3>,3> tri);

template <size_t dim>
Vec<double,dim-1> closest_pt_facet(const Vec<double,dim>& pt,
    const Vec<Vec<double,dim>,dim> tri);

} // end namespace tbem

#endif
