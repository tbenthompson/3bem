#ifndef ZXCVBNMLKJHGFDSAQWE_GTE_WRAPPER_H
#define ZXCVBNMLKJHGFDSAQWE_GTE_WRAPPER_H

#include <vector>
#include "vec.h"

namespace tbem {

template <size_t dim> 
std::vector<Vec<double,dim>> seg_facet_intersection(
    const Vec<Vec<double,dim>,dim>& f, const Vec<Vec<double,dim>,2>& seg);

template <size_t dim> 
std::vector<Vec<double,dim>> facet_facet_intersection(
    const Vec<Vec<double,dim>,dim>& fA, const Vec<Vec<double,dim>,dim>& fB);



template <size_t dim>
struct NearestPoint {
    const Vec<double,dim-1> ref_pt;
    const Vec<double,dim> pt;
    const double distance;
};
template <size_t dim>
NearestPoint<dim> closest_pt_facet(const Vec<double,dim>& pt,
    const Vec<Vec<double,dim>,dim> tri);


template <size_t dim> struct Box;
template <size_t dim> struct Ball;

bool is_point_in_polygon(const Vec<double,2>& pt,
    const std::vector<Vec<double,2>>& poly);

template <size_t dim>
bool is_intersection_box_ball(const Box<dim>& box, const Ball<dim>& ball);

template <size_t dim>
bool is_point_in_ball(const Vec<double,dim>& pt, const Ball<dim>& ball);

bool is_point_in_triangle(const Vec<double,2>& pt, const Vec<Vec<double,2>,3>& tri);

}

#endif
