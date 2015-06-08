#ifndef TBEMTEST_SHARED_H
#define TBEMTEST_SHARED_H

#include <vector>
#include <chrono>
#include "vec.h"

namespace tbem {

#define TIC\
    std::chrono::high_resolution_clock::time_point start =\
        std::chrono::high_resolution_clock::now();\
    int time_ms;
#define TIC2\
    start = std::chrono::high_resolution_clock::now();
#define TOC(name)\
    time_ms = std::chrono::duration_cast<std::chrono::milliseconds>\
                (std::chrono::high_resolution_clock::now() - start).count();\
    std::cout << name << " took "\
              << time_ms\
              << "ms.\n";

std::vector<double> random_list(size_t N, double a = 0, double b = 1);

template <size_t dim>
Vec<double,dim> random_pt();

template <typename T>
T random(T min, T max);

template <size_t dim>
std::vector<Vec<double,dim>> random_pts(size_t N, double a = 0, double b = 1);

//TODO: Should this be here? Or in geometry.h
template <size_t dim> struct Ball;
template <size_t dim>
std::vector<Ball<dim>> random_balls(size_t n, double r_max);

} //END namespace tbem
#endif
