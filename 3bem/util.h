#ifndef __TEST_SHARED_H
#define __TEST_SHARED_H

#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <string>
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

std::array<std::vector<double>,3> three_pts();

std::vector<double> random_list(int N);
std::array<std::vector<double>,3> random_pts(int N);
double random_val();
Vec3<double> random_pt();
Vec3<double> spherify_pt(Vec3<double> pt, Vec3<double> c, double r);
template <int dim>
Vec<double,dim> random_pt_sphere(Vec<double,dim> c, double r);

double error_inf(const std::vector<double>& a, const std::vector<double>& b);

} //END namespace tbem
#endif
