#ifndef __TEST_SHARED_H
#define __TEST_SHARED_H

#include <memory>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include "geom.h"

#define TIC\
    std::chrono::high_resolution_clock::time_point start =\
        std::chrono::high_resolution_clock::now();
#define TIC2\
    start = std::chrono::high_resolution_clock::now();
#define TOC(name)\
    std::cout << name << " took "\
              << std::chrono::duration_cast<std::chrono::milliseconds>\
                (std::chrono::high_resolution_clock::now() - start).count()\
              << "ms.\n";


std::unique_ptr<std::vector<Vec<3> > > three_pts();

template<int dim>
std::unique_ptr<std::vector<Vec<dim> > > random_pts(int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    std::unique_ptr<std::vector<Vec<dim> > > es(new std::vector<Vec<dim> >);;
    for (int i = 0; i < N; ++i) {
        Vec<dim> v;
        for (int d = 0; d < dim; ++d) {
            v.loc[d] = dis(gen); 
        }
        es->push_back(v);
    }
    return es;
}


#endif
