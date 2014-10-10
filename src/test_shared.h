#ifndef __TEST_SHARED_H
#define __TEST_SHARED_H

#include <memory>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>

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


inline std::vector<std::array<double,3>> three_pts() {
    std::vector<std::array<double,3>> es(3);
    es[0] = {1.0, 2.0, 0.0};
    es[1] = {-1.0, 0.0, -3.0};
    es[2] = {0.0, -2.0, 3.0};
    return es;
}

inline std::vector<double> random_list(int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<double> es(N);
    for (int i = 0; i < N; ++i) {
        es[i] = dis(gen);
    }
    return es;
}

inline std::array<std::vector<double>,3> random_pts(int N) {
    std::array<std::vector<double>,3> locs = 
        {random_list(N), random_list(N), random_list(N)};
    return locs;
}

class Mesh;

//TODO: Move these out to a mesh_gen module
Mesh square_mesh();
Mesh refined_square_mesh(int levels);
Mesh circle_mesh(std::array<double, 2> center, double r, int n_segments);

#endif
