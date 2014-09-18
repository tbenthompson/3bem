#ifndef __TEST_SHARED_H
#define __TEST_SHARED_H

#include <memory>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>

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


inline std::array<std::vector<double>,3> three_pts() {
    std::array<std::vector<double>,3> es;
    es[0] = {1.0, -1.0, 0.0};
    es[1] = {2.0, 0.0, -2.0};
    es[2] = {0.0, -3.0, 3.0};
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

#endif
