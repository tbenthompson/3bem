#ifndef __TEST_SHARED_H
#define __TEST_SHARED_H

#include <memory>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include "bem.h"

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

inline Mesh square_mesh() {
    std::vector<std::array<double, 2>> vertices = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0},
    };
    
    std::vector<std::array<int, 2>> segs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };

    Mesh m = {vertices, segs};
    return m;
}

inline Mesh refined_square_mesh(int levels) {

    Mesh m = square_mesh();
    for (int i = 0; i < levels; i++) {
        m = refine_mesh(m, naturals(m.segments.size()));
    }
    return m;
}

inline Mesh circle_mesh(std::array<double, 2> center, double r, int n_segments) {
    std::vector<std::array<double, 2>> vertices(n_segments);
    std::vector<std::array<int, 2>> segs(n_segments);

    for (int i = 0; i < n_segments; i++) {
        double theta = (2 * M_PI) * (i / (double)n_segments);
        vertices[i] = {center[0] + r * cos(theta), center[1] + r * sin(theta)};
        segs[i] = {i, i + 1};
    }
    segs[n_segments - 1][1] = 0;

    Mesh m = {vertices, segs};
    return m;
}

#endif
