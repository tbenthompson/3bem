#ifndef __TEST_SHARED_H
#define __TEST_SHARED_H

#include <vector>
#include <array>
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

inline std::array<std::vector<double>,3> random_pts(int N) {
    std::array<std::vector<double>,3> locs = 
        {random_list(N), random_list(N), random_list(N)};
    return locs;
}

inline std::array<double,3> random_pt() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    return {
        dis(gen), dis(gen), dis(gen)
    };
}

inline std::array<double,3> spherify_pt(std::array<double,3> pt, 
                                        std::array<double,3> c, double r) {
    double rf = pt[0];
    double tf = pt[1] * 2 * M_PI;
    double pf = pt[2] * M_PI;
    std::array<double,3> ret = {
        c[0] + r * rf * std::cos(tf) * std::sin(pf),
        c[1] + r * rf * std::sin(tf) * std::sin(pf),
        c[2] + r * rf * std::cos(pf)
    };
    return ret;
}

inline std::array<double,3> random_pt_sphere(std::array<double,3> c, double r) {
    auto pt = random_pt();
    return spherify_pt(pt, c, r);
}

class Mesh;

//TODO: Move these out to a mesh_gen module
Mesh square_mesh();
Mesh refined_square_mesh(int levels);
Mesh circle_mesh(std::array<double, 2> center, double r, int n_segments);

//TODO: Move these to a separate kernels module
inline double BEMone(double r2, std::array<double,3> delta,
                         std::array<double,3> n) {
    return 1.0;
}

inline double BEMlaplace_single(double r2,
                             std::array<double,3> delta,
                             std::array<double,3> nsrc) {
    return 1.0 / (4.0 * M_PI * std::sqrt(r2));
}

inline double BEMlaplace_double(double r2,
                             std::array<double,3> delta,
                             std::array<double,3> nsrc) {
    return -(nsrc[0] * delta[0] + nsrc[1] * delta[1] + nsrc[2] * delta[2]) / 
           (4 * M_PI * pow(r2, 1.5));
}

#endif
