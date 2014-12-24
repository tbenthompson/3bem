#include "util.h"

namespace tbem {

/* equivalent to range(min, max) in python */
std::vector<int> naturals(int min, int max) {
    std::vector<int> indices(max);
    std::iota(indices.begin(), indices.end(), min);
    return indices;
}

/* equivalent to range(0, max) in python */
std::vector<int> naturals(int max) {
    return naturals(0, max);
}

std::array<std::vector<double>,3> three_pts() {
    std::array<std::vector<double>,3> es;
    es[0] = {1.0, -1.0, 0.0};
    es[1] = {2.0, 0.0, -2.0};
    es[2] = {0.0, -3.0, 3.0};
    return es;
}

std::vector<double> random_list(int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<double> es(N);
    for (int i = 0; i < N; ++i) {
        es[i] = dis(gen);
    }
    return es;
}

std::array<std::vector<double>,3> random_pts(int N) {
    std::array<std::vector<double>,3> locs = 
        {random_list(N), random_list(N), random_list(N)};
    return locs;
}

double random_val() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
}

Vec3<double> random_pt() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    return {
        dis(gen), dis(gen), dis(gen)
    };
}

Vec3<double> spherify_pt(Vec3<double> pt, Vec3<double> c, double r) {
    double rf = pt[0];
    double tf = pt[1] * 2 * M_PI;
    double pf = pt[2] * M_PI;
    Vec3<double> ret = {
        c[0] + r * rf * std::cos(tf) * std::sin(pf),
        c[1] + r * rf * std::sin(tf) * std::sin(pf),
        c[2] + r * rf * std::cos(pf)
    };
    return ret;
}

template <>
Vec3<double> random_pt_sphere<3>(Vec3<double> c, double r) {
    auto pt = random_pt();
    return spherify_pt(pt, c, r);
}

template <>
Vec2<double> random_pt_sphere<2>(Vec2<double> c, double r) {
    auto sphere_pt = random_pt_sphere<3>({c[0], c[1], 0.0}, r);
    return {sphere_pt[0], sphere_pt[1]};
}

double error_inf(const std::vector<double>& a, 
                        const std::vector<double>& b) {
    double error = 0.0;
    for (unsigned int i = 0; i < a.size(); i++) {
        error = std::max(std::fabs(a[i] - b[i]), error);
    }
    return error;
}

} //END NAMESPACE TBEM
