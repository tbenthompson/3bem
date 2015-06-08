#include "util.h"
#include <random>
#include "geometry.h"

namespace tbem {

std::vector<double> random_list(size_t N, double a, double b) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(a, b);
    std::vector<double> es(N);
    for (size_t i = 0; i < N; ++i) {
        es[i] = dis(gen);
    }
    return es;
}

template <size_t dim>
std::vector<Vec<double,dim>> random_pts(size_t N, double a, double b) 
{
    std::vector<Vec<double,dim>> out(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(a, b);
    for (size_t i = 0; i < N; i++) {
        for (size_t d = 0; d < dim; d++) {
            out[i][d] = dis(gen);
        }
    }
    return out;
}

template <>
double random<double>(double min, double max) 
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

template <>
size_t random<size_t>(size_t min, size_t max) 
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dis(min, max);
    return dis(gen);
}

template std::vector<Vec<double,2>> random_pts<2>(size_t N, double a, double b); 
template std::vector<Vec<double,3>> random_pts<3>(size_t N, double a, double b);

template <size_t dim>
Vec<double,dim> random_pt() 
{
    return random_pts<dim>(1)[0];
}
template Vec<double,2> random_pt<2>();
template Vec<double,3> random_pt<3>();

template <size_t dim>
std::vector<Ball<dim>> random_balls(size_t n, double r_max)
{
    auto pts = random_pts<dim>(n);
    auto r = random_list(n, 0, r_max);
    std::vector<Ball<dim>> balls(n);
    for (size_t i = 0; i < n; i++) {
        balls[i] = {pts[i], r[i]};
    }
    return balls;
}

template std::vector<Ball<2>> random_balls(size_t n, double r_max);
template std::vector<Ball<3>> random_balls(size_t n, double r_max);

} //END NAMESPACE TBEM
