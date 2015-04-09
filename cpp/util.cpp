#include "util.h"
#include <random>

namespace tbem {

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

template <size_t dim>
std::vector<Vec<double,dim>> random_pts(size_t N) 
{
    std::vector<Vec<double,dim>> out(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    for (size_t i = 0; i < N; i++) {
        for (size_t d = 0; d < dim; d++) {
            out[i][d] = dis(gen);
        }
    }
    return out;
}
template std::vector<Vec<double,2>> random_pts<2>(size_t N); 
template std::vector<Vec<double,3>> random_pts<3>(size_t N);

template <size_t dim>
Vec<double,dim> random_pt() 
{
    return random_pts<dim>(1)[0];
}
template Vec<double,2> random_pt<2>();
template Vec<double,3> random_pt<3>();

} //END NAMESPACE TBEM
