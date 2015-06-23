#ifndef TBEMQUAD_RULE_H
#define TBEMQUAD_RULE_H
#include <array>
#include <vector>

namespace tbem {

template <size_t dim>
struct QuadPt {
    const std::array<double,dim> x_hat;
    const double w;
};

template <size_t dim>
using QuadRule = std::vector<QuadPt<dim>>;

} // end namespace tbem

#endif
