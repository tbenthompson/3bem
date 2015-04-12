#ifndef __ddddddddddddEQwQWEA_LAPLACE_KERNELS_H
#define __ddddddddddddEQwQWEA_LAPLACE_KERNELS_H

#include <exception>
#include "vec_ops.h"
#include "kernel.h"

namespace tbem {

/* Laplace/Poisson equation kernels. */
template <size_t dim>
struct LaplaceSingle;
template <size_t dim>
struct LaplaceDouble;
template <size_t dim>
struct LaplaceHypersingular;

template <>
struct LaplaceSingle<3>: Kernel<3,1,1> 
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,3>& delta,
        const Vec<double,3>& nobs, const Vec<double,3>& nsrc) const 
    {
        return {{{1.0 / (4.0 * M_PI * std::sqrt(r2))}}};
    }
};

template <>
struct LaplaceDouble<3>: Kernel<3,1,1> 
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,3>& delta,
        const Vec<double,3>& nobs, const Vec<double,3>& nsrc) const 
    {
        return {{{dot_product(nsrc, delta) / (4 * M_PI * r2 * std::sqrt(r2))}}};
    }
};

template <>
struct LaplaceHypersingular<3>: Kernel<3,1,1>
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,3>& delta,
        const Vec<double,3>& nobs, const Vec<double,3>& nsrc) const 
    {
        throw std::exception();
    }
};

template <>
struct LaplaceSingle<2>: Kernel<2,1,1>
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,2>& delta,
        const Vec<double,2>& nobs, const Vec<double,2>& nsrc) const 
    {
        return {{{std::log(std::sqrt(r2)) / (2 * M_PI)}}};
    }
};

template <>
struct LaplaceDouble<2>: Kernel<2,1,1>
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,2>& delta,
        const Vec<double,2>& nobs, const Vec<double,2>& nsrc) const 
    {
        return {{{dot_product(nsrc, delta) / (2 * M_PI * r2)}}};
    }
};

template <>
struct LaplaceHypersingular<2>: Kernel<2,1,1>
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,2>& delta,
        const Vec<double,2>& nobs, const Vec<double,2>& nsrc) const 
    {
        return {{{((-dot_product(nobs, nsrc) / r2) + 
                ((2 * dot_product(nsrc, delta) * dot_product(nobs, delta)) / (r2 * r2)))
            / (2 * M_PI)}}};
    }
};

} // END namespace tbem

#endif
