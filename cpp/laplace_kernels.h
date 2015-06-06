#ifndef TBEMddddddddddddEQwQWEA_LAPLACE_KERNELS_H
#define TBEMddddddddddddEQwQWEA_LAPLACE_KERNELS_H

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
        (void)delta; (void)nobs; (void)nsrc;
        return {{{1.0 / (4.0 * M_PI * std::sqrt(r2))}}};
    }

    virtual std::unique_ptr<Kernel<3,1,1>> clone() const
    {
        return std::unique_ptr<Kernel<3,1,1>>(new LaplaceSingle<3>());
    }
};

template <>
struct LaplaceDouble<3>: Kernel<3,1,1> 
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,3>& delta,
        const Vec<double,3>& nobs, const Vec<double,3>& nsrc) const 
    {
        (void)nobs;
        return {{{dot_product(nsrc, delta) / (4 * M_PI * r2 * std::sqrt(r2))}}};
    }

    virtual std::unique_ptr<Kernel<3,1,1>> clone() const
    {
        return std::unique_ptr<Kernel<3,1,1>>(new LaplaceDouble<3>());
    }
};

template <>
struct LaplaceHypersingular<3>: Kernel<3,1,1>
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,3>& delta,
        const Vec<double,3>& nobs, const Vec<double,3>& nsrc) const 
    {
        (void)r2; (void) delta; (void)nobs; (void)nsrc;
        return {{0}};
    }

    virtual std::unique_ptr<Kernel<3,1,1>> clone() const
    {
        return std::unique_ptr<Kernel<3,1,1>>(new LaplaceHypersingular<3>());
    }
};

template <>
struct LaplaceSingle<2>: Kernel<2,1,1>
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,2>& delta,
        const Vec<double,2>& nobs, const Vec<double,2>& nsrc) const 
    {
        (void)delta; (void)nobs; (void)nsrc;
        return {{{std::log(std::sqrt(r2)) / (2 * M_PI)}}};
    }

    virtual std::unique_ptr<Kernel<2,1,1>> clone() const
    {
        return std::unique_ptr<Kernel<2,1,1>>(new LaplaceSingle<2>());
    }
};

template <>
struct LaplaceDouble<2>: Kernel<2,1,1>
{
    virtual Vec1<Vec1<double>> call(double r2, const Vec<double,2>& delta,
        const Vec<double,2>& nobs, const Vec<double,2>& nsrc) const 
    {
        (void)nobs;
        return {{{dot_product(nsrc, delta) / (2 * M_PI * r2)}}};
    }

    virtual std::unique_ptr<Kernel<2,1,1>> clone() const
    {
        return std::unique_ptr<Kernel<2,1,1>>(new LaplaceDouble<2>());
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

    virtual std::unique_ptr<Kernel<2,1,1>> clone() const
    {
        return std::unique_ptr<Kernel<2,1,1>>(new LaplaceHypersingular<2>());
    }
};

} // END namespace tbem

#endif
