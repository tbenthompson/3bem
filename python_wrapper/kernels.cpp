#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "elastic_kernels.h"
#include "laplace_kernels.h"
#include "gravity_kernels.h"

template <size_t dim>
void export_kernels() {
    using namespace boost::python;
    using namespace tbem;
    class_<Kernel<dim,1,1>, boost::noncopyable>("Kernel", no_init);
    class_<Kernel<dim,dim,dim>, boost::noncopyable>("Kernel", no_init);
    class_<LaplaceSingle<dim>, bases<Kernel<dim,1,1>>>("LaplaceSingle");
    class_<LaplaceDouble<dim>, bases<Kernel<dim,1,1>>>("LaplaceDouble");
    class_<LaplaceHypersingular<dim>, bases<Kernel<dim,1,1>>>("LaplaceHypersingular");
    class_<ElasticDisplacement<dim>, bases<Kernel<dim,dim,dim>>>(
        "ElasticDisplacement", init<double,double>());
    class_<ElasticTraction<dim>, bases<Kernel<dim,dim,dim>>>(
        "ElasticTraction", init<double,double>());
    class_<ElasticAdjointTraction<dim>, bases<Kernel<dim,dim,dim>>>(
        "ElasticAdjointTraction", init<double,double>());
    class_<ElasticHypersingular<dim>, bases<Kernel<dim,dim,dim>>>(
        "ElasticHypersingular", init<double,double>());
    class_<GravityDisplacement<dim>, bases<Kernel<dim,dim,dim>>>(
        "GravityDisplacement", init<double,double,Vec<double,dim>>());
    class_<GravityTraction<dim>, bases<Kernel<dim,dim,dim>>>(
        "GravityTraction", init<double,double,Vec<double,dim>>());
}
template void export_kernels<2>();
template void export_kernels<3>();
