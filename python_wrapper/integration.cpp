#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "integral_term.h"

template <size_t dim>
void export_integration() {
    using namespace boost::python;
    using namespace tbem;
    class_<IntegrationStrategy<dim,1,1>>
        ("IntegrationStrategyScalar");
    class_<IntegrationStrategy<dim,dim,dim>>
        ("IntegrationStrategyTensor");

    def("make_adaptive_integrator", make_adaptive_integrator<dim,1,1>);
    def("make_adaptive_integrator", make_adaptive_integrator<dim,dim,dim>);

    def("make_sinh_integrator", make_sinh_integrator<dim,1,1>);
    def("make_sinh_integrator", make_sinh_integrator<dim,dim,dim>);
}
template void export_integration<2>();
template void export_integration<3>();
