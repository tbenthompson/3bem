#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "integral_term.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
struct IntegrationMethodIWrap: public IntegrationMethodI<dim,R,C>,
    boost::python::wrapper<IntegrationMethodI<dim,R,C>> 
{
    virtual Vec<Vec<Vec<double,C>,R>,dim>
    compute_singular(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const {
        return this->get_override("compute_singular")();
    }

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_nearfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const {
        return this->get_override("compute_nearfield")();
    }

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_farfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const {
        return this->get_override("compute_farfield")();
    }

    virtual double far_threshold() const {
        return this->get_override("far_threshold")();
    }
    virtual QuadRule<dim-1> get_obs_quad() const {
        return this->get_override("get_obs_quad")();
    }
};

} // end namepsace tbem

template <size_t dim>
void export_integration() {
    using namespace boost::python;
    using namespace tbem;
    class_<IntegrationMethodIWrap<dim,1,1>, boost::noncopyable>
        ("IntegrationMethodIScalar");
    class_<IntegrationMethodIWrap<dim,dim,dim>, boost::noncopyable>
        ("IntegrationMethodITensor");

    class_<AdaptiveIntegrationMethod<dim,1,1>, bases<IntegrationMethodI<dim,1,1>>>
        ("AdaptiveIntegrationMethodScalar", no_init);
    class_<AdaptiveIntegrationMethod<dim,dim,dim>, bases<IntegrationMethodI<dim,dim,dim>>>
        ("AdaptiveIntegrationMethodTensor", no_init);
    def("make_adaptive_integration_mthd", make_adaptive_integration_mthd<dim,1,1>);
    def("make_adaptive_integration_mthd", make_adaptive_integration_mthd<dim,dim,dim>);

    class_<SinhIntegrationMethod<dim,1,1>, bases<IntegrationMethodI<dim,1,1>>>
        ("SinhIntegrationMethodScalar", no_init);
    class_<SinhIntegrationMethod<dim,dim,dim>, bases<IntegrationMethodI<dim,dim,dim>>>
        ("SinhIntegrationMethodTensor", no_init);
    def("make_sinh_integration_mthd", make_sinh_integration_mthd<dim,1,1>);
    def("make_sinh_integration_mthd", make_sinh_integration_mthd<dim,dim,dim>);
}
template void export_integration<2>();
template void export_integration<3>();
