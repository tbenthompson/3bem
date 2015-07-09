#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include "iterable_converter.h"
#include "op_wrap.h"

#include "continuity_builder.h"
#include "integral_operator.h"
#include "mass_operator.h"
#include "basis.h"
namespace p = boost::python;
namespace np = boost::numpy;

namespace tbem {

template <size_t dim>
std::vector<double>
interpolate_wrapper(const Mesh<dim>& mesh, const boost::python::object& fnc) 
{
    return interpolate<dim>(mesh,
        [&](const Vec<double,dim>& x) {
            double res = boost::python::extract<double>(fnc(x));
            return res;
        });
}

template <size_t dim>
std::vector<ConstraintEQ> interpolate_bc_constraints_wrapper(
    const Mesh<dim>& m,
    const std::vector<size_t>& which_dofs,
    const boost::python::object& fnc) 
{
    return interpolate_bc_constraints<dim>(m, which_dofs,
        [&](const Vec<double,dim>& x) {
            double res = boost::python::extract<double>(fnc(x));
            return res;
        });
}


} //end namespace tbem

template <size_t dim>
void export_mesh();
template <size_t dim>
void export_integration();
template <size_t dim>
void export_kernels();

template <size_t dim>
void export_dimension() {
    using namespace boost::python;
    using namespace tbem;

    export_mesh<dim>(); 

    class_<OverlapMap<dim>>("OverlapMap")
        .def("size", &OverlapMap<dim>::size);
    def("mesh_continuity", mesh_continuity<dim>);
    def("cut_at_intersection", cut_at_intersection<dim>);
    def("normal_constraints", normal_constraints<dim>);
    def("form_neighbor_bcs", form_neighbor_bcs<dim>);
    def("bc_constraints", bc_constraints<dim>);

    def("convert_to_constraints", convert_to_constraints<dim>);
    def("interpolate_bc_constraints", interpolate_bc_constraints_wrapper<dim>);

    def("interpolate", interpolate_wrapper<dim>); 

    export_kernels<dim>();
    export_integration<dim>();

    auto integral_op_scalar = class_<
            IntegralOperator<dim,1,1>
        >("IntegralOperatorScalar", no_init)
        .def_readonly("nearfield", &IntegralOperator<dim,1,1>::nearfield);
    export_operator<IntegralOperator<dim,1,1>>(integral_op_scalar);

    auto integral_op_tensor = class_<IntegralOperator<dim,dim,dim>>(
            "IntegralOperatorTensor", no_init)
        .def_readonly("nearfield", &IntegralOperator<dim,dim,dim>::nearfield);
    export_operator<IntegralOperator<dim,dim,dim>>(integral_op_tensor);

    def("boundary_operator", boundary_operator<dim,1,1>);
    def("boundary_operator", boundary_operator<dim,dim,dim>);
    def("dense_boundary_operator", dense_boundary_operator<dim,1,1>);
    def("dense_boundary_operator", dense_boundary_operator<dim,dim,dim>);

    def("dense_interior_operator", dense_interior_operator<dim,1,1>);
    def("dense_interior_operator", dense_interior_operator<dim,dim,dim>);

    // mass_operator cannot be overloaded like boundary_operator, because there
    // is no kernel parameter to decide on the tensor shape
    def("mass_operator_scalar", mass_operator<dim,1,1>);
    def("mass_operator_tensor", mass_operator<dim,dim,dim>);
}
template void export_dimension<2>();
template void export_dimension<3>();
