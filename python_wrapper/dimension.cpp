#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "iterable_converter.h"

#include "mesh.h"
#include "continuity_builder.h"
#include "quadrature.h"
#include "vertex_iterator.h"
#include "obs_pt.h"
#include "dense_builder.h"
#include "integral_operator.h"
#include "mass_operator.h"
#include "basis.h"
namespace p = boost::python;
namespace np = boost::numpy;

namespace tbem {

template <size_t dim>
VectorX interpolate_wrapper(const Mesh<dim>& mesh, const boost::python::object& fnc) 
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
    const std::vector<int>& which_dofs,
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
void export_integration();
template <size_t dim>
void export_kernels();

template <size_t dim>
void export_dimension() {
    using namespace boost::python;
    using namespace tbem;

    class_<VertexIterator<dim>>("VertexIterator", no_init);
    class_<Mesh<dim>>("Mesh", init<std::vector<Facet<dim>>>())
        .add_property("facets", 
            make_getter(&Mesh<dim>::facets, return_value_policy<return_by_value>()))
        .def("get_vertex", &Mesh<dim>::get_vertex,
               return_value_policy<reference_existing_object>())
        .def("refine_repeatedly", &Mesh<dim>::refine_repeatedly)
        .def("begin", &Mesh<dim>::begin)
        .def("n_facets", &Mesh<dim>::n_facets)
        .def("n_dofs", &Mesh<dim>::n_dofs)
        .def("create_union", &Mesh<dim>::create_union).staticmethod("create_union");
    VectorFromIterable().from_python<std::vector<Mesh<dim>>>();
    

    class_<QuadStrategy<dim>>("QuadStrategy", init<int,int,double,double>());

    class_<OverlapMap<dim>>("OverlapMap")
        .def("size", &OverlapMap<dim>::size);
    def("mesh_continuity", mesh_continuity<dim>);
    def("cut_at_intersection", cut_at_intersection<dim>);

    def("convert_to_constraints", convert_to_constraints<dim>);
    def("interpolate_bc_constraints", interpolate_bc_constraints_wrapper<dim>);

    def("interpolate", interpolate_wrapper<dim>); 

    export_kernels<dim>();
    export_integration<dim>();

    class_<BlockIntegralOperator<dim,1,1>, bases<BlockOperatorI>>
        ("BlockIntegralOperatorScalar", no_init);
    class_<BlockIntegralOperator<dim,dim,dim>, bases<BlockOperatorI>>
        ("BlockIntegralOperatorTensor", no_init);

    def("integral_operator", integral_operator<dim,1,1>);
    def("integral_operator", integral_operator<dim,dim,dim>);
    def("dense_integral_operator", dense_integral_operator<dim,1,1>);
    def("dense_integral_operator", dense_integral_operator<dim,dim,dim>);

    //TODO: don't expose this...??
    //TODO: Instead have a wrapper for mesh_to_points that takes a
    //numpy array for each of the constructor inputs to ObsPt
    class_<ObsPt<dim>>("ObsPt", 
        init<double, Vec<double,dim>, Vec<double,dim>, Vec<double,dim>>())
        .add_property("loc", 
            make_getter(&ObsPt<dim>::loc, return_value_policy<return_by_value>()));
    VectorFromIterable().from_python<std::vector<ObsPt<dim>>>();
    def("mesh_to_points_operator", mesh_to_points_operator<dim,1,1>);


    class_<BlockMassOperator<dim>, bases<BlockOperatorI>>("BlockMassOperator", no_init);
    def("mass_operator_scalar", mass_operator<dim,1,1>);
    def("mass_operator_tensor", mass_operator<dim,dim,dim>);
}
template void export_dimension<2>();
template void export_dimension<3>();
