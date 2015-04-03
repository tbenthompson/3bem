#ifndef __ALSDJLJLA_DIMENSION_H
#define __ALSDJLJLA_DIMENSION_H
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "mesh.h"
#include "constraint_builder.h"
#include "quadrature.h"
#include "vertex_iterator.h"
#include "obs_pt.h"
#include "dense_builder.h"

#include "identity_kernels.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"

#include "basis.h"

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
}

template <size_t dim>
void export_dimension() {
    using namespace boost::python;
    using namespace tbem;
    class_<VertexIterator<dim>>("VertexIterator", no_init);
    class_<Mesh<dim>>("Mesh")
        .def("get_vertex", &Mesh<dim>::get_vertex,
               return_value_policy<reference_existing_object>())
        .def("refine_repeatedly", &Mesh<dim>::refine_repeatedly)
        .def("begin", &Mesh<dim>::begin)
        .def("n_facets", &Mesh<dim>::n_facets)
        .def("n_dofs", &Mesh<dim>::n_dofs);
    

    class_<QuadStrategy<dim>>("QuadStrategy", init<int,int,int,double,double>());

    class_<OverlapMap<dim>>("OverlapMap")
        .def("size", &OverlapMap<dim>::size);
    def("mesh_continuity", mesh_continuity<dim>);
    def("cut_at_intersection", cut_at_intersection<dim>);

    def("convert_to_constraints", convert_to_constraints<dim>);
    def("interpolate_bc_constraints", interpolate_bc_constraints_wrapper<dim>);

    def("interpolate", interpolate_wrapper<dim>); 

    class_<Kernel<dim,1,1>, boost::noncopyable>("Kernel", no_init);
    class_<Kernel<dim,dim,dim>, boost::noncopyable>("Kernel", no_init);
    class_<IdentityScalar<dim>, bases<Kernel<dim,1,1>>>("IdentityScalar");
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

    class_<BoundaryIntegral<dim,1,1>>("BoundaryIntegralScalar", no_init);
    class_<BoundaryIntegral<dim,dim,dim>>("BoundaryIntegralTensor", no_init);
    def("make_boundary_integral", make_boundary_integral<dim,1,1>);
    def("make_boundary_integral", make_boundary_integral<dim,dim,dim>);

    def("mesh_to_mesh_operator", mesh_to_mesh_operator<dim,1,1>);
    def("mesh_to_mesh_operator", mesh_to_mesh_operator<dim,dim,dim>);
    class_<ObsPt<dim>>("ObsPt", 
        init<double, Vec<double,dim>, Vec<double,dim>, Vec<double,dim>>())
        .def_readonly("loc", &ObsPt<dim>::loc);
    VectorFromIterable().from_python<std::vector<ObsPt<dim>>>();
    def("mesh_to_points_operator", mesh_to_points_operator<dim,1,1>);
    def("mass_operator", mass_operator<dim,1,1>);
}
#endif
