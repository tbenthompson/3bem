#include <boost/python.hpp>
#include "iterable_converter.h"
#include "vertex_iterator.h"
#include "mesh.h"
#include "mesh_preprocess.h"
#include "mesh_converter.h"
#include "interior_mesh.h"
namespace p = boost::python;

template <size_t dim>
void export_mesh() {
    using namespace tbem;
    p::class_<VertexIterator<dim>>("VertexIterator", p::no_init);
    p::class_<Mesh<dim>>("Mesh", p::init<std::vector<Facet<dim>>>())
        .add_property("facets", 
            make_getter(&Mesh<dim>::facets,
                p::return_value_policy<p::return_by_value>()))
        .def("get_vertex", &Mesh<dim>::get_vertex,
               p::return_value_policy<p::reference_existing_object>())
        .def("refine_repeatedly", &Mesh<dim>::refine_repeatedly)
        .def("refine", &Mesh<dim>::refine)
        .def("remove_facets", &Mesh<dim>::remove_facets)
        .def("begin", &Mesh<dim>::begin)
        .def("n_facets", &Mesh<dim>::n_facets)
        .def("n_dofs", &Mesh<dim>::n_dofs)
        .def("create_union", &Mesh<dim>::create_union).staticmethod("create_union");
    VectorFromIterable().from_python<std::vector<Mesh<dim>>>();

    p::class_<std::vector<FacetIntersection<dim>>>(
        "ArrayOfFacetIntersection", p::no_init);
    p::class_<MeshPreprocessor<dim>>("MeshPreprocessor")
        .def("find_intersections", &MeshPreprocessor<dim>::find_intersections)
        .def("split_facets_at_intersections",
            &MeshPreprocessor<dim>::split_facets_at_intersections);

    p::class_<PtIndexMesh<dim>>("PtIndexMesh")
        .add_property("points", 
            make_getter(&PtIndexMesh<dim>::points,
                p::return_value_policy<p::return_by_value>()))
        .add_property("facets", 
            make_getter(&PtIndexMesh<dim>::facets,
                p::return_value_policy<p::return_by_value>()));
    p::def("convert_facet_to_pt_index", convert_facet_to_pt_index<dim>);
    p::def("identify_regions", identify_regions);
}

template void export_mesh<2>();
template void export_mesh<3>(); 
