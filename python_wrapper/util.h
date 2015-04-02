#ifndef __QQQQQQQQQQQQQQQ_UTIL_H
#define __QQQQQQQQQQQQQQQ_UTIL_H
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>
#include "iterable_converter.h"

#include "mesh_gen.h"

template <typename T, size_t dim> 
T get_item_from_std_array(const std::array<T,dim>& A, const boost::python::object& o) 
{
    using namespace boost::python;
    size_t idx = extract<size_t>(o); 
    assert(idx < dim);
    return A[idx];
}

void export_util() {
    using namespace boost::python;
    using namespace tbem;

    ArrayFromIterable().from_python<double,2>();
    ArrayFromIterable().from_python<double,3>();

    class_<std::array<double,2>>("ArrayOf2Doubles")
        .def("__getitem__", get_item_from_std_array<double,2>);
    class_<std::array<double,3>>("ArrayOf3Doubles")
        .def("__getitem__", get_item_from_std_array<double,3>);

    VectorFromIterable().from_python<std::vector<double>>();
    class_<std::vector<double>>("VectorOfDoubles")
        .def(vector_indexing_suite<std::vector<double>>());
    class_<std::vector<Vec<double,2>>>("VectorOf2Doubles")
        .def(vector_indexing_suite<std::vector<Vec<double,2>>>());
    class_<std::vector<Vec<double,3>>>("VectorOf3Doubles")
        .def(vector_indexing_suite<std::vector<Vec<double,3>>>());

    def("line_mesh", line_mesh);
    def("circle_mesh", circle_mesh);
    def("sphere_mesh", sphere_mesh);
    def("rect_mesh", rect_mesh);
}

#endif
