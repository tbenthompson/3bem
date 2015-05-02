#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>
#include "iterable_converter.h"

#include "mesh_gen.h"
#include "mesh.h"
namespace p = boost::python;
namespace np = boost::numpy;

namespace tbem {

template <size_t dim>
struct TensorsToNPArray
{
    static PyObject* convert(const std::vector<Vec<Vec<double,dim>,dim>>& tensors)
    {
        size_t bytes = sizeof(double);
        auto out = np::from_data(
            reinterpret_cast<const double*>(tensors.data()),
            np::dtype::get_builtin<double>(),
            p::make_tuple(tensors.size(), dim, dim),
            p::make_tuple(dim * dim * bytes, dim * bytes, bytes),
            p::object()
        );

        return p::incref(out.copy().ptr());
    }
};

template <typename T>
struct VectorToNPArray
{
    static PyObject* convert(const std::vector<T>& array) 
    {
        size_t bytes = sizeof(T);
        auto out = np::from_data(
            reinterpret_cast<const T*>(array.data()),
            np::dtype::get_builtin<T>(),
            p::make_tuple(array.size()),
            p::make_tuple(bytes),
            p::object()
        );

        return p::incref(out.copy().ptr());
    }
};

template <typename T, size_t dim>
struct ArrayToNPArray
{
    static PyObject* convert(const std::array<T,dim>& array) 
    {
        size_t bytes = sizeof(T);
        auto out = np::from_data(
            reinterpret_cast<const T*>(array.data()),
            np::dtype::get_builtin<T>(),
            p::make_tuple(dim),
            p::make_tuple(bytes),
            p::object()
        );

        return p::incref(out.copy().ptr());
    }
};

template <typename T, size_t dim>
struct VectorOfArraysToNPArray
{
    static PyObject* convert(const std::vector<std::array<T,dim>>& vec) 
    {
        size_t bytes = sizeof(T);
        auto out = np::from_data(
            reinterpret_cast<const T*>(vec.data()),
            np::dtype::get_builtin<T>(),
            p::make_tuple(vec.size(), dim),
            p::make_tuple(dim * bytes, bytes),
            p::object()
        );

        return p::incref(out.copy().ptr());
    }
};

} //end namespace tbem

void export_util() {
    using namespace boost::python;
    using namespace tbem;

    to_python_converter<std::vector<std::array<std::array<double,2>,2>>, 
                        TensorsToNPArray<2>>();
    to_python_converter<std::vector<std::array<std::array<double,3>,3>>, 
                        TensorsToNPArray<3>>();
    to_python_converter<std::vector<double>, VectorToNPArray<double>>();
    to_python_converter<std::vector<int>, VectorToNPArray<int>>();
    to_python_converter<std::vector<size_t>, VectorToNPArray<size_t>>();
    to_python_converter<std::array<double,2>, ArrayToNPArray<double,2>>();
    to_python_converter<std::array<double,3>, ArrayToNPArray<double,3>>();
    to_python_converter<std::vector<std::array<double,2>>,
                        VectorOfArraysToNPArray<double,2>>();
    to_python_converter<std::vector<std::array<double,3>>,
                        VectorOfArraysToNPArray<double,3>>();

    ArrayFromIterable().from_python<double,2>();
    ArrayFromIterable().from_python<double,3>();

    VectorFromIterable().from_python<std::vector<int>>();
    VectorFromIterable().from_python<std::vector<double>>();
    VectorFromIterable().from_python<std::vector<size_t>>();

    def("line_mesh", line_mesh);
    def("circle_mesh", circle_mesh);
    def("sphere_mesh", sphere_mesh);
    def("rect_mesh", rect_mesh);
}
