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

template <typename T, size_t dim>
struct NPArrayToVectorOfArrays {
    typedef std::vector<std::array<T,dim>> OutType;

    NPArrayToVectorOfArrays() {
        p::converter::registry::push_back(
            &convertible,
            &construct,
            p::type_id<OutType>()
        );
    }

    static void* convertible(PyObject* p) {
        try {
            p::object obj(p::handle<>(p::borrowed(p)));
            std::unique_ptr<np::ndarray> array(
                new np::ndarray(
                    np::from_object(
                        obj, np::dtype::get_builtin<T>(), 2,
                        2, np::ndarray::V_CONTIGUOUS
                    )
                )
            );
            if (array->shape(1) != dim) {
                return 0;
            }
            return array.release();
        } catch (p::error_already_set & err) {
            p::handle_exception();
            return 0;
        }
    }

    static void construct(PyObject* obj,
        p::converter::rvalue_from_python_stage1_data* data) 
    {
        // Extract the array we passed out of the convertible() member function.
        std::unique_ptr<np::ndarray> array(
            reinterpret_cast<np::ndarray*>(data->convertible)
        );
        // Find the memory block Boost.Python has prepared for the result.
        typedef p::converter::rvalue_from_python_storage<T> storage_t;
        storage_t * storage = reinterpret_cast<storage_t*>(data);
        // Use placement new to initialize the result.
        OutType* out = new (storage->storage.bytes) OutType(array->shape(0));

        // Fill the result with the values from the NumPy array.
        auto strides = array->get_strides();
        for (size_t i = 0; i < array->shape(0); i++) {
            for (size_t d1 = 0; d1 < dim; d1++) {
                auto offset = i * strides[0] + d1 * strides[1];
                auto byte_ptr = array->get_data() + offset;
                (*out)[i][d1] = *reinterpret_cast<T*>(byte_ptr);
            }
        }
        // Finish up.
        data->convertible = storage->storage.bytes;
    }
};

template <typename T, size_t dim>
struct NPArrayToVectorOfTensors {
    typedef std::vector<std::array<std::array<T,dim>,dim>> OutType;

    NPArrayToVectorOfTensors() {
        p::converter::registry::push_back(
            &convertible,
            &construct,
            p::type_id<OutType>()
        );
    }

    static void* convertible(PyObject* p) {
        try {
            p::object obj(p::handle<>(p::borrowed(p)));
            std::unique_ptr<np::ndarray> array(
                new np::ndarray(
                    np::from_object(
                        obj, np::dtype::get_builtin<T>(), 3,
                        3, np::ndarray::V_CONTIGUOUS
                    )
                )
            );
            if (array->shape(1) != dim) {
                return 0;
            }
            if (array->shape(2) != dim) {
                return 0;
            }
            return array.release();
        } catch (p::error_already_set & err) {
            p::handle_exception();
            return 0;
        }
    }

    static void construct(PyObject* obj,
        p::converter::rvalue_from_python_stage1_data* data) 
    {
        // Extract the array we passed out of the convertible() member function.
        std::unique_ptr<np::ndarray> array(
            reinterpret_cast<np::ndarray*>(data->convertible)
        );
        // Find the memory block Boost.Python has prepared for the result.
        typedef p::converter::rvalue_from_python_storage<T> storage_t;
        storage_t * storage = reinterpret_cast<storage_t*>(data);
        // Use placement new to initialize the result.
        OutType* out = new (storage->storage.bytes) OutType(array->shape(0));

        // Fill the result with the values from the NumPy array.
        auto strides = array->get_strides();
        for (size_t i = 0; i < array->shape(0); i++) {
            for (size_t d1 = 0; d1 < dim; d1++) {
                for (size_t d2 = 0; d2 < dim; d2++) {
                    auto offset = i * strides[0] + d1 * strides[1] + d2 * strides[2];
                    auto byte_ptr = array->get_data() + offset;
                    (*out)[i][d1][d2] = *reinterpret_cast<T*>(byte_ptr);
                }
            }
        }
        // Finish up.
        data->convertible = storage->storage.bytes;
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
    to_python_converter<std::vector<std::array<size_t,2>>,
                        VectorOfArraysToNPArray<size_t,2>>();
    to_python_converter<std::vector<std::array<size_t,3>>,
                        VectorOfArraysToNPArray<size_t,3>>();

    NPArrayToVectorOfArrays<double,2>();
    NPArrayToVectorOfArrays<double,3>();
    NPArrayToVectorOfArrays<size_t,2>();
    NPArrayToVectorOfArrays<size_t,3>();
    NPArrayToVectorOfTensors<double,2>();
    NPArrayToVectorOfTensors<double,3>();

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
