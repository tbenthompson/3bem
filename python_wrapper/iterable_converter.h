#ifndef __lkLADJSLKAJH_ARRAY_FROM_LIST_H
#define __lkLADJSLKAJH_ARRAY_FROM_LIST_H
#include <array>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

struct ArrayFromIterable 
{
    template <typename T, size_t dim>
    ArrayFromIterable& from_python() {
        boost::python::converter::registry::push_back(
                &ArrayFromIterable::is_convertible,
                &ArrayFromIterable::construct<T,dim>, 
                boost::python::type_id<std::array<T,dim>>());
        return *this;
    }

    static void* is_convertible(PyObject* object)
    {
        return PyObject_GetIter(object) ? object : NULL;
    }

    template <typename T, size_t dim>
    static void construct(PyObject* object,
        boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        namespace bp = boost::python;
        // Object is a borrowed reference, so create a handle indicating 
        // it is borrowed for proper reference counting.
        bp::object obj(bp::handle<>(bp::borrowed(object)));

        // Obtain a handle to the memory block that the converter 
        // has allocated for the C++ type.
        typedef bp::converter::rvalue_from_python_storage<std::array<T,dim>> storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        bp::stl_input_iterator<T> it(obj);

        std::array<T,dim>& out = *(new (storage) std::array<T,dim>());
        for (size_t i = 0; i < dim; i++) {
            out[i] = *it;
            it++;
        }

        data->convertible = storage;
    }
};

struct VectorFromIterable 
{
    template <typename Container>
    VectorFromIterable& from_python() {
        boost::python::converter::registry::push_back(
                &VectorFromIterable::is_convertible,
                &VectorFromIterable::construct<Container>, 
                boost::python::type_id<Container>());
        return *this;
    }

    static void* is_convertible(PyObject* object)
    {
        auto convertible = PyObject_GetIter(object) ? object : NULL;
        return convertible;
    }

    template <typename Container>
    static void construct(PyObject* object,
        boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        namespace bp = boost::python;

        // Object is a borrowed reference, so create a handle indicating 
        // it is borrowed for proper reference counting.
        bp::handle<> handle(bp::borrowed(object));

        // Obtain a handle to the memory block that the converter 
        // has allocated for the C++ type.
        typedef bp::converter::rvalue_from_python_storage<Container> storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        typedef bp::stl_input_iterator<typename Container::value_type> iterator;

        // Allocate the C++ type into the converter's memory block, and 
        // assign its handle to the converter's convertible variable.  The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.
        data->convertible = new (storage) Container(
            iterator(bp::object(handle)), // begin
            iterator());                      // end
    }
};

#endif
