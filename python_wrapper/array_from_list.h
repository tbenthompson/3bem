#ifndef __lkLADJSLKAJH_ARRAY_FROM_LIST_H
#define __lkLADJSLKAJH_ARRAY_FROM_LIST_H
#include <array>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

namespace tbem {

struct ArrayFromList 
{
    template <typename T, size_t dim>
    ArrayFromList& from_python() {
        boost::python::converter::registry::push_back(
                &ArrayFromList::is_convertible,
                &ArrayFromList::construct<T,dim>, 
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

} // end namespace tbem

#endif
