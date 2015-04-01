#ifndef _LHJASDLHJJAH123__ITERABLE_CONVERTER_H
#define _LHJASDLHJJAH123__ITERABLE_CONVERTER_H

/// I (Ben) took this from 
// http://stackoverflow.com/questions/15842126/feeding-a-python-list-into-a-function-taking-in-a-vector-with-boost-python

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

namespace bp = boost::python;

// Type that allows for registration of conversions from python iterable types.
struct iterable_converter
{
    // Registers converter from a python interable type to the provided type.
    template <typename Container> 
    iterable_converter& from_python()
    {
        bp::converter::registry::push_back(
                &iterable_converter::is_convertible,
                &iterable_converter::construct<Container>, 
                bp::type_id<Container>());
        return *this;
    }

    static void* is_convertible(PyObject* object)
    {
        return PyObject_GetIter(object) ? object : NULL;
    }

    // Convert iterable PyObject to C++ container type.
    //
    // Container Concept requirements:
    //   * Container::value_type is CopyConstructable.
    //   * Container can be constructed and populated with two iterators.
    //     I.e. Container(begin, end)
    template <typename Container>
    static void construct(PyObject* object,
        bp::converter::rvalue_from_python_stage1_data* data)
    {
        namespace python = bp;
        // Object is a borrowed reference, so create a handle indicting 
        // it is borrowed for proper reference counting.
        bp::handle<> handle(bp::borrowed(object));

        // Obtain a handle to the memory block that the converter 
        // has allocated for the C++ type.
        typedef bp::converter::rvalue_from_python_storage<Container>
                                                     storage_type;
        void* storage = 
            reinterpret_cast<storage_type*>(data)->storage.bytes;

        typedef bp::stl_input_iterator<typename Container::value_type>
                                                         iterator;

        // Allocate the C++ type into the converter's memory block, and 
        // assign its handle to the converter's convertible variable.  The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.

        // TODO: POSSIBLE MEMORY LEAK?
        data->convertible = new (storage) Container(
            iterator(bp::object(handle)),
            iterator()
        );
    }
};
#endif
