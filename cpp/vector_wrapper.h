#ifndef __12312312312_VECTOR_WRAPPER_H
#define __12312312312_VECTOR_WRAPPER_H

#include <vector>

namespace tbem {

template <typename T>
struct VectorWrapper: public std::vector<T>
{
    typedef typename std::_Vector_base<T,std::allocator<T>>::_Vector_impl* BasePtr;
    VectorWrapper(T* source, size_t source_size) 
    {
        BasePtr ptr = BasePtr((void*)this);
        ptr->_M_start = source;
        ptr->_M_finish = ptr->_M_start + source_size;
        ptr->_M_end_of_storage = ptr->_M_finish;
    }

    ~VectorWrapper() 
    {
        BasePtr ptr = BasePtr((void*)this);
        ptr->_M_start = nullptr;
        ptr->_M_finish = nullptr; 
        ptr->_M_end_of_storage = nullptr;
    }
};

template <typename T>
VectorWrapper<T> make_vector_wrapper(T* source, size_t source_size) 
{
    return VectorWrapper<T>(source, source_size);
}

} // end namespace tbem

#endif
