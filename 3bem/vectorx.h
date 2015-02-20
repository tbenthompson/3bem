#ifndef __rH1JH11LHLHAHAh_VEcTORX_H
#define __rH1JH11LHLHAHAh_VEcTORX_H

#include <vector>
#include <iosfwd>
#include "fwd_vectorx.h"

namespace tbem {

template <typename T>
struct InternalVec {
    typedef T ValueType;
    typedef InternalVec<T> MyType;
    typedef std::vector<T> ContainerType;

    ContainerType storage;

    InternalVec(); 
    explicit InternalVec(size_t n_elements);
    InternalVec(size_t n_elements, const T& value);
    InternalVec(const ContainerType& data);
    InternalVec(std::initializer_list<T> s);

    void resize(size_t new_size);

    ValueType& operator[] (size_t idx);
    const ValueType& operator[] (size_t idx) const;
    typename ContainerType::iterator begin();
    typename ContainerType::const_iterator begin() const;
    typename ContainerType::iterator end();
    typename ContainerType::const_iterator end() const;

    ValueType* data();
    const ValueType* data() const;

    size_t size() const;

    MyType& operator+=(const MyType& b);
    MyType& operator+=(double b);
    MyType& operator-=(const MyType& b);
    MyType& operator-=(double b);
    MyType& operator*=(const MyType& b);
    MyType& operator*=(double b);

    MyType operator-();

    MyType operator+(const MyType& rhs);
    MyType operator+(double rhs);
    MyType operator-(const MyType& rhs);
    MyType operator-(double rhs);
    MyType operator*(const MyType& rhs);
    MyType operator*(double rhs);

    bool operator==(const MyType& b) const;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const InternalVec<T>& a);

} // end namespace tbem

#endif
