#ifndef __rH1JH11LHLHAHAh_FUNCTION_H
#define __rH1JH11LHLHAHAh_FUNCTION_H

#include <vector>
#include <iosfwd>

namespace tbem {

template <typename T>
struct InternalFnc {
    typedef T ValueType;
    typedef InternalFnc<T> MyType;
    typedef std::vector<T> ContainerType;

    ContainerType _data;

    InternalFnc(); 
    explicit InternalFnc(size_t n_elements);
    InternalFnc(size_t n_elements, const T& value);
    InternalFnc(const ContainerType& data);
    InternalFnc(std::initializer_list<T> s);

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

typedef InternalFnc<double> Function;
typedef InternalFnc<Function> BlockFunction;

template <typename T>
std::ostream& operator<<(std::ostream& os, const InternalFnc<T>& a);


} // end namespace tbem

#endif
