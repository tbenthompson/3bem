#include "function.h"
#include <cassert>
#include <iostream>

namespace tbem {

template <typename T>
using ValueType = T;
template <typename T>
using MyType = InternalFnc<T>;
template <typename T>
using ContainerType = std::vector<T>;

template <typename T>
InternalFnc<T>::InternalFnc() {}

template <typename T>
InternalFnc<T>::InternalFnc(size_t n_elements):
    _data(n_elements)
{}

template <typename T>
InternalFnc<T>::InternalFnc(size_t n_elements, const T& value):
    _data(n_elements, value)
{}

template <typename T>
InternalFnc<T>::InternalFnc(const ContainerType& data):
    _data(data)
{}

template <typename T>
InternalFnc<T>::InternalFnc(std::initializer_list<T> s):
    _data(s)
{}

template <typename T>
void InternalFnc<T>::resize(size_t new_size) {
    _data.resize(new_size);
}

template <typename T>
ValueType<T>& InternalFnc<T>::operator[] (size_t idx) {
    return _data[idx];
}

template <typename T>
const ValueType<T>& InternalFnc<T>::operator[] (size_t idx) const {
    return _data[idx];
}

template <typename T>
typename ContainerType<T>::iterator InternalFnc<T>::begin() {
    return _data.begin();
}

template <typename T>
typename ContainerType<T>::const_iterator InternalFnc<T>::begin() const {
    return _data.begin();
}

template <typename T>
typename ContainerType<T>::iterator InternalFnc<T>::end() {
    return _data.end();
}

template <typename T>
typename ContainerType<T>::const_iterator InternalFnc<T>::end() const {
    return _data.end();
}

template <typename T>
ValueType<T>* InternalFnc<T>::data() {
    return _data.data();
}

template <typename T>
const ValueType<T>* InternalFnc<T>::data() const {
    return _data.data();
}

template <typename T>
size_t InternalFnc<T>::size() const {
    return _data.size();
}

template <typename T>
MyType<T>& InternalFnc<T>::operator+=(const MyType& b) {
    assert(size() == b.size());
    for (size_t i = 0; i < size(); i++) {
        _data[i] += b[i];
    }
    return *this;
}

template <typename T>
MyType<T>& InternalFnc<T>::operator+=(double b) {
    for (size_t i = 0; i < size(); i++) {
        _data[i] += b;
    }
    return *this;
}

template <typename T>
MyType<T>& InternalFnc<T>::operator-=(const MyType& b) {
    assert(size() == b.size());
    for (size_t i = 0; i < size(); i++) {
        _data[i] -= b[i];
    }
    return *this;
}

template <typename T>
MyType<T>& InternalFnc<T>::operator-=(double b) {
    for (size_t i = 0; i < size(); i++) {
        _data[i] -= b;
    }
    return *this;
}

template <typename T>
MyType<T>& InternalFnc<T>::operator*=(const MyType& b) {
    assert(size() == b.size());
    for (size_t i = 0; i < size(); i++) {
        _data[i] *= b[i];
    }
    return *this;
}

template <typename T>
MyType<T>& InternalFnc<T>::operator*=(double b) {
    for (size_t i = 0; i < size(); i++) {
        _data[i] *= b;
    }
    return *this;
}

template <typename T>
MyType<T> InternalFnc<T>::operator-() {
    auto out = *this;
    for (size_t i = 0; i < out.size(); i++) {
        out[i] = -out[i];
    }
    return out;
}

template <typename T>
bool InternalFnc<T>::operator==(const InternalFnc<T>& b) const {
    if (size() != b.size()) {
        return false;
    }
    bool result = true;
    for (size_t i = 0; i < size(); i++) {
        result = result && (_data[i] == b[i]);
    }
    return result;
}


template <typename T>
MyType<T> InternalFnc<T>::operator+(const MyType& rhs) {
    auto out = *this;
    return out += rhs;
}
template <typename T>
MyType<T> InternalFnc<T>::operator+(double rhs) {
    auto out = *this;
    return out += rhs;
}
template <typename T>
MyType<T> InternalFnc<T>::operator-(const MyType& rhs) {
    auto out = *this;
    return out -= rhs;
}
template <typename T>
MyType<T> InternalFnc<T>::operator-(double rhs) {
    auto out = *this;
    return out -= rhs;
}
template <typename T>
MyType<T> InternalFnc<T>::operator*(const MyType& rhs) {
    auto out = *this;
    return out *= rhs;
}
template <typename T>
MyType<T> InternalFnc<T>::operator*(double rhs) {
    auto out = *this;
    return out *= rhs;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const InternalFnc<T>& a) {
    os << "[";
    for (size_t i = 0; i < a.size(); i++) {
        os << a[i];
        if (i != a.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

template 
std::ostream& operator<<(std::ostream& os, const Function& a);
template 
std::ostream& operator<<(std::ostream& os, const BlockFunction& a);

template class InternalFnc<double>;
template class InternalFnc<InternalFnc<double>>;
} //end namespace tbem
