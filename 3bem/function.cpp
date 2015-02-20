#include "function.h"
#include <cassert>
#include <iostream>

namespace tbem {

template <typename T>
using ValueType = T;
template <typename T>
using MyType = InternalVec<T>;
template <typename T>
using ContainerType = std::vector<T>;

template <typename T>
InternalVec<T>::InternalVec() {}

template <typename T>
InternalVec<T>::InternalVec(size_t n_elements):
    _data(n_elements)
{}

template <typename T>
InternalVec<T>::InternalVec(size_t n_elements, const T& value):
    _data(n_elements, value)
{}

template <typename T>
InternalVec<T>::InternalVec(const ContainerType& data):
    _data(data)
{}

template <typename T>
InternalVec<T>::InternalVec(std::initializer_list<T> s):
    _data(s)
{}

template <typename T>
void InternalVec<T>::resize(size_t new_size) {
    _data.resize(new_size);
}

template <typename T>
ValueType<T>& InternalVec<T>::operator[] (size_t idx) {
    return _data[idx];
}

template <typename T>
const ValueType<T>& InternalVec<T>::operator[] (size_t idx) const {
    return _data[idx];
}

template <typename T>
typename ContainerType<T>::iterator InternalVec<T>::begin() {
    return _data.begin();
}

template <typename T>
typename ContainerType<T>::const_iterator InternalVec<T>::begin() const {
    return _data.begin();
}

template <typename T>
typename ContainerType<T>::iterator InternalVec<T>::end() {
    return _data.end();
}

template <typename T>
typename ContainerType<T>::const_iterator InternalVec<T>::end() const {
    return _data.end();
}

template <typename T>
ValueType<T>* InternalVec<T>::data() {
    return _data.data();
}

template <typename T>
const ValueType<T>* InternalVec<T>::data() const {
    return _data.data();
}

template <typename T>
size_t InternalVec<T>::size() const {
    return _data.size();
}

template <typename T>
MyType<T>& InternalVec<T>::operator+=(const MyType& b) {
    assert(size() == b.size());
    for (size_t i = 0; i < size(); i++) {
        _data[i] += b[i];
    }
    return *this;
}

template <typename T>
MyType<T>& InternalVec<T>::operator+=(double b) {
    for (size_t i = 0; i < size(); i++) {
        _data[i] += b;
    }
    return *this;
}

template <typename T>
MyType<T>& InternalVec<T>::operator-=(const MyType& b) {
    assert(size() == b.size());
    for (size_t i = 0; i < size(); i++) {
        _data[i] -= b[i];
    }
    return *this;
}

template <typename T>
MyType<T>& InternalVec<T>::operator-=(double b) {
    for (size_t i = 0; i < size(); i++) {
        _data[i] -= b;
    }
    return *this;
}

template <typename T>
MyType<T>& InternalVec<T>::operator*=(const MyType& b) {
    assert(size() == b.size());
    for (size_t i = 0; i < size(); i++) {
        _data[i] *= b[i];
    }
    return *this;
}

template <typename T>
MyType<T>& InternalVec<T>::operator*=(double b) {
    for (size_t i = 0; i < size(); i++) {
        _data[i] *= b;
    }
    return *this;
}

template <typename T>
MyType<T> InternalVec<T>::operator-() {
    auto out = *this;
    for (size_t i = 0; i < out.size(); i++) {
        out[i] = -out[i];
    }
    return out;
}

template <typename T>
bool InternalVec<T>::operator==(const InternalVec<T>& b) const {
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
MyType<T> InternalVec<T>::operator+(const MyType& rhs) {
    auto out = *this;
    return out += rhs;
}
template <typename T>
MyType<T> InternalVec<T>::operator+(double rhs) {
    auto out = *this;
    return out += rhs;
}
template <typename T>
MyType<T> InternalVec<T>::operator-(const MyType& rhs) {
    auto out = *this;
    return out -= rhs;
}
template <typename T>
MyType<T> InternalVec<T>::operator-(double rhs) {
    auto out = *this;
    return out -= rhs;
}
template <typename T>
MyType<T> InternalVec<T>::operator*(const MyType& rhs) {
    auto out = *this;
    return out *= rhs;
}
template <typename T>
MyType<T> InternalVec<T>::operator*(double rhs) {
    auto out = *this;
    return out *= rhs;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const InternalVec<T>& a) {
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

template class InternalVec<double>;
template class InternalVec<InternalVec<double>>;
} //end namespace tbem
