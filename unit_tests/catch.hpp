#ifndef __3BEM_WRAPPER_OF_CATCH_HPP
#define __3BEM_WRAPPER_OF_CATCH_HPP

#include "_catch.hpp"
#include <cmath>
#include "vec_ops.h"

namespace tbem {

template <typename T1, typename T2>
void REQUIRE_ARRAY_EQUAL(const T1& a1, const T2& a2, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        REQUIRE(a1[i] == a2[i]);
    }
}

template <typename T1, typename T2>
void REQUIRE_CLOSE(const T1& a, const T2& b, double epsilon) 
{
    REQUIRE(all(fabs(a - b) < epsilon));
}

template <typename T1, typename T2>
void REQUIRE_ARRAY_CLOSE(const T1& a1, const T2& a2, size_t n, double epsilon)
{
    for (size_t i = 0; i < n; i++) {
        REQUIRE_CLOSE(a1[i], a2[i], epsilon);
    }
}

template <typename T1, typename T2>
void CHECK_ARRAY_EQUAL(const T1& a1, const T2& a2, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        CHECK(a1[i] == a2[i]);
    }
}

template <typename T1, typename T2>
void CHECK_CLOSE(const T1& a, const T2& b, double epsilon) 
{
    CHECK(all(fabs(a - b) < epsilon));
}

template <typename T1, typename T2>
void CHECK_ARRAY_CLOSE(const T1& a1, const T2& a2, size_t n, double epsilon)
{
    for (size_t i = 0; i < n; i++) {
        CHECK_CLOSE(a1[i], a2[i], epsilon);
    }
}

} //end namespace tbem


#endif
