#ifndef __ASDKJWEORIOQ_VEC_ARRAYS_H 
#define __ASDKJWEORIOQ_VEC_ARRAYS_H 
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include "numbers.h"

namespace tbem {

/* [Note: Justification of overloading std::array operations]
 * Having a 3-vector class is useful because it reduces instances where
 * the same code has to be replicated 3 times for each element.
 *
 * Writing a new 3-vector class is not hard at all, but it can't hurt
 * to rely on an existing standard class. 
 *
 * This avoids large libraries like Eigen just to have a 3-vector class.
 *
 * This can later be extended to work on std::vector or arbitrary length 
 * std::array.
 *
 * C++11 templated typedefs allow this Vec3 "class" to be templated on float
 * vs double.
 *
 * It would be interesting to extend this to use expression templates for
 * lazy evaluation of results.
 *
 * Extending this to use recursion to do 1D,2D,3D all in one.
 */

template <typename T, size_t dim>
using Vec = std::array<T,dim>;

template <typename T>
using Vec3 = Vec<T,3>;

template <typename T>
using Vec2 = Vec<T,2>;

} // END namespace tbem

#endif
