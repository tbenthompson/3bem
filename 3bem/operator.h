#ifndef ADKJASdLL_OPERATOR_H
#define ADKJASdLL_OPERATOR_H

#include "vectorx.h"

namespace tbem {


/* This NotImplemented exception class is justified for implemented an interface
 * in order to avoid using an abstract base class(ABC). ABCs force you to
 * use a pointer or reference in order to refer to them, meaning that any
 * usage of an interface requires pointers. This is unpleasant as pointers should
 * be avoided as much as possible. To get around this, the base class should 
 * simply not be abstract. In order to be clear that the base class should not
 * actually be used, NotImplemented exceptions are thrown.
 */
struct NotImplemented {};

struct OperatorShape {
    const size_t n_rows;
    const size_t n_cols;
};

struct OperatorI {
    virtual size_t n_rows() const {throw NotImplemented();};
    virtual size_t n_cols() const {throw NotImplemented();};
    virtual VectorX apply(const VectorX& x) const {throw NotImplemented();};
};

} // end namespace tbem

#endif
