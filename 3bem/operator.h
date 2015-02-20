#ifndef ADKJASdLL_OPERATOR_H
#define ADKJASdLL_OPERATOR_H

namespace tbem {

template <typename T>
struct InternalVec;
typedef InternalVec<double> VectorX;
typedef InternalVec<InternalVec<double>> BlockVectorX;

struct OperatorI {
    virtual size_t n_rows() const = 0;
    virtual size_t n_cols() const = 0;
    virtual VectorX apply(const VectorX& x) const = 0;
};

struct BlockOperatorI {
    virtual size_t n_block_rows() const = 0;
    virtual size_t n_block_cols() const = 0;
    virtual size_t n_total_rows() const = 0;
    virtual size_t n_total_cols() const = 0;
    virtual BlockVectorX apply(const BlockVectorX& x) const = 0;
};

} // end namespace tbem

#endif
