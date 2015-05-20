#ifndef __1231231231898972_SPARSE_OPERATOR_H
#define __1231231231898972_SPARSE_OPERATOR_H
#include <vector>
#include <assert.h>
#include <iostream>
#include "operator.h"

namespace tbem {

struct DenseOperator;

struct MatrixEntry 
{
    size_t loc[2];
    double value;

    size_t row() const {return loc[0];}
    size_t col() const {return loc[1];}
};

struct SparseOperator: public OperatorI
{
    const OperatorShape shape;
    const std::vector<double> values;
    const std::vector<size_t> column_indices;
    const std::vector<size_t> row_ptrs;

    SparseOperator(size_t n_rows, size_t n_cols,
        const std::vector<double>& values,
        const std::vector<size_t>& column_indices,
        const std::vector<size_t>& row_ptrs);

    virtual size_t n_rows() const {return shape.n_rows;}
    virtual size_t n_cols() const {return shape.n_cols;}
    size_t nnz() const {return row_ptrs.back();}

    virtual std::vector<double> apply(const std::vector<double>& x) const;
    DenseOperator to_dense() const; 

    /* Left multiply a sparse matrix by a dense matrix getting a dense matrix
     * as output:
     * (dense) * (sparse) = out
     * Simple dense matrix multiplication is an O(n^3) operation.
     * These are implemented as O(n^2 + nZ) where Z is the nnz of the sparse matrix.
     */
    DenseOperator left_multiply_with_dense(const DenseOperator& other) const;

    /* Right multiply -- see comments on left_multiply_with_dense
     */
    DenseOperator right_multiply_with_dense(const DenseOperator& other) const;

    DenseOperator add_with_dense(const DenseOperator& other) const;

    SparseOperator right_multiply(const SparseOperator& other) const;

    static SparseOperator csr_from_coo(size_t n_rows, size_t n_cols,
        const std::vector<MatrixEntry>& entries);
};


} // end namespace tbem

#endif
