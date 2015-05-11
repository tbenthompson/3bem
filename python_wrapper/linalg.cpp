#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "iterable_converter.h"
#include "op_wrap.h"

#include "dense_operator.h"
#include "sparse_operator.h"
namespace p = boost::python;

void export_linalg() {
    using namespace tbem;

    export_operator<DenseOperator>(
        p::class_<DenseOperator>("DenseOperator",
            p::init<size_t, size_t, std::vector<double>>())
        .def("data", &DenseOperator::data,
             p::return_value_policy<p::return_by_value>())
    );

    export_operator<SparseOperator>(
        p::class_<SparseOperator>("SparseOperator", p::no_init)
        .def("nnz", &SparseOperator::nnz)
        .add_property("values", make_getter(
                &SparseOperator::values, 
                p::return_value_policy<p::return_by_value>()))
        .add_property("column_indices", make_getter(
                &SparseOperator::column_indices, 
                p::return_value_policy<p::return_by_value>()))
        .add_property("row_ptrs", make_getter(
                &SparseOperator::row_ptrs, 
                p::return_value_policy<p::return_by_value>()))
    );
}
