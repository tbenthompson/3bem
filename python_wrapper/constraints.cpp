#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "iterable_converter.h"
#include "constraint_matrix.h"
#include "row_zero_distributor.h"
#include "op_wrap.h"

namespace p = boost::python;

void export_constraints() {
    using namespace tbem;

    p::class_<LinearTerm>("LinearTerm", p::init<size_t,double>())
        .def_readonly("dof", &LinearTerm::dof)
        .def_readonly("weight", &LinearTerm::weight);
    p::class_<std::vector<LinearTerm>>("VectorOfLinearTerm", p::no_init)
        .def(p::vector_indexing_suite<std::vector<LinearTerm>>());
    VectorFromIterable().from_python<std::vector<LinearTerm>>();

    p::class_<ConstraintEQ>("ConstraintEQ", p::init<std::vector<LinearTerm>,double>())
        .def_readonly("terms", &ConstraintEQ::terms)
        .def_readonly("rhs", &ConstraintEQ::rhs);
    p::class_<std::vector<ConstraintEQ>>("VectorOfConstraints", p::no_init)
        .def(p::vector_indexing_suite<std::vector<ConstraintEQ>>());
    VectorFromIterable().from_python<std::vector<ConstraintEQ>>();

    p::class_<ConstraintMatrix>("ConstraintMatrix", p::no_init);
    p::def("shift_constraints", &shift_constraints);
    p::def("homogenize_constraints", &homogenize_constraints);
    p::def("from_constraints", &from_constraints);
    p::def("condense_vector", &condense_vector);

    DenseOperator (*condense_dense)(const ConstraintMatrix&,
        const ConstraintMatrix&, const DenseOperator&) = &condense_matrix;
    p::def("condense_matrix", condense_dense);
    SparseOperator (*condense_sparse)(const ConstraintMatrix&,
        const ConstraintMatrix&, const SparseOperator&) = &condense_matrix;
    p::def("condense_matrix", condense_sparse);

    p::def("distribute_vector", &distribute_vector);

    p::def("identify_ignored_dofs", &identify_ignored_dofs); 
    p::def("distribute_row_zeros", &distribute_row_zeros);

    auto rzd_wrap = p::class_<RowZeroDistributor, p::bases<OperatorI>, boost::noncopyable>(
        "RowZeroDistributor", 
        p::init<const ConstraintMatrix&, const OperatorI&>()
    );
    export_operator<RowZeroDistributor>(rzd_wrap);
}
