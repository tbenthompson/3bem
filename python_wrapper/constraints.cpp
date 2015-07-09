#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "iterable_converter.h"
#include "constraint_matrix.h"

void export_constraints() {
    using namespace boost::python;
    using namespace tbem;

    class_<LinearTerm>("LinearTerm", init<size_t,double>())
        .def_readonly("dof", &LinearTerm::dof)
        .def_readonly("weight", &LinearTerm::weight);
    class_<std::vector<LinearTerm>>("VectorOfLinearTerm", no_init)
        .def(vector_indexing_suite<std::vector<LinearTerm>>());
    VectorFromIterable().from_python<std::vector<LinearTerm>>();

    class_<ConstraintEQ>("ConstraintEQ", init<std::vector<LinearTerm>,double>())
        .def_readonly("terms", &ConstraintEQ::terms)
        .def_readonly("rhs", &ConstraintEQ::rhs);
    class_<std::vector<ConstraintEQ>>("VectorOfConstraints", no_init)
        .def(vector_indexing_suite<std::vector<ConstraintEQ>>());
    VectorFromIterable().from_python<std::vector<ConstraintEQ>>();

    class_<ConstraintMatrix>("ConstraintMatrix", no_init);
    def("shift_constraints", shift_constraints);
    def("homogenize_constraints", homogenize_constraints);
    def("from_constraints", from_constraints);
    def("condense_vector", &condense_vector);

    DenseOperator (*condense_dense)(const ConstraintMatrix&,
        const ConstraintMatrix&, const DenseOperator&) = &condense_matrix;
    def("condense_matrix", condense_dense);
    SparseOperator (*condense_sparse)(const ConstraintMatrix&,
        const ConstraintMatrix&, const SparseOperator&) = &condense_matrix;
    def("condense_matrix", condense_sparse);

    def("distribute_vector", &distribute_vector);
}
