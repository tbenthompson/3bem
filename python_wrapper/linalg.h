#ifndef __LKJJJJJAJJJAJJ_LINALG_H
#define __LKJJJJJAJJJAJJ_LINALG_H
#include <boost/python.hpp>

#include "constraint_matrix.h"
#include "vectorx.h"
#include "dense_operator.h"
#include "block_operator.h"

void export_linalg() {
    using namespace boost::python;
    using namespace tbem;
    class_<VectorX>("VectorX", init<std::vector<double>>())
        .def_readonly("storage", &VectorX::storage)
        .def("size", &VectorX::size)
        .def(self + self)
        .def(self + double())
        .def(self - self)
        .def(self - double())
        .def(self * self)
        .def(self * double())
        .def(self += self)
        .def(self += double())
        .def(self -= self)
        .def(self -= double())
        .def(self *= self)
        .def(self *= double())
        .def(-self);

    class_<DenseOperator>("DenseOperator", no_init)
        .def("apply", &DenseOperator::apply)
        .def("data", &DenseOperator::data,
             return_value_policy<reference_existing_object>());
    class_<BlockDenseOperator>("BlockDenseOperator", no_init)
        .def("apply_scalar", &BlockDenseOperator::apply_scalar)
        .def("get_block", &BlockDenseOperator::get_block,
             return_value_policy<reference_existing_object>());

    class_<std::vector<ConstraintEQ>>("VectorOfConstraints");
    class_<ConstraintMatrix>("ConstraintMatrix");
    def("from_constraints", from_constraints);
    def("condense_vector", &condense_vector);
    def("condense_matrix", &condense_matrix);
    def("distribute_vector", &distribute_vector);
}

#endif
