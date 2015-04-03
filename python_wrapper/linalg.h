#ifndef __LKJJJJJAJJJAJJ_LINALG_H
#define __LKJJJJJAJJJAJJ_LINALG_H
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "iterable_converter.h"

#include "constraint_matrix.h"
#include "vectorx.h"
#include "dense_operator.h"
#include "block_operator.h"
#include "block_dof_map.h"

template <typename T>
void export_internal_vec(std::string name) {
    using namespace boost::python;
    class_<T>(name.c_str(), init<std::vector<typename T::ValueType>>())
        .def(init<size_t,typename T::ValueType>())
        .def_readonly("storage", &T::storage)
        .def("size", &T::size)
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
}

void export_linalg() {
    using namespace boost::python;
    using namespace tbem;
    export_internal_vec<VectorX>("VectorX");

    class_<std::vector<VectorX>>("VectorOfVectorX")
        .def(vector_indexing_suite<std::vector<VectorX>>());
    VectorFromIterable().from_python<std::vector<VectorX>>();

    export_internal_vec<BlockVectorX>("BlockVectorX");

    class_<DenseOperator>("DenseOperator", no_init)
        .def("apply", &DenseOperator::apply)
        .def("data", &DenseOperator::data,
             return_value_policy<reference_existing_object>());
    class_<BlockDenseOperator>("BlockDenseOperator", no_init)
        .def("apply_scalar", &BlockDenseOperator::apply_scalar)
        .def("apply", &BlockDenseOperator::apply)
        .def("get_block", &BlockDenseOperator::get_block,
             return_value_policy<reference_existing_object>());

    class_<ConstraintEQ>("ConstraintEQ", no_init);
    class_<std::vector<ConstraintEQ>>("VectorOfConstraints", no_init)
        .def(vector_indexing_suite<std::vector<ConstraintEQ>>());
    class_<ConstraintMatrix>("ConstraintMatrix", no_init);
    def("from_constraints", from_constraints);
    def("condense_vector", &condense_vector);
    def("condense_matrix", &condense_matrix);
    def("distribute_vector", &distribute_vector);

    class_<BlockDOFMap>("BlockDOFMap", no_init);
    def("block_dof_map_from_functions", block_dof_map_from_functions);
    def("concatenate", concatenate);
    def("expand", expand);
}

#endif
