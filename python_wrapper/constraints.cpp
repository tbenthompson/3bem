#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "iterable_converter.h"
#include "constraint_matrix.h"

void export_constraints() {
    using namespace boost::python;
    using namespace tbem;

    class_<ConstraintEQ>("ConstraintEQ", no_init);
    class_<std::vector<ConstraintEQ>>("VectorOfConstraints", no_init)
        .def(vector_indexing_suite<std::vector<ConstraintEQ>>());
    VectorFromIterable().from_python<std::vector<ConstraintEQ>>();

    class_<ConstraintMatrix>("ConstraintMatrix", no_init);
    def("shift_constraints", shift_constraints);
    def("from_constraints", from_constraints);
    def("condense_vector", &condense_vector);
    def("condense_matrix", &condense_matrix);
    def("distribute_vector", &distribute_vector);
}
