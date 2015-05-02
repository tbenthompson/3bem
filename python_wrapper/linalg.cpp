#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "iterable_converter.h"

#include "vectorx.h"
#include "dense_operator.h"
#include "sparse_operator.h"
#include "block_operator.h"
#include "block_dof_map.h"
#include "dense_builder.h"

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

namespace tbem {
    struct OperatorIWrap: OperatorI, boost::python::wrapper<OperatorI> 
    {
        virtual size_t n_rows() const {return this->get_override("n_rows")();}
        virtual size_t n_cols() const {return this->get_override("n_cols")();}

        virtual VectorX apply(const VectorX& x) const 
        {
            return this->get_override("apply")();
        }
    };

    struct BlockOperatorIWrap: BlockOperatorI, boost::python::wrapper<BlockOperatorI> 
    {
        virtual size_t n_total_rows() const 
        {
            return this->get_override("n_total_rows")();
        }
        virtual size_t n_total_cols() const 
        {
            return this->get_override("n_total_cols")();
        }
        virtual size_t n_block_rows() const 
        {
            return this->get_override("n_block_rows")();
        }
        virtual size_t n_block_cols() const 
        {
            return this->get_override("n_block_cols")();
        }

        virtual VectorX apply(const VectorX& x) const 
        {
            return this->get_override("apply")();
        }
        virtual BlockVectorX apply(const BlockVectorX& x) const 
        {
            return this->get_override("apply")();
        }
    };
}

void export_linalg() {
    using namespace boost::python;
    using namespace tbem;
    export_internal_vec<VectorX>("VectorX");

    class_<std::vector<VectorX>>("VectorOfVectorX")
        .def(vector_indexing_suite<std::vector<VectorX>>());
    VectorFromIterable().from_python<std::vector<VectorX>>();

    export_internal_vec<BlockVectorX>("BlockVectorX");

    class_<OperatorIWrap, boost::noncopyable>("OperatorI")
        .def("apply", pure_virtual(&OperatorI::apply))
        .def("n_rows", pure_virtual(&OperatorI::n_rows))
        .def("n_cols", pure_virtual(&OperatorI::n_cols));
    BlockVectorX (BlockOperatorI::*apply1)(const BlockVectorX&) const =
        &BlockOperatorI::apply;
    VectorX (BlockOperatorI::*apply2)(const VectorX&) const =
        &BlockOperatorI::apply;
    class_<BlockOperatorIWrap, boost::noncopyable>("BlockOperatorI")
        .def("apply", pure_virtual(apply1))
        .def("apply", pure_virtual(apply2))
        .def("n_total_rows", pure_virtual(&BlockOperatorI::n_total_rows))
        .def("n_total_cols", pure_virtual(&BlockOperatorI::n_total_cols))
        .def("n_block_rows", pure_virtual(&BlockOperatorI::n_block_rows))
        .def("n_block_cols", pure_virtual(&BlockOperatorI::n_block_cols));

    class_<DenseOperator, bases<OperatorI>>("DenseOperator", no_init)
        .def("data", &DenseOperator::data,
             return_value_policy<reference_existing_object>());
    class_<SparseOperator, bases<OperatorI>>("SparseOperator", no_init);
    class_<BlockDenseOperator, bases<BlockOperatorI>>("BlockDenseOperator", no_init)
        .def("get_block", &BlockDenseOperator::get_block,
             return_value_policy<reference_existing_object>());
    class_<BlockSparseOperator, bases<BlockOperatorI>>("BlockSparseOperator", no_init)
        .def("get_block", &BlockSparseOperator::get_block,
             return_value_policy<reference_existing_object>());

    class_<BlockDOFMap>("BlockDOFMap", no_init);
    def("build_block_dof_map", build_block_dof_map);
    def("block_dof_map_from_functions", block_dof_map_from_functions);
    def("concatenate", concatenate);
    def("expand", expand);
}
