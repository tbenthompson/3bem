#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "iterable_converter.h"

#include "vectorx.h"
#include "dense_operator.h"
#include "sparse_operator.h"
#include "block_operator.h"
#include "block_dof_map.h"
#include "dense_builder.h"
namespace p = boost::python;

template <typename T>
void export_internal_vec(std::string name) {
    p::class_<T>(name.c_str(), p::init<std::vector<typename T::ValueType>>())
        .def(p::init<size_t,typename T::ValueType>())
        .add_property("storage", p::make_getter(
            &T::storage, p::return_value_policy<p::return_by_value>()
        ))
        .def("size", &T::size)
        .def(p::self + p::self)
        .def(p::self + double())
        .def(p::self - p::self)
        .def(p::self - double())
        .def(p::self * p::self)
        .def(p::self * double())
        .def(p::self += p::self)
        .def(p::self += double())
        .def(p::self -= p::self)
        .def(p::self -= double())
        .def(p::self *= p::self)
        .def(p::self *= double())
        .def(-p::self);
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
    using namespace tbem;
    export_internal_vec<VectorX>("VectorX");

    p::class_<std::vector<VectorX>>("VectorOfVectorX")
        .def(p::vector_indexing_suite<std::vector<VectorX>>());
    VectorFromIterable().from_python<std::vector<VectorX>>();

    export_internal_vec<BlockVectorX>("BlockVectorX");

    p::class_<OperatorIWrap, boost::noncopyable>("OperatorI")
        .def("apply", p::pure_virtual(&OperatorI::apply))
        .def("n_rows", p::pure_virtual(&OperatorI::n_rows))
        .def("n_cols", p::pure_virtual(&OperatorI::n_cols));
    BlockVectorX (BlockOperatorI::*apply1)(const BlockVectorX&) const =
        &BlockOperatorI::apply;
    VectorX (BlockOperatorI::*apply2)(const VectorX&) const =
        &BlockOperatorI::apply;
    p::class_<BlockOperatorIWrap, boost::noncopyable>("BlockOperatorI")
        .def("apply", p::pure_virtual(apply1))
        .def("apply", p::pure_virtual(apply2))
        .def("n_total_rows", p::pure_virtual(&BlockOperatorI::n_total_rows))
        .def("n_total_cols", p::pure_virtual(&BlockOperatorI::n_total_cols))
        .def("n_block_rows", p::pure_virtual(&BlockOperatorI::n_block_rows))
        .def("n_block_cols", p::pure_virtual(&BlockOperatorI::n_block_cols));

    p::class_<DenseOperator, p::bases<OperatorI>>("DenseOperator", p::no_init)
        .def("data", &DenseOperator::data,
             p::return_value_policy<p::reference_existing_object>());
    p::class_<SparseOperator, p::bases<OperatorI>>("SparseOperator", p::no_init);
    p::class_<BlockDenseOperator, p::bases<BlockOperatorI>>(
        "BlockDenseOperator", p::no_init
    ).def("get_block", &BlockDenseOperator::get_block,
        p::return_value_policy<p::reference_existing_object>());
    p::class_<BlockSparseOperator, p::bases<BlockOperatorI>>(
        "BlockSparseOperator", p::no_init
    ).def("get_block", &BlockSparseOperator::get_block,
        p::return_value_policy<p::reference_existing_object>());

    p::class_<BlockDOFMap>("BlockDOFMap", p::no_init)
        .def_readonly("n_components", &BlockDOFMap::n_components)
        .def_readonly("n_dofs", &BlockDOFMap::n_dofs)
        .add_property("start_positions", p::make_getter(
            &BlockDOFMap::start_positions, p::return_value_policy<p::return_by_value>()
        ));
    p::def("build_block_dof_map", build_block_dof_map);
    p::def("block_dof_map_from_functions", block_dof_map_from_functions);
    p::def("concatenate", concatenate);
    p::def("expand", expand);
}
