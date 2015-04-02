#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "array_from_list.h"

#include "mesh.h"
#include "quadrature.h"
#include "mesh_gen.h"
#include "constraint_builder.h"
#include "constraint_matrix.h"
#include "basis.h"

#include "laplace_kernels.h"
#include "identity_kernels.h"

#include "dense_builder.h"

#include "petsc_facade.h"

using namespace boost::python;
using namespace tbem;

template <size_t dim>
VectorX interpolate_wrapper(const Mesh<dim>& mesh, const object& fnc) 
{
    return interpolate<2>(mesh,
        [&](const Vec<double,dim>& x) {
            double res = extract<double>(fnc(x[0], x[1]));
            return res;
        });
}

VectorX solve_system_wrapper(const VectorX& rhs, double tolerance, const object& fnc) 
{
    return solve_system(rhs, tolerance,
        [&](std::vector<double>& x, std::vector<double>& y) {
            VectorX out = extract<VectorX>(fnc(VectorX(x)));
            std::copy(out.storage.begin(), out.storage.end(), y.begin());
        });
}

BOOST_PYTHON_MODULE(tbempy)
{
    ArrayFromIterable().from_python<double,2>();
    VectorFromIterable().from_python<std::vector<double>>();
    class_<std::vector<double>>("VectorOfDoubles")
        .def(vector_indexing_suite<std::vector<double>>());

    class_<VertexIterator<2>>("VertexIterator", no_init);
    class_<Mesh<2>>("Mesh")
        .def("begin", &Mesh<2>::begin)
        .def("n_dofs", &Mesh<2>::n_dofs);

    def("line_mesh", line_mesh);
    def("circle_mesh", circle_mesh);
    def("sphere_mesh", sphere_mesh);
    def("rect_mesh", rect_mesh);

    class_<QuadStrategy<2>>("QuadStrategy", init<int,int,int,double,double>());

    class_<OverlapMap<2>>("OverlapMap"); 
    def("mesh_continuity", mesh_continuity<2>);

    class_<std::vector<ConstraintEQ>>("VectorOfConstraints");
    def("convert_to_constraints", convert_to_constraints<2>);

    class_<ConstraintMatrix>("ConstraintMatrix");
    def("from_constraints", from_constraints);
    def("condense_vector", &condense_vector);
    def("condense_matrix", &condense_matrix);
    def("distribute_vector", &distribute_vector);

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
        

    def("interpolate", interpolate_wrapper<2>); 

    class_<Kernel<2,1,1>, boost::noncopyable>("Kernel", no_init);
    class_<IdentityScalar<2>, bases<Kernel<2,1,1>>>("IdentityScalar");
    class_<LaplaceSingle<2>, bases<Kernel<2,1,1>>>("LaplaceSingle");
    class_<LaplaceDouble<2>, bases<Kernel<2,1,1>>>("LaplaceDouble");

    class_<BoundaryIntegral<2,1,1>>("BoundaryIntegralScalar", no_init);
    def("make_boundary_integral", make_boundary_integral<2,1,1>);

    class_<DenseOperator>("DenseOperator", no_init)
        .def("apply", &DenseOperator::apply)
        .def("data", &DenseOperator::data,
             return_value_policy<reference_existing_object>());
    class_<BlockDenseOperator>("BlockDenseOperator", no_init)
        .def("apply_scalar", &BlockDenseOperator::apply_scalar)
        .def("get_block", &BlockDenseOperator::get_block,
             return_value_policy<reference_existing_object>());
    def("mesh_to_mesh_operator", mesh_to_mesh_operator<2,1,1>);
    def("mass_operator", mass_operator<2,1,1>);


    def("solve_system", solve_system_wrapper);
}
