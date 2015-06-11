#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <iostream>

namespace p = boost::python;
namespace np = boost::numpy;

void export_util();
void export_linalg();
void export_constraints();
void export_fmm_config();
template <size_t dim>
void export_dimension();

boost::python::scope start_module(std::string name) {
    using namespace boost::python;
    std::string fullname = extract<std::string>(scope().attr("__name__") + "." + name);
    object module(handle<>(borrowed(PyImport_AddModule(fullname.c_str()))));
    scope().attr(name.c_str()) = module;
    return module;
}

BOOST_PYTHON_MODULE(_tbempy)
{
    np::initialize();

    export_util();
    export_linalg();
    export_constraints();
    export_fmm_config();

    {
        p::scope scope2D = start_module("TwoD");
        scope2D.attr("dim") = 2;
        export_dimension<2>();
    }

    {
        p::scope scope3D = start_module("ThreeD");
        scope3D.attr("dim") = 3;
        export_dimension<3>();
    }
}
