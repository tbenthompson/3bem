#include <boost/python.hpp>
// #include "iterable_converter.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "array_from_list.h"

using namespace boost::python;
using namespace tbem;

BOOST_PYTHON_MODULE(lib3bem)
{
    class_<Mesh<2>>("Mesh2D");

    def("line_mesh", line_mesh);
    def("circle_mesh", circle_mesh);
    def("sphere_mesh", sphere_mesh);
    def("rect_mesh", rect_mesh);

    ArrayFromList().from_python<double,2>();
    ArrayFromList().from_python<double,3>();
}
