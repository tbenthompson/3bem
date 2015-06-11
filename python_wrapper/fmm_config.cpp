#include <boost/python.hpp>
#include "fmm.h"

void export_fmm_config() {
    using namespace boost::python;
    using namespace tbem;
    class_<FMMConfig>("FMMConfig", init<double,size_t,size_t,double,bool>());
}
