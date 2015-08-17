#ifndef __ASDKJLAKSDJLJASD_OP_PY_WRAP_H
#define __ASDKJLAKSDJLJASD_OP_PY_WRAP_H

template <typename T, typename BPObj>
void export_operator(BPObj& class_obj)
{
   auto c = class_obj.def("apply", &T::apply)
        .def("n_rows", &T::n_rows)
        .def("n_cols", &T::n_cols)
        .def("clone", &T::clone);
}

#endif
