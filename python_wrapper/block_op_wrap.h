#ifndef __ASDKJLAKSDJLJASD_BLOCK_OP_PY_WRAP_H
#define __ASDKJLAKSDJLJASD_BLOCK_OP_PY_WRAP_H

template <typename T, typename BPObj>
void export_block_operator(BPObj& class_obj)
{
   class_obj.def("apply", &T::apply)
        .def("n_total_rows", &T::n_total_rows)
        .def("n_total_cols", &T::n_total_cols)
        .def("n_block_rows", &T::n_block_rows)
        .def("n_block_cols", &T::n_block_cols);
}

#endif
