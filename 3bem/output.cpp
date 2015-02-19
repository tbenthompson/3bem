#include <cassert>
#include "output.h"
#include <hdf5.h>
#include "mesh.h"

namespace tbem {

hid_t default_file(const std::string& filename) {
    return H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

void close_file(hid_t file_id) {
    H5Fclose(file_id);
}

HDFOutputter::HDFOutputter(const std::string& filename):
    file_id(default_file(filename))
{}
    
HDFOutputter::~HDFOutputter() {
    close_file(file_id);
}

void HDFOutputter::write_locations(int dim1, int dim2, const void* data_ptr) const {
    // Create the data space for the vertices dataset.
    hsize_t data_dims[2];
    data_dims[0] = dim1;
    data_dims[1] = dim2;
    hid_t locs_dataspace_id = H5Screate_simple(2, data_dims, NULL);

    // Create the dataset.
    hid_t locs_dataset_id = H5Dcreate2(file_id, "/locations", H5T_NATIVE_DOUBLE,
            locs_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write the data.
    H5Dwrite(locs_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             data_ptr);

    // End access to the dataset and release resources used by it. 
    H5Dclose(locs_dataset_id);

    // Terminate access to the data space. 
    H5Sclose(locs_dataspace_id);
}

void HDFOutputter::write_values(const std::string& name, 
        int n_vars, const Function& data) const {
    assert(data.size() % n_vars == 0);

    int n_dofs = data.size() / n_vars;

    std::string dataset_name = "/";
    dataset_name += name;

    hsize_t data_dims[2];
    data_dims[0] = n_dofs;
    data_dims[1] = n_vars;
    hid_t values_dataspace_id = H5Screate_simple(2, data_dims, NULL); 
    hid_t values_dataset_id = H5Dcreate2(file_id, dataset_name.c_str(),
            H5T_NATIVE_DOUBLE, values_dataspace_id, H5P_DEFAULT, H5P_DEFAULT,
            H5P_DEFAULT);
    H5Dwrite(values_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             data.data());
    H5Dclose(values_dataset_id);
    H5Sclose(values_dataspace_id);
}


} // END namespace tbem
