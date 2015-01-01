#include <cassert>
#include "output.h"
#include "hdf5.h"
#include "mesh.h"

namespace tbem {

hid_t default_file(const std::string& filename) {
    return H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

void close_file(hid_t file_id) {
    H5Fclose(file_id);
}

void write_locations(hid_t file_id, int data_dim_1, int data_dim_2,
                     const void* data_ptr) {
    // Create the data space for the vertices dataset.
    hsize_t data_dims[2];
    data_dims[0] = data_dim_1;
    data_dims[1] = data_dim_2;
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

void write_values(hid_t file_id, int n_vars,
                  const std::vector<double>& data) {
    assert(data.size() % n_vars == 0);

    int n_dofs = data.size() / n_vars;

    hsize_t data_dims[2];
    data_dims[0] = n_dofs;
    data_dims[1] = n_vars;
    hid_t values_dataspace_id = H5Screate_simple(2, data_dims, NULL); 
    hid_t values_dataset_id = H5Dcreate2(file_id, "/values", H5T_NATIVE_DOUBLE,
            values_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(values_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             data.data());
    H5Dclose(values_dataset_id);
    H5Sclose(values_dataspace_id);
}

std::vector<double> zip(const std::vector<std::vector<double>>& data) {
    int n_vars = data.size();
    int n_dofs = data[0].size();
    std::vector<double> zipped_data(n_vars * n_dofs);
    for (int i = 0; i < n_dofs; i++) {
        for (int j = 0; j < n_vars; j++) {
            zipped_data[i * n_vars + j] = data[j][i];
        }
    }
    return zipped_data;
}

template <int dim>
void hdf_out_surface(const std::string& filename,
                     const Mesh<dim>& mesh,
                     const std::vector<std::vector<double>>& data) {

    // Create a new file using default properties.
    hid_t file_id = default_file(filename);

    write_locations(file_id, mesh.facets.size(), dim * dim, mesh.facets.data());

    auto zipped_values = zip(data);

    write_values(file_id, data.size(), zipped_values);

    // Close the file.
    close_file(file_id);
}

template <int space_dim, int data_dim>
void hdf_out_surface(const std::string& filename, const Mesh<space_dim>& mesh,
                     const std::vector<Vec<double,data_dim>>& data) {
    // Create a new file using default properties.
    hid_t file_id = default_file(filename);

    write_locations(file_id, mesh.facets.size(),
                    space_dim * space_dim, mesh.facets.data());

    write_values(file_id, data_dim, reinterpret_vector<double>(data));

    // Close the file.
    close_file(file_id);
}


template <int dim> 
void hdf_out_volume(const std::string& filename,
                    const std::vector<Vec<double,dim>>& points,
                    const std::vector<std::vector<double>>& data) {

    // Create a new file using default properties.
    hid_t file_id = default_file(filename);

    write_locations(file_id, points.size(), dim, points.data());

    auto zipped_values = zip(data);

    write_values(file_id, data.size(), zipped_values);

    // Close the file.
    close_file(file_id);
}

template
void hdf_out_surface<2>(const std::string&, const Mesh<2>&,
                        const std::vector<std::vector<double>>&);
template
void hdf_out_surface<3>(const std::string&, const Mesh<3>&, 
                        const std::vector<std::vector<double>>&);

template 
void hdf_out_surface<2,1>(const std::string& filename, const Mesh<2>& mesh,
                          const std::vector<Vec<double,1>>& data);
template 
void hdf_out_surface<2,2>(const std::string& filename, const Mesh<2>& mesh,
                          const std::vector<Vec<double,2>>& data);
template
void hdf_out_surface<3,1>(const std::string& filename, const Mesh<3>& mesh,
                          const std::vector<Vec<double,1>>& data);
template 
void hdf_out_surface<3,3>(const std::string& filename, const Mesh<3>& mesh,
                          const std::vector<Vec<double,3>>& data);

template
void hdf_out_volume<2>(const std::string&, const std::vector<Vec<double,2>>&,
                       const std::vector<std::vector<double>>&);
template
void hdf_out_volume<3>(const std::string&, const std::vector<Vec<double,3>>&,
                       const std::vector<std::vector<double>>&);

} // END namespace tbem
