#include "util.h"
#include "hdf5.h"
#include "mesh.h"

/* equivalent to range(min, max) in python */
std::vector<int> naturals(int min, int max) {
    std::vector<int> indices(max);
    std::iota(indices.begin(), indices.end(), min);
    return indices;
}

/* equivalent to range(0, max) in python */
std::vector<int> naturals(int max) {
    return naturals(0, max);
}

//TODO: Allow multiple output vectors.
template <int dim>
void hdf_out(const std::string& filename, const Mesh<dim>& mesh,
             const std::vector<double>& data) {

    /* Create a new file using default properties. */
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create the data space for the vertices dataset. */
    hsize_t data_dims[2];
    data_dims[0] = mesh.facets.size();
    data_dims[1] = dim * dim;
    hid_t facets_dataspace_id = H5Screate_simple(2, data_dims, NULL);

    data_dims[0] = mesh.facets.size() * dim;
    data_dims[1] = 1;
    hid_t values_dataspace_id = H5Screate_simple(2, data_dims, NULL); 

    /* Create the dataset. */
    hid_t facets_dataset_id = H5Dcreate2(file_id, "/facets", H5T_NATIVE_DOUBLE,
            facets_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t values_dataset_id = H5Dcreate2(file_id, "/values", H5T_NATIVE_DOUBLE,
            values_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write the data.
    H5Dwrite(facets_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             mesh.facets.data());
    H5Dwrite(values_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             data.data());

    /* End access to the dataset and release resources used by it. */
    H5Dclose(facets_dataset_id);
    H5Dclose(values_dataset_id);

    /* Terminate access to the data space. */ 
    H5Sclose(facets_dataspace_id);
    H5Sclose(values_dataspace_id);

    /* Close the file. */
    H5Fclose(file_id);
}

template
void hdf_out<2>(const std::string&, const Mesh<2>&, const std::vector<double>&);
template
void hdf_out<3>(const std::string&, const Mesh<3>&, const std::vector<double>&);
