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
void hdf_out(const std::string& filename, const Mesh<3>& mesh,
             const std::vector<double>& data) {

    /* Create a new file using default properties. */
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create the data space for the vertices dataset. */
    hsize_t dims[2];
    dims[0] = mesh.facets.size();
    dims[1] = 9;
    hid_t facets_dataspace_id = H5Screate_simple(2, dims, NULL);

    dims[0] = mesh.facets.size() * 3;
    dims[1] = 1;
    hid_t values_dataspace_id = H5Screate_simple(2, dims, NULL); 

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
