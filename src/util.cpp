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

void hdf_out(const std::string& filename, const Mesh& mesh,
             const std::vector<double>& data) {

    /* Create a new file using default properties. */
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create the data space for the vertices dataset. */
    hsize_t dims[2];
    dims[0] = mesh.faces.size();
    dims[1] = 3;
    hid_t faces_dataspace_id = H5Screate_simple(2, dims, NULL);

    dims[0] = mesh.vertices.size(); 
    dims[1] = 3; 
    hid_t verts_dataspace_id = H5Screate_simple(2, dims, NULL);

    dims[1] = 1;
    hid_t values_dataspace_id = H5Screate_simple(2, dims, NULL); 



    /* Create the dataset. */
    hid_t faces_dataset_id = H5Dcreate2(file_id, "/faces", H5T_NATIVE_INT,
            faces_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t verts_dataset_id = H5Dcreate2(file_id, "/vertices", H5T_NATIVE_DOUBLE,
            verts_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t values_dataset_id = H5Dcreate2(file_id, "/values", H5T_NATIVE_DOUBLE,
            values_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write the data.
    H5Dwrite(faces_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             mesh.faces.data());
    H5Dwrite(verts_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             mesh.vertices.data());
    H5Dwrite(values_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             data.data());

    /* End access to the dataset and release resources used by it. */
    H5Dclose(faces_dataset_id);
    H5Dclose(verts_dataset_id);
    H5Dclose(values_dataset_id);

    /* Terminate access to the data space. */ 
    H5Sclose(faces_dataspace_id);
    H5Sclose(verts_dataspace_id);
    H5Sclose(values_dataspace_id);

    /* Close the file. */
    H5Fclose(file_id);
}
