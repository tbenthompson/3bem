#ifndef __QMSJOIRNXLKSPAPPZPP_OUTPUT_H
#define __QMSJOIRNXLKSPAPPZPP_OUTPUT_H

#include <array>
#include <vector>

namespace tbem {

//Forward declaration of Vec
template <typename T, unsigned long dim>
using Vec = std::array<T,dim>;

//Forward declaration of Mesh.
template <typename T, int dim>
struct MeshField;
template <int dim>
using Mesh = MeshField<Vec<double,dim>,dim>;
// 
// class Outputter {
// public:
//     virtual void write_locations(int dim_1, int dim_2, const void* data_ptr) = 0;
//     virtual void write_values(int n_vars, const std::vector<double>& data) = 0;
// };
// 
// class HDFOutputter {
// public:
//     HDFOutputter(const std::string& filename);
// 
//     virtual void write_locations(int dim_1, int dim_2, const void* data_ptr);
//     virtual void write_values(int n_vars, const std::vector<double>& data);
// 
//     int file_id;
// };

template <int dim>
void hdf_out_surface(const std::string& filename, const Mesh<dim>& mesh,
                     const std::vector<std::vector<double>>& data);

template <int space_dim, int data_dim>
void hdf_out_surface(const std::string& filename, const Mesh<space_dim>& mesh,
                     const std::vector<Vec<double,data_dim>>& data);

template <int dim> 
void hdf_out_volume(const std::string& filename,
                    const std::vector<Vec<double,dim>>& points,
                    const std::vector<std::vector<double>>& data);

} // END namespace tbem

#endif
