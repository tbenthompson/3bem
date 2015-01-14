#ifndef __QMSJOIRNXLKSPAPPZPP_OUTPUT_H
#define __QMSJOIRNXLKSPAPPZPP_OUTPUT_H

#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include "mesh.h"
#include "vec_ops.h"

namespace tbem {

inline bool does_file_exist(const std::string& filename)
{
    std::ifstream infile(filename);
    return infile.good();
}

class Outputter {
public:
    virtual void write_locations(int dim1, int dim2, const void* data_ptr) const = 0;
    virtual void write_values(int n_vars, const std::vector<double>& data) const = 0;
};

class HDFOutputter: public Outputter {
public:
    HDFOutputter(const std::string& filename);
    ~HDFOutputter();

    virtual void write_locations(int dim1, int dim2, const void* data_ptr) const;
    virtual void write_values(int n_vars, const std::vector<double>& data) const;

    const int file_id;
};

template <size_t dim, typename T>
void out_surface(const Outputter& o, const Mesh<dim>& mesh,
                 const std::vector<T> data, int n_vars) {
    o.write_locations(mesh.facets.size(), dim * dim, mesh.facets.data());
    o.write_values(n_vars, reinterpret_vector<double>(data));
}

template <size_t dim, typename T> 
void out_volume(const Outputter& o, const std::vector<Vec<double,dim>>& points,
                const std::vector<T>& data, int n_vars) {
    o.write_locations(points.size(), dim, points.data());
    o.write_values(n_vars, reinterpret_vector<double>(data));
}

} // END namespace tbem

#endif
