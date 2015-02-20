#ifndef __QMSJOIRNXLKSPAPPZPP_OUTPUT_H
#define __QMSJOIRNXLKSPAPPZPP_OUTPUT_H

#include <cassert>
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include "mesh.h"
#include "function.h"
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
    virtual void write_values(const std::string& name, int n_vars,
        const VectorX& data) const = 0;
};

class HDFOutputter: public Outputter {
public:
    HDFOutputter(const std::string& filename);
    ~HDFOutputter();

    virtual void write_locations(int dim1, int dim2, const void* data_ptr) const;
    virtual void write_values(const std::string& name, int n_vars,
        const VectorX& data) const;

    const int file_id;
};

template <size_t dim>
void out_surface(const Outputter& o, const Mesh<dim>& mesh,
                 const BlockVectorX& data) {
    o.write_locations(mesh.facets.size(), dim * dim, mesh.facets.data());
    for (size_t i = 0; i < data.size(); i++) {
        assert(data[i].size() == mesh.n_dofs());
        o.write_values("values" + std::to_string(i), 1, data[i]);
    }
}

template <size_t dim> 
void out_volume(const Outputter& o, const std::vector<Vec<double,dim>>& points,
                const BlockVectorX& data) {
    o.write_locations(points.size(), dim, points.data());
    for (size_t i = 0; i < data.size(); i++) {
        assert(data[i].size() == points.size());
        o.write_values("values" + std::to_string(i), 1, data[i]);
    }
}

} // END namespace tbem

#endif
