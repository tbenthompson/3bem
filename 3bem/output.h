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

    template <int dim>
    void hdf_out_surface(const std::string& filename, const Mesh<dim>& mesh,
                         const std::vector<std::vector<double>>& data);

    template <int dim> 
    void hdf_out_volume(const std::string& filename,
                        const std::vector<Vec<double,dim>>& points,
                        const std::vector<std::vector<double>>& data);
} // END namespace tbem

#endif
