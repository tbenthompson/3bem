#ifndef __FMM_H
#define __FMM_H

#include <vector>

class FMMInfo {

};

class Octree;

std::vector<double> P2M(Octree& oct, int order, std::vector<double> p_values);

#endif
