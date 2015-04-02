#ifndef __ALGS_H
#define __ALGS_H
#include <vector>
#include <algorithm>
#include <iostream>

namespace tbem {

template <typename T, typename Compare>
std::vector<int> sort_permutation(std::vector<T> const& vec,
                                  Compare compare)
{
    std::vector<int> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&](int i, int j){return compare(vec[i], vec[j]);});
    return p;
}

template <typename T>
std::vector<T> apply_permutation(std::vector<T> const& vec,
                                 std::vector<int> const& p)
{
    std::vector<T> sorted(p.size());
    std::transform(p.begin(), p.end(), sorted.begin(), [&](int i){return vec[i];});
    return sorted;
}

// method to seperate bits from a given integer 3 positions apart
inline uint64_t split_by_3(unsigned int a){
    uint64_t x = a & 0x1fffff; // we only look at the first 21 bits
    // shift left 32 bits, OR with self, and 
    // 00011111000000000000000000000000000000001111111111111111
    x = (x | x << 32) & 0x1f00000000ffff;  
    // shift left 32 bits, OR with self, and 
    // 00011111000000000000000011111111000000000000000011111111
    x = (x | x << 16) & 0x1f0000ff0000ff;  
    // shift left 32 bits, OR with self, and
    // 0001000000001111000000001111000000001111000000001111000000000000
    x = (x | x << 8) & 0x100f00f00f00f00f;
    // shift left 32 bits, OR with self, and 
    // 0001000011000011000011000011000011000011000011000011000100000000
    x = (x | x << 4) & 0x10c30c30c30c30c3; 
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}

inline uint64_t morton_encode(unsigned int x, unsigned int y, unsigned int z) {
    uint64_t answer = 0 | split_by_3(x) | split_by_3(y) << 1 | split_by_3(z) << 2;
    return answer;
}

} // END namespace tbem
#endif
