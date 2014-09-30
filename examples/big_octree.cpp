#include "octree.h"
#include "test_shared.h"

int main() {
    const int n = 3e7;
    auto pts = random_pts(n);
    TIC
    Octree octree(pts, 200); 
    TOC("Octree assembly");
    //computing the bounding box = 0.25s
    //computing morton codes = 0.75s
    //sorting = 8sec
    //building tree = 0.05s
}
