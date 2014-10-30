#include "quadrature.h"
#include <iostream>
#include <cmath>

int main() {
    double b_start = 0.5;
    int start = 32; //sufficient for 1e-13 error -- for double precision
    int step = 19; //total of 727 pts
    double b_min = 0.002;

    // int start = 15; //sufficient for 1e-6 error -- for single precision
    // int step = 10; //total of 175 pts
    // double b_min = 0.03;

    // int start = 10; //sufficient for 1e-3 error -- fast
    // int step = 6; //total of 54 pts
    // double b_min = 0.1;

    int total_pts = 0;
    for (double b = b_start; b > b_min; b *= 0.5) {
        int dili_n = (int)(start - step * std::log(b));
        total_pts += dili_n;
        auto dili_map = diligenti_mapping(dili_n, 0.0, 7);
        // double y_dili = integrate(dili_map, [&] (double x) {
        //     return 1.0 / std::pow(x * x + b * b, 0.5);
        // });
        // double y_exact = -log(-1 + sqrt(1 + b * b)) + log(1 + sqrt(1 + b * b));
        double y_dili = integrate(dili_map, [&] (double x) {
            return 1.0 / std::pow(x * x + b * b, 1.5);
        });
        double y_exact = 2 / ((b * b) * sqrt(1 + b * b));
        double error = std::fabs((y_dili - y_exact) / y_exact);
        std::cout.precision(17);
        std::cout << "b = " << b << std::endl;
        std::cout << "error = " << error << std::endl;
        std::cout << "with pts = " << dili_n << std::endl;
    }
    std::cout << "\ntotal pts = " << total_pts << std::endl;
}
