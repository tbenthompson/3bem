#ifndef __CHEB_H
#define __CHEB_H

inline double to_reference_coords(double x, double a, double b) {
    return 2 * ((x - a) / (b - a)) - 1;
}

std::vector<double> cheb_polys(double x_hat, int n_max) {
    std::vector<double> res(n_max + 1);
    res[0] = 1.0; 
    if (n_max == 0) {
        return res;
    }
    res[1] = x_hat; 
    if (n_max == 1) {
        return res;
    }
    for (int i = 2; i < n_max + 1; i++) {
        res[i] = 2 * x_hat * res[i - 1] - res[i - 2];
    }
    return res;
}

double s_n(double x_hat, double y_hat, int n) {
    auto x_cheb = cheb_polys(x_hat, n - 1);
    auto y_cheb = cheb_polys(y_hat, n - 1);
    double result = 0.0;
    for (int i = 1; i < n; i++) {
        result += x_cheb[i] * y_cheb[i];
    }
    result *= (2.0 / n); 
    result += (1.0 / n);
    return result;
}

std::vector<double> cheb_pts_first_kind(unsigned int n_pts) {
    std::vector<double> x(n_pts);
    for (unsigned int i = 1; i <= n_pts; i++) {
        x[i - 1] = std::cos((2 * i - 1) * M_PI / (2 * n_pts));
    }
    return x;
}

std::vector<double> cheb_pts_second_kind(unsigned int n_pts) {
    std::vector<double> x(n_pts);
    for (unsigned int i = 0; i < n_pts; i++) {
        x[i] = std::cos(i * M_PI / (n_pts - 1));
    }
    return x;
}
#endif
