#include "UnitTest++.h"
#include "richardson.h"

using namespace tbem;

void richardson_test(double ratio) {
    // For n values in a sequence, Richardson extrapolation should be 
    // exact for polynomial sequences with degree n - 1.
    int n = 4;
    double allowed_error = 1e-6;
    std::vector<double> x(n);
    for (int i = 0; i < n; i++) {
        x[i] = std::pow(ratio, -i);
    }

    for (int i = 1; i < n; i++) {
        std::vector<double> input;
        for (int j = 0; j < n; j++) {
            input.push_back(std::pow(x[j], i) - 1.0);
        }
        double result = richardson_limit(ratio, input);
        CHECK_CLOSE(result, -1.0, allowed_error);
    }
}

TEST(RichardsonExtrapolate) {
    richardson_test(2.0);
    richardson_test(10.0);
}

int main() {
    return UnitTest::RunAllTests();
}
