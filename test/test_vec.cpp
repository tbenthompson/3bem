#include "UnitTest++.h"
#include "vec.h"
#include "numerics.h"
#include "util.h"

struct Data {
    Vec3<double> a = {{1.0, 1.0, 2.0}};
    Vec3<double> b = {{2.0, 0.5, -1.0}};
};

TEST_FIXTURE(Data, VecAdd) {
    auto c = a + b; 
    b += a;
    double exact[3] = {3.0, 1.5, 1.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-12);
    CHECK_ARRAY_CLOSE(b, exact, 3, 1e-12);
}

TEST_FIXTURE(Data, VecSub) {
    auto c = a - b; 
    double exact[3] = {-1.0, 0.5, 3.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-12);
}

TEST_FIXTURE(Data, VecMul) {
    auto c = a * b; 
    double exact[3] = {2.0, 0.5, -2.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-12);
}

TEST_FIXTURE(Data, VecDiv) {
    auto c = a / b; 
    double exact[3] = {0.5, 2.0, -2.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-12);
}

TEST_FIXTURE(Data, VecNorm) {
    auto c = normalized(b);
    normalize(b);
    double m = std::sqrt(5.25);
    double exact[3] = {2 / m, 0.5 / m, -1 / m};
    CHECK_ARRAY_CLOSE(b, exact, 3, 1e-6);
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-6);
}

TEST_FIXTURE(Data, VecNegate) {
    auto c = -b;
    double exact[3] = {-2.0, -0.5, 1.0};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-6);
}

TEST_FIXTURE(Data, VecCross) {
    auto c = cross(a, b);
    double exact[3] = {-2, 5, -1.5};
    CHECK_ARRAY_CLOSE(c, exact, 3, 1e-6);
}

TEST(VecPrint) {
    Vec3<double> a = {1.0, 2.0, 3.0};
    std::cout << a << std::endl;
}


int main() {
    return UnitTest::RunAllTests();
}
