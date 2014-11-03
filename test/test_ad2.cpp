#include "UnitTest++.h"
#include "ad2.h"
#include "util.h"
#include <iostream>
#include "libtaylor/taylor.hpp"

// TEST(Create) {
//     auto t = var(-1.5);
//     auto b = t + 5 + t - 3 - t + (-2);
//     CHECK_ARRAY_CLOSE(b.get(), t.get(), 2, 1e-15);
// }
// 
// TEST(Cube) {
//     //TODO: Move over some tests.
//     TIC
//     auto t = var(-1.5);
//     int n = 6;
//     std::vector<Dd> derivs = {t}; 
//     for (int i = 1; i < n; i++) {
//         derivs.push_back(derivs[i - 1] * derivs[i - 1]);
//     }
//     TOC("Mult");
//     std::cout << derivs[n - 1] << std::endl;
// }
// 
// TEST(Recip) {
//     auto t = var<double,5>(1.0);
//     auto t_r = recip(t);
//     auto t2 = recip(t_r);
//     CHECK_ARRAY_CLOSE(t.v,t2.v,5,1e-15);
// }
// 
// TEST(Log) {
//     auto t = var<double,5>(1.0);
//     auto log_t = log(t);
//     auto t2 = exp(log_t);
//     CHECK_ARRAY_CLOSE(t.v,t2.v,5,1e-15);
// }
// 
// TEST(Sqrt) {
//     auto t = var<double,5>(1.0);
//     auto sqrt_t = sqrt(t);
//     auto t2 = sqrt_t * sqrt_t;
//     CHECK_ARRAY_CLOSE(t.v,t2.v,5,1e-15);
// }
// 
// TEST(SinInv) {
//     auto t = var<double,5>(0.5);
//     auto asin_t = asin(t);
//     auto t2 = sin(asin_t);
//     CHECK_ARRAY_CLOSE(t.v,t2.v,5,1e-15);
// }

using namespace std;
template<class T>
T f(const T &x, const T &y)
{
  return sin(log(7*x)+exp(y))+9;
}

TEST(ABC) {
  // Compute a directional derivative of f(x,y) in the direction
  // (1,2), at point (x,y) = (3,4). This is equivalent to computing the
  // taylor expansion of f(3+1*eps, 4+2*eps) in the variable eps.
  const int Ndeg = 30; // Order of expansion.
  const int Nvar = 1; // Only one differentiation variable in this example
  TIC
    taylor<double,Nvar,Ndeg> fexpansion;
  for (int i = 0; i < 100000; i++) {
      taylor<double,Nvar,Ndeg> eps(0,0); // Set up seed variable.
        fexpansion = f(3+eps,4+2*eps);
      // Now fexpansion contains Taylor coefficients. If we want
      // derivatives we have to multiply by the appropriate factorials
      // fexpansion.deriv_facs();
  }
  TOC("libtaylor");
  cout << "Directional derivative: " << fexpansion << endl;
}

int main() {
    return UnitTest::RunAllTests();
}
