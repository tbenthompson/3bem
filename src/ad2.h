#ifndef QWEKJHWELR_AD2_H
#define QWEKJHWELR_AD2_H
#include <vector>
#include <functional>
#include <cmath>
#include <iostream>
#include <cassert>
#include <memory>

// template <typename T>
// std::unique_ptr<T> maybe(std::unique_ptr<T>& t) {
//     if (t == nullptr) {
//         return nullptr;
//     }
//     return std::move(t);
// }
// 
// template <typename T>
// struct D {
//     typedef D<T> type;
//     typedef std::unique_ptr<D<T>> type_ptr;
// 
//     D(D<T>& c): v(c.v), d(maybe(c.d)) { }
//     D(const T& x): v(x), d(nullptr) {}
//     D(const T& x, std::unique_ptr<type>& p_d): v(x), d(maybe(p_d)) {}
// 
//     const T v;
//     const std::unique_ptr<type> d;
// // 
// //     type chain(std::function<double(double)> f,
// //                std::function<type(const type&)> fp) const {
// //         std::array<T,M> res{};    
// //         res[0] = f(v[0]);
// //         if (degree == 0) {
// //             return type{res, 0};
// //         }
// //         auto c = fp(init()) * tail();
// //         for (int i = 0; i < c.degree; i++) {
// //             res[i + 1] += c.v[i];
// //         }
// //         return type{res, degree};
// //     }
// // 
//     type operator*(const type& b) const {
//         if(d == nullptr && b.d == nullptr) {
//             return type{v * b.v, nullptr};
//         } else if (b.d == nullptr) {
//             return type{v * b.v, d};
//         } else if (d == nullptr) {
//             return type{v * b.v, b.d};
//         }
//         return type{v * b.v, type_ptr(new type((*this) * *b.d + *d * b))};
//     }
// 
//     type operator+(const type& b) const {
//         if (d == nullptr && b.d == nullptr) {
//             return type{v + b.v, nullptr};
//         } else if (b.d == nullptr) {
//             return type{v + b.v, d};
//         } else if (d == nullptr) {
//             return type{v + b.v, b.d};
//         }
//         return type{v + b.v, type_ptr(new type(*d + *b.d))};
//     }
// 
//     type operator-() const {
//         if (d == nullptr) {
//             return type{-v, nullptr};
//         }
//         return type {-v, type_ptr(new type(-*d))};
//     }
// 
//     type operator-(const type& b) const {
//         return (*this) + (-b);
//     }
// 
//     std::vector<double> get() {
//         std::vector<double> res{v};
//         if (d == nullptr) {
//             return res;
//         }
//         auto tail = d->get();
//         res.insert(res.end(), tail.begin(), tail.end());
//         return res;
//     }
// 
//     friend std::ostream& operator<<(std::ostream& os, const type& a) {
//         os << "(" << a.v;
//         if (a.d != nullptr) {
//             os << ", " << *a.d;
//         }
//         os << ")";
//         return os;
//     }
// };
// 
// template <typename T>
// D<T> constant(const T& x) {
//     return D<T>(x);
// }
// 
// template <typename T>
// D<T> var(const T& x) {
//     return D<T>(x, std::unique_ptr<D<T>>(new D<T>(constant(1.0))));
// }
// 
// // 
// // 
// // template <typename T, int M>
// // D<T,M> exp(const D<T,M>& a) {
// //     return a.chain((double(*)(double))&std::exp, &exp<T,M>);
// // }
// // 
// // template <typename T, int M>
// // D<T,M> sin(const D<T,M>& a);
// // 
// // template <typename T, int M>
// // D<T,M> cos(const D<T,M>& a) {
// //     return a.chain((double(*)(double))&std::cos, [](const D<T,M>& x) {
// //         return -sin<T,M>(x);
// //     });
// // }
// // 
// // template <typename T, int M>
// // D<T,M> sin(const D<T,M>& a) {
// //     return a.chain((double(*)(double))&std::sin, &cos<T,M>);
// // }
// // 
// // template <typename T>
// // T recip(const T& x) {return 1.0 / x;}
// // 
// // template <typename T, int M>
// // D<T,M> recip(const D<T,M>& a) {
// //     return a.chain(recip<T>, [](const D<T,M>& x) {
// //         auto r = recip(x); 
// //         return (-r) * r;
// //     });
// // }
// // 
// // template <typename T, int M>
// // D<T,M> log(const D<T,M>& a) {
// //     return a.chain((double(*)(double))&std::log, &recip<T,M>);
// // }
// // 
// // template <typename T, int M>
// // D<T,M> sqrt(const D<T,M>& a) {
// //     return a.chain((double(*)(double))&std::sqrt, [](const D<T,M>& x) {
// //         return recip(sqrt(x) * 2.0);
// //     });
// // }
// // 
// // template <typename T, int M>
// // D<T,M> asin(const D<T,M>& a) {
// //     return a.chain((double(*)(double))&std::asin, [](const D<T,M>& x) {
// //         return recip(sqrt(-(x * x - 1.0)));
// //     });
// // }
// // 
// 
// typedef D<double> Dd;
// 
#endif
