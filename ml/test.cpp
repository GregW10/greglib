#include <iostream>
#include "gregvct.hpp"

int main() {
    gml::vector<long double> vec = {{1, 2, 3, 4, 5, 6}};
    std::cout << "vec:\n" << vec << std::endl;
    gml::vector<long double> vec2 = {{8},
                                     {1},
                                     {3},
                                     {9},
                                     {3},
                                     {2},
                                     {3},
                                     {7},
                                     {4}};
    std::cout << "vec2:\n" << vec2 << "\n---------\n" << std::endl;
    gml::matrix<long double> mat = {{5, 3, 2, 5, 3, 3, 1, 7, 3},
                                    {9, 8, 7, 6, 5, 4, 8, 2, 5},
                                    {1, 12898, 5, 2, 8, 2, 1, 7, 3},
                                    {8, 5, 2, 3, 1, 2, 8, 9, 4},
                                    {1, 2, 3, 4, 5, 6, 1, 7, -193829},
                                    {5, 8, 2, 4, 1, 9, 6, 0, 1}};
    long double accum = 0;
    auto sum = vec.fold(0, +[](long double &accum, const long double &val){accum += val;});
    auto prod = vec.fold(1, +[](const long double &accum, const long double &val){return accum*val;});
    std::cout << "Sum of\n" << vec << "\nis: " << sum << std::endl;
    std::cout << "Prod. of\n" << vec << "\nis: " << prod << std::endl;
    gml::vector<int> vint;
    auto s = vint.fold(0, +[](int &accum, const int &val){accum += val;});
    auto p = vint.fold(1, +[](const int &accum, const int &val){return accum + val;});
    std::cout << "s: " << s << std::endl;
    std::cout << "p: " << p << std::endl;
    std::cout << "max of \n" << mat << "\n= " << mat.max() << std::endl;
    std::cout << "min of \n" << mat << "\n= " << mat.min() << std::endl;
    // auto x = vint.max();
    auto y = vint.min();
    return 0;
}
