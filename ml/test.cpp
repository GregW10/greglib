#include <iostream>
#include <concepts>
#include <sstream>
#include <iomanip>
#include "gregvct.hpp"

int main() {
    gml::vector<long double> vec = {1, 2, 3, 4, 5, 6};
    std::cout << "vec:\n" << vec << std::endl;
    gml::vector<long double, false> vec2 = {1, 2, 3, 4, 5, 6};
    std::cout << "vec2:\n" << vec2 << "\n---------\n" << std::endl;
    gml::matrix<long double> mat = {{5, 3, 2, 5, 3, 3, 1, 7, 3},
                                    {9, 8, 7, 6, 5, 4, 8, 2, 5},
                                    {1, 1, 5, 2, 8, 2, 1, 7, 3},
                                    {8, 5, 2, 3, 1, 2, 8, 9, 4},
                                    {1, 2, 3, 4, 5, 6, 1, 7, 0},
                                    {5, 8, 2, 4, 1, 9, 6, 0, 1}};
    // std::cout << mat << "\nx\n" << vec << "\n=\n" << vec.apply(mat) << "\n--------\n" << std::endl;
    std::cout << vec2 << "\nx\n" << mat << "\n=\n" << vec2.apply(mat) << "\n--------\n" << std::endl;
    std::cout << vec2 << "\nx\n" << mat << "\n=\n" << mat*vec2 << "\n--------\n" << std::endl;
    std::cout << vec << "\nx\n" << vec2 << "\n=\n" << vec*vec2 << "\n-----------\n" << std::endl;
    gml::matrix<long double> m;
    m = std::move(mat);
    std::cout << mat.shape() << std::endl;
    gml::tensor<long double> tens = m;
    std::cout << tens << std::endl;
    tens = mat;
    std::cout << tens << std::endl;
    tens = m;
    mat = tens;
    std::cout << mat << std::endl;
    mat = vec;
    std::cout << mat << std::endl;
    vec = mat;
    std::cout << vec << std::endl;
    return 0;
}
