#include <iostream>
#include "gregvct.hpp"

template <typename T>
class A {
protected:
    int x{};
public:
    A() = default;
    A(int _x) : x{_x} {}
    template <typename ...Args>
    int get(Args ...args) {
        if constexpr (!sizeof...(args))
            return 0;
        return (x + (args + ...));
    }
    void print(const char *str) const noexcept {
        std::cout << str << ", " << x << std::endl;
    }
};

template <typename T>
class B : public A<T> {
public:
    B() = default;
    B(int _x) : A<T>{_x} {}
    using A<T>::get;
    using A<T>::print;
    int get(const std::initializer_list<T> &list) {
        return A<T>::x;
    }
    void print(int a) {
        std::cout << a << ", " << A<T>::x << std::endl;
    }
};

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
                                    {1, 1, 5, 2, 8, 2, 1, 7, 3},
                                    {8, 5, 2, 3, 1, 2, 8, 9, 4},
                                    {1, 2, 3, 4, 5, 6, 1, 7, 0},
                                    {5, 8, 2, 4, 1, 9, 6, 0, 1}};
    // std::cout << mat << "\nx\n" << vec << "\n=\n" << vec.apply(mat) << "\n--------\n" << std::endl;
    std::cout << mat << "\nx\n" << vec2 << "\n=\n" << vec2.apply(mat) << "\n--------\n" << std::endl;
    std::cout << vec2.transpose() << "\nx\n" << mat << "\n=\n" << vec2*mat << "\n--------\n" << std::endl;
    std::cout << "vec: " << vec << "\n----------\n" << std::endl;
    std::cout << vec.transpose() << "\nx\n" << vec2 << "\n=\n" << vec*vec2 << "\n-----------\n" << std::endl;
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
    std::cout << mat.at(1, 0) << std::endl;
    std::cout << vec.at(1, 0) << std::endl;
    gml::vector<long double> cvec = {{4},
                                     {3},
                                     {9},
                                     {9},
                                     {7},
                                     {2},
                                     {1}};
    std::cout << "Column vector:\n" << cvec << std::endl;
    std::cout << "dotted: " << gml::dot(vec, vec2) << std::endl;
    B<int> b;
    // std::cout << b.get(2, 3) << std::endl;
    b.print(4);
    b.print("hi");
    return 0;
}
