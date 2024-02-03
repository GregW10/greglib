#include "gregbase.hpp"

int main() {
    // gml::tens_ld::precision = 30;
    gml::tens_ld t{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1957483748375.28329292138}, {2, 2, 4}};
    std::cout << t.at(1, 0, 2) << std::endl;
    t.reshape({2, 8});
    std::cout << t.at(1, 1) << std::endl;
    t.reshape({2, 4, 2});
    std::cout << t << std::endl;
    gml::tens_ld t2;
    std::cout << std::boolalpha << "t.empty(): " << t.empty() << std::endl;
    std::cout << std::boolalpha << "t2.empty(): " << t2.empty() << std::endl;
    std::cout << "t just before writing:\n" << t << std::endl;
    t.to_tsr("out.tsr");
    t2.to_tsr("empty.tsr");
    gml::tens_ld t3{"out.tsr"};
    std::cout << "t3 after reading in:\n" << t3 << std::endl;
    std::cout << "t == t3:" << (t == t3) << std::endl;
    gml::tens_i i1 = {{1, 2, 3, 4}, {4}};
    gml::tens_f i2 = {{1, 2, 3, 4}, {4}};
    std::cout << "i1 == i2: " << (i1 == i2) << std::endl;
    return 0;
}
