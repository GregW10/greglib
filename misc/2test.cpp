#include "gregnum.hpp"
#include <iostream>

int main() {

    uint64_t val = 0x9f38ed82ac;

    std::cout << std::hex << val << std::dec << std::endl;

    gtd::big_integer i(0x0fffffffffffffefull);
    gtd::big_integer j(0x0fffffffffffffffull);

    uint64_t by = 123;

    std::cout << i << " * " << by << " =" << std::endl;

    i *= by;

    std::cout << i << std::endl;

    return 0;
}
