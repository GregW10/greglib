#include "gregnum.hpp"
#include <iostream>

int main() {

    gtd::big_integer i(0x0fffffffffffffffull);
    gtd::big_integer j(0x0fffffffffffffffull);

    std::cout << j << std::endl;

    j += 1;

    std::cout << j << std::endl;

    j.shift_left();

    std::cout << j << std::endl;

        j.shift_left();

    std::cout << j << std::endl;

        j.shift_left();

    std::cout << j << std::endl;

        j.shift_left();

    std::cout << j << std::endl;

        j.shift_left();

    std::cout << j << std::endl;

    char n = 49;
/*
    printf("n = %x\n", n);

    std::cout << i << std::endl;

    unsigned char a = 255;
    unsigned char b = 255;

    unsigned char c = a + b;

    printf("%d\n", c);
*/
    unsigned int c = 1'000;

    while (c --> 0) {
        std::cout << c << std::flush << "\r";
        j += i;

    }

    std::cout << std::endl;

    std::cout << j << std::endl;

    uint64_t x = 2;
    uint64_t y = 3;

    uint64_t res = x - y;

    std::cout << "What the fuck" << std::endl;

    std::cout << std::dec << "res: " << res << std::endl;

    x = 0b11111111;
    y = 0b1000000000;

    res = y - x;

    std::cout << res << std::endl;

    return 0;
}
