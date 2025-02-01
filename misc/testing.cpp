#include "gregnum.hpp"
#include <iostream>

int main() {

    gtd::big_integer i(0x0fffffffffffffffull);
    gtd::big_integer j(0x0fffffffffffffffull);

    char n = 49;
/*
    printf("n = %x\n", n);

    std::cout << i << std::endl;

    unsigned char a = 255;
    unsigned char b = 255;

    unsigned char c = a + b;

    printf("%d\n", c);
*/
    std::cout << j << " + " << i << " = ";

    j += i;

    std::cout << j << std::endl;

        std::cout << j << " + " << i << " = ";

    j += i;

    std::cout << j << std::endl;

        std::cout << j << " + " << i << " = ";

    j += i;

    std::cout << j << std::endl;

        std::cout << j << " + " << i << " = ";

    j += i;

    std::cout << j << std::endl;

        std::cout << j << " + " << i << " = ";

    j += i;

    std::cout << j << std::endl;

    return 0;
}
