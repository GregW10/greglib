#include <fstream>
#include <iostream>
#include "gregsha256.hpp"

int main(int argc, char **argv) {
    if (argc != 2)
        return 1;
    uint64_t len = strlen(*(argv + 1));
    gtd::sha256 sha(*(argv + 1), len);
    char *rev = new char[len + 1];
    char *ptr = *(argv + 1) + len;
    uint64_t counter = len;
    while (counter --> 0)
        *rev++ = *--ptr;
    *rev = 0;
    rev -= len;
    sha.add_data(rev, len);
    delete [] rev;
    const char *digest = sha.digest();
    gtd::print_hex(digest, 32);
    return 0;
}
