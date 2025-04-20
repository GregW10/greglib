#include <fstream>
#include <iostream>
#include "gregsha256.hpp"

int main(int argc, char **argv) {
    if (argc != 2)
        return 1;
    gtd::sha256 sha(*(argv + 1), strlen(*(argv + 1)));
    std::ofstream out{"shit.dat", std::ios_base::binary | std::ios_base::trunc};
    if (!out.good())
        return 2;
    const char *digest = sha.digest();
    out.write(digest, 32);
    out.close();
    return 0;
}
