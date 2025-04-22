#include <fstream>
#include <cstdlib>
#include <iostream>
#include "gregsha256.hpp"

uint64_t to_uint(const char *s) {
    if (!s)
        return -1;
    uint64_t val = 0;
    goto start;
    while (*s) {
        val *= 10;
        start:
        if (*s < 48 || *s > 57) {
            std::cerr << "Error: non-numeric character \"" << *s << "\" found.\n";
            exit(1);
        }
        val += *s++ - 48;
    }
    return val;
}

int main(int argc, char **argv) {
    if (argc != 2)
        return 1;
    uint64_t n = to_uint(*(argv + 1));
    // std::cout << n << std::endl;
    uint64_t header[10];
    std::ifstream in{"/dev/random", std::ios_base::in | std::ios_base::binary};
    if (!in.good())
        return 2;
    in.read((char *) header, 72);
    in.close();
    header[9] = 0;
    uint64_t &nonce = header[9];
    gtd::sha256 lowest((char *) header, 80);
    lowest.digest();
    uint64_t lnonce = 0;
    gtd::sha256 sha1;
    gtd::sha256 sha2;
    while (n --> 0) {
        ++nonce;
        sha1.reset();
        sha2.reset();
        sha1.add_data((char *) header, 80);
        sha2.add_data(sha1.digest(), 32);
        sha2.digest();
        if (sha2 < lowest) {
            lowest = sha2;
            lnonce = nonce;
        }
    }
    std::cout << "Lowest:\n";
    gtd::print_hex(lowest.digest(), 32);
    printf("Nonce: %" PRIu64 "\n", lnonce);
    header[9] = lnonce;
    std::ofstream out{"lowest-header.dat", std::ios_base::binary | std::ios_base::trunc};
    if (!out.good())
        return 3;
    out.write((char *) header, 80);
    out.close();
    return 0;
}
