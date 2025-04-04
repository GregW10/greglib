#include "gregfile.hpp"
#include "gregrand.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>

int main() {
	gtd::file f{"out.txt"};
	std::cout << std::boolalpha << f.write("hi", 2) << std::endl;;
	f.close();
    gtd::file r{"out.txt", O_RDONLY};
    char buff[8192]{};
    std::cout << std::boolalpha << r.read(buff) << std::endl;;
    perror("Error: ");
    r.close();
    std::cout << buff << std::endl;
    gtd::random<uint64_t> rr{};
    uint64_t i = 1'000'000;
    std::ofstream out{"rand.dat", std::ios_base::out | std::ios_base::trunc | std::ios_base::binary};
    if (!out.good())
        return 1;
    while (i --> 0) {
        uint64_t n = rr.next();
        out.write((char *) &n, sizeof(uint64_t));
    }
    constexpr uint64_t nels = 1'000'000;
    uint64_t bufff[nels];
    rr.fill(bufff, nels);
    out.write((char *) bufff, nels*sizeof(uint64_t));
    out.close();
	return 0;
}

