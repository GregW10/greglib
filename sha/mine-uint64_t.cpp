#include "gregsha256.hpp"
#include <csignal>
#include <fstream>

bool bye = false;

void handler(int sig) {
    bye = true;
}


int main(int argc, char **argv) {
    signal(SIGINT , &handler);
    signal(SIGTERM, &handler);
    const uint64_t start = 1024ull*1024*1024*1024;
    uint64_t counter = start;
    gtd::sha256 sha;
    gtd::sha256 low("hello", 5, true);
    std::cout << "Starting at: " << low << std::endl;
    while (true) {
        sha.add_data(&counter, sizeof(uint64_t));
        sha.hash();
        if (sha < low) {
            low = sha;
            std::cout << low << " at " << counter << '\n';
        }
        sha.reset();
        if (bye)
            break;
        ++counter;
    }
    std::cout << "Best reached at: " << counter << std::endl;
    std::ofstream out{"best.bin", std::ios_base::out | std::ios_base::trunc | std::ios_base::binary};
    if (!out.good())
        return 1;
    out.write((char *) &counter, sizeof(uint64_t));
    out.close();
    return 0;
}
