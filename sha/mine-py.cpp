#include "gregsha256.hpp"
#include <csignal>
#include <fstream>
#include <sys/stat.h>
#include <mutex>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>

#define NONCE_SIZE 24

bool bye = false;

std::mutex mtx;

gtd::sha256 low("Highest hash to begin with! (Probably)", 38, true);
uint64_t    val = 0;

void handler(int sig) {
    bye = true;
}

void worker(const unsigned char *data, uint64_t dsize, uint64_t tid, uint64_t interval) {
    std::cout << dsize << ", " << tid << ", " << interval << std::endl;
    gtd::sha256 local_low("local", 5, true);
    uint64_t local_best_val = 0;
    gtd::sha256 base(data, dsize, false);
    gtd::sha256 sha;
    char buff[NONCE_SIZE]{'\n'};
    // buff[0] = '\n';
    buff[1] = '#';
    buff[2] = ' ';
    uint64_t _i = tid;
    char *ptr = ((char*) &buff) + 3;
    while (!bye) {
        snprintf(ptr, 21, "%020x", _i);
        buff[23] = '\n';
        // sha.add_data(data, dsize);
        sha = base;
        sha.add_data(buff, NONCE_SIZE);
        sha.hash();
        // std::cout << sha << std::endl;
        if (sha < local_low) {
            local_low = sha;
            local_best_val = _i;
            std::cout << sha << " at " << _i << "\n";
        }
        sha.reset();
        _i += interval;
    }
    std::lock_guard<std::mutex> guard{mtx};
    if (local_low < low) {
        low = local_low;
        val = local_best_val;
    }
}

uint64_t to_uint(const char *str) {
    if (!str)
        return -1;
    uint64_t res = 0;
    while (*str) {
        if (*str < 48 || *str > 57) {
            std::cerr << "Error: non-numeric character '" << *str << "' found.\n";
            exit(1);
        }
        res *= 10;
        res += *str++ - 48;
    }
    return res;
}


int main(int argc, char **argv) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: ./" << *argv << " <in_fname> <out_fname> [num_threads]\n";
        return 1;
    }
    uint64_t hconc = std::thread::hardware_concurrency();
    uint64_t nthreads;
    if (argc == 4) {
        nthreads = to_uint(*(argv + 3));
        if (!nthreads || nthreads > hconc)
            nthreads = 1;
    } else
        nthreads = hconc;
    signal(SIGINT , &handler);
    signal(SIGTERM, &handler);
    struct stat buff{};
    if (stat(*(argv + 1), &buff) == -1) {
        std::cerr << "Error: couldn't obtain file information for \"" << *(argv + 1) << "\".\n";
        return 1;
    }
    if (!buff.st_size) {
        std::cerr << "File \"" << *(argv + 1) << "\" is empty. Nothing to do.\n";
        return 2;
    }
    std::ifstream in{*(argv + 1), std::ios_base::in | std::ios_base::binary};
    if (!in.good()) {
        std::cerr << "Error: couldn't open \"" << *(argv + 1) << "\".\n";
        return 1;
    }
    const unsigned char *data = new unsigned char[buff.st_size + NONCE_SIZE];
    in.read((char *) data, buff.st_size);
    in.close();
    uint64_t tid = 0;
    std::vector<std::thread> threads;
    std::cout << low << std::endl;
    while (tid < nthreads)
        threads.emplace_back(worker, data, buff.st_size, tid++, nthreads);
    for (std::thread &t : threads)
        t.join();
    std::cout << "Best reached at: " << val << std::endl;
    std::ofstream out{*(argv + 2), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary};
    if (!out.good())
        return 1;
    char buf[NONCE_SIZE]{'\n'};
    buf[1] = '#';
    buf[2] = ' ';
    snprintf(&buf[3], 20, "%020x", val);
    buf[23] = '\n';
    out.write((char *) data, buff.st_size);
    out.write((char *) buf, NONCE_SIZE);
    out.close();
    return 0;
}
