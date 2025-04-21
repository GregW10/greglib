#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <concepts>
#include <iomanip>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bit>

uint64_t to_uint(const char *str) {
    if (!str || !*str)
        return -1;
    uint64_t num = 0;
    goto start;
    while (*str) {
        num *= 10;
        start:
        if (*str < 48 || *str > 57) {
            fprintf(stderr, "Error: non-numeric character \"%c\" found.\n", *str);
            exit(1);
        }
        num += *str++ - 48;
    }
    return num;
}

void memcopy(char *dst, const char *src, uint64_t count) {
    if (!dst || !src || !count)
        return;
    while (count --> 0)
        *dst++ = *src++;
}

void output(uint64_t i) {
    double _d;
    // char _di[8];
    uint64_t _di = 0;
    _d = std::cbrt((double) i);
    _d -= std::floor(_d);
    _d += 1;
    memcopy((char *) &_di, (char *) &_d, 7);
    _di &= 0x000fffffffffffff;
    _di >>= 20;
    uint32_t cb = _di;
    _di = 0;
    _d = std::sqrt((double) i);
    _d -= std::floor(_d);
    _d += 1;
    memcopy((char *) &_di, (char *) &_d, 7);
    _di &= 0x000fffffffffffff;
    _di >>= 20;
    uint32_t sq = _di;
    /* _di[6] <<= 4;
    _di[6] &= (_di[5] >> 4);
    _di[5] <<= 4;
    _di[5] &= (_di[4] >> 4);
    _di[4] <<= 4;
    _di[4] &= (_di[3] >> 4);
    _di[3] <<= 4;
    _di[3] &= (_di[2] >> 4);
    uint32_t v = _di[3] + (_di[4] << 8) + (_di[5] << 16) + (_di[6] << 24); */
    // std::cout << i << ": " << _d << " -> " << std::hex << _di << " -> " << v << std::dec << '\n';
    // std::cout << i << ": " << std::hex << v << std::dec << '\n';
    if constexpr (std::endian::native == std::endian::little) {
        uint32_t temp;
        temp = (cb & 0x000000ff);
        temp <<= 8;
        temp += (cb & 0x0000ff00) >> 8;
        temp <<= 8;
        temp += (cb & 0x00ff0000) >> 16;
        temp <<= 8;
        temp += cb >> 24;
        cb = temp;
        temp = (sq & 0x000000ff);
        temp <<= 8;
        temp += (sq & 0x0000ff00) >> 8;
        temp <<= 8;
        temp += (sq & 0x00ff0000) >> 16;
        temp <<= 8;
        temp += sq >> 24;
        sq = temp;
    }
    printf("%08" PRIx32 " %08" PRIx32 "\n", cb, sq);
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: ./sieve <num_vals>\n");
        return 1;
    }
    uint64_t N = to_uint(*(argv + 1));
    if (N < 2)
        return 0;
    if (N == 2) {
        output(2);
        return 0;
    }
    if (N == 3) {
        output(2);
        output(3);
        return 0;
    }
    uint64_t i = 3;
    double broot = sqrt((double) N);
    uint64_t broot_u = broot;
    broot_u += broot_u % 2; // made sure the root-value is even as this avoids having to switch between odd/even sizes of "marked" in the loop
    uint64_t size = broot_u/2;//broot_u % 2 ? broot_u / 2 : broot_u / 2 - 1;
    std::vector<bool> marked(size, false);
    uint64_t multiple;
    output(2); // only even prime
    uint64_t two_i;
    uint64_t idx = 0;
    uint64_t midx;
    std::vector<uint64_t> primes;
    // std::cout << "broot: " << broot << ", broot_u: " << broot_u << ", marked.size(): " << marked.size() << std::endl;
    static_assert(sizeof(double) == 8);
    do {
        multiple = i*i;
        two_i = i << 1;
        midx = multiple/2 - 1;
        while (multiple <= broot_u) {
            marked[midx] = true; // marked as composite
            multiple += two_i;
            midx += i;
        }
        // std::cout << i << '\n';
        output(i);
        primes.push_back(i);
        i += 2;
        while (++idx < size && marked[idx]) i += 2;
    } while (i <= broot_u);
    // primes.push_back(i);
    uint64_t prev = broot_u + 1;
    uint64_t next = broot_u << 1;
    uint64_t maxw = broot_u*broot_u;
    maxw = maxw < N ? maxw : N;
    // std::vector<bool> window(broot_u);
    // uint64_t psq;
    // uint64_t diff;
    uint64_t p2;
    // for (const auto &prime : primes)
        // std::cout << "prime: " << prime << std::endl;
    uint64_t rem;
    // std::cout << "marked.size(): " << marked.size() << std::endl;
    while (next <= maxw) {
        idx = 0;
        while (idx < size) {
            marked[idx++] = false;
        }
        for (const uint64_t &prime : primes) {
            i = prev + ((rem = prev % prime) ? prime - rem : 0);
            if (!(i % 2))
                i += prime;
            idx = (i - prev)/2;
            p2 = prime << 1;
            while (i <= next) {
                /* if (idx == 49)
                    std::cout << "Yes! idx: " << idx << ", i: " << i << ", prime: " << prime << std::endl; */
                marked[idx] = true;
                idx += prime;
                i += p2;
            }
        }
        i = prev + 1 - (prev % 2);
        std::for_each(marked.begin(), marked.end(), [&i](bool b){if (!b) output(i); i += 2;});
        /* for (bool b : marked) {
            std::cout << std::boolalpha << "b: " << b << ", i: " << i << std::endl;
            i += 2;
        } */
        prev = next + 1;
        next += broot_u;
    }
    if (prev > N)
        return 0;
    uint64_t left = N - prev;
    left = left/2 + ((prev % 2) || (N % 2));
    // std::cout << "left: " << left << ", N: " << N << ", prev: " << prev << std::endl;
    if (left) {
        marked.resize(left);
        idx = left;
        while (idx --> 0)
            marked[idx] = false;
        for (const uint64_t &prime : primes) {
            i = prev + ((rem = prev % prime) ? prime - rem : 0);
            if (!(i % 2))
                i += prime;
            idx = (i - prev)/2;
            p2 = prime << 1;
            while (i <= N) {
                marked[idx] = true;
                idx += prime;
                i += p2;
            }
        }
        i = prev % 2 ? prev : prev + 1;
        std::for_each(marked.begin(), marked.end(), [&i](bool b){if (!b) output(i); i += 2;});
    }
    return 0;
}
