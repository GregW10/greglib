#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cinttypes>
#include "../misc/gregparse.hpp"
#include "greggen.hpp"
#include <vector>

uint64_t to_uint64(const char *str) {
    if (!str)
        exit(1);
    uint64_t tot = 0;
    goto start;
    while (*str) {
        tot *= 10;
        start:
        if (*str < 48 || *str > 57) {
            fprintf(stderr, "\033[1m\033[31mError: \033[32mcharacter '%c' is non-numeric.\n", *str);
            exit(1);
        }
        tot += *str++ - 48;
    }
    return tot;
}

bool equal(const char *s1, const char *s2) {
    if (!s1 || !s2)
        return false;
    while (*s1 || *s2)
        if (*s1++ != *s2++)
            return false;
    return true;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "\033[1m\033[31mError: \033[38;5;11minvalid number of command-line arguments.\n");
        return 1;
    }
    gtd::parser parser{argc, argv};
    const char *shape = parser.get_arg(std::regex{R"(^\(\s*\d{1,18}\s*(\s*,\s*\d{1,18}\s*)*\)$)"});
    if (!shape) {
        fprintf(stderr, "\033[1m\033[31mError: \033[38;5;11mno valid shape provided.\n");
        return 1;
    }
    const char *T_size = parser.get_arg(std::regex{R"(^\d{1,18}$)"});
    if (!shape) {
        fprintf(stderr, "\033[1m\033[31mError: \033[38;5;11mno valid sizeof(T) provided.\n");
        return 1;
    }
    const char *endptr;
    uint64_t T_s = strtoull(T_size, const_cast<char**>(&endptr), 10);
    if (T_size == endptr) {
        fprintf(stderr, "\033[1m\033[31mError: \033[38;5;11mstrtoull error.\n");
        return 1;
    }
    uint64_t rank = 0;
    std::vector<uint64_t> dims;
    while (*shape) {
        if (*shape < 48 || *shape > 57) {
            ++shape;
            continue;
        }
        ++rank;
        endptr = shape;
        dims.push_back(strtoull(shape, const_cast<char**>(&endptr), 10));
        if (endptr == shape) {
            fprintf(stderr, "\033[1m\033[31mError: \033[38;5;11mstrtoull error.\n");
            return 1;
        }
        shape = endptr; // to point to one past the last char that was converted
    }
    if (rank == 2) {
        uint64_t fsize;
        printf("\033[38;5;93m.tsr file size \033[38;5;226m=\033[1m\033[38;5;10m %" PRIu64
        "\033[0m \033[38;5;14mbytes\n", (fsize = 4*(5 + 2*rank) + gml::gen::product(dims.begin(), dims.end())*T_s));
        printf("\033[38;5;93m.mtsr file size \033[38;5;226m=\033[1m\033[38;5;10m %" PRIu64
        "\033[0m \033[38;5;14mbytes\n", fsize - 8);
        return 0;
    }
    printf("\033[38;5;93m.tsr file size \033[38;5;226m=\033[1m\033[38;5;10m %" PRIu64 "\033[0m \033[38;5;14mbytes\n",
            4*(5 + 2*rank) + gml::gen::product(dims.begin(), dims.end())*T_s);
    return 0;
}
