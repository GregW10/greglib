#include "../misc/gregfile.hpp"
#include "gregsha256.hpp"
#include <sys/stat.h>
#include <vector>

/* A quick-and-dirty SHA-256 command-line program, mainly for testing my implementation. */

#define SBUFF       65'536 // stack buff size
#define HBUFF1   1'048'576 // heap  buff size 1
#define HBUFF2  16'777'216 // heap  buff size 2
#define HBUFF3 268'435'456 // heap  buff size 3

uint64_t buff_size = SBUFF;
char *buffer;


#define MAXVSIZE 268'435'456 // max. size of vector allowed

std::vector<char> v;

template <bool binary>
void hash_file(const char *fpath) {
    if (!fpath)
        return;
    gtd::file f{fpath, O_RDONLY};
    if (!f.good()) {
        fprintf(stderr, "Error: could not open \"%s\".\n", fpath);
        perror("Cause");
        return;
    }
    struct stat buff{};
    if (stat(fpath, &buff) == -1) {
        fprintf(stderr, "Error: could not obtain \"%s\" file info.\n", fpath);
        f.close();
        return;
    }
    if (buff.st_size <= HBUFF3) {
        if (buff.st_size > buff_size) {
            if (buff.st_size <= HBUFF1) {
                buffer = new char[HBUFF1];
                buff_size = HBUFF1;
                //std::cout << "Not in here?\n";
                printf("Hex: %p\n", buffer);
            } else if (buff.st_size <= HBUFF2) {
                if (buff_size == HBUFF1)
                    delete [] buffer;
                buffer = new char[HBUFF2];
                buff_size = HBUFF2;
                //std::cout << "NOOOO" << std::endl;
            } else {
                if (buff_size == HBUFF1 || buff_size == HBUFF2)
                    delete [] buffer;
                buffer = new char[HBUFF3];
                buff_size = HBUFF3;
            }
        }
        if (f.read(buffer, buff.st_size) == -1) {
            fprintf(stderr, "Error: could not read from \"%s\".\n", fpath);
            perror("Reason");
            f.close();
            return;
        }
        gtd::sha256 sha(buffer, buff.st_size, true);
        sha.hash();
        if constexpr (binary)
            std::cout.write(sha.digest(), 32);
        else
            std::cout << sha << "  " << fpath << '\n';
    } else {
        if (buff_size != HBUFF3) {
            if (buff_size != SBUFF)
                delete [] buffer;
            buffer = new char[HBUFF3];
            buff_size = HBUFF3;
        }
        int res;
        gtd::sha256 sha;
        while (1) {
            if ((res = f.read(buffer, HBUFF3)) == -1) {
                fprintf(stderr, "Error: could not read from \"%s\".\n", fpath);
                return;
            }
            sha.add_data(buffer, res);
            if (res != HBUFF3)
                break;
        }
        sha.hash();
        if constexpr (binary)
            std::cout.write(sha.digest(), 32);
        else
            std::cout << sha << "  " << fpath << '\n';
    }
    f.close();
}

int main(int argc, char **argv) {
    char stack[SBUFF];
    buffer = stack;
    bool binary = false;
    if (argc > 1) {
        char **aptr = argv + 1;
        if (**(argv + 1) == '-' && *(*(argv + 1) + 1) == 'b') {
            binary = true;
            if (argc == 2)
                goto from_stdin;
            ++aptr;
            --argc;
        }
        --argc;
        if (binary)
            while (argc --> 0)
                hash_file<true>(*aptr++);
        else
            while (argc --> 0)
                hash_file<false>(*aptr++);
        if (buff_size > SBUFF)
            delete [] buffer;
        return 0;
    }
    from_stdin:
    gtd::file f(STDIN_FILENO, true);
    ssize_t res = f.read(buffer, SBUFF);
    if (res == -1) {
        errbye:
        fprintf(stderr, "Error: could not read from stdin.\n");
        if (buff_size > SBUFF)
            delete [] buffer;
        return 1;
    }
    gtd::sha256 sha(buffer, res);
    //std::cout << "First res: " << res << std::endl;
    if (res < SBUFF)
        goto bye;
    buffer = new char[HBUFF1];
    buff_size = HBUFF1;
    if ((res = f.read(buffer, HBUFF1)) == -1)
        goto errbye;
    //std::cout << "res: " << res << std::endl;
    sha.add_data(buffer, res);
    if (res < HBUFF1)
        goto bye;
    //------------------------
    delete [] buffer;
    buffer = new char[HBUFF2];
    buff_size = HBUFF2;
    if ((res = f.read(buffer, HBUFF2)) == -1)
        goto errbye;
    sha.add_data(buffer, res);
    if (res < HBUFF2)
        goto bye;
    //------------------------
    delete [] buffer;
    buffer = new char[HBUFF3];
    buff_size = HBUFF3;
    if ((res = f.read(buffer, HBUFF3)) == -1)
        goto errbye;
    sha.add_data(buffer, res);
    if (res < HBUFF3)
        goto bye;
    while (1) {
        if ((res = f.read(buffer, HBUFF3)) == -1) {
            fprintf(stderr, "Error: could not read from stdin\n");
            return 1;
        }
        sha.add_data(buffer, res);
        if (res != HBUFF3)
            break;
    }
    bye:
    sha.hash();
    if (binary)
        std::cout.write(sha.digest(), 32);
    else
        std::cout << sha << "  -" << '\n';
    if (buff_size > SBUFF)
        delete [] buffer;
    return 0;
}
