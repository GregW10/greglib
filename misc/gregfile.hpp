#ifndef GREGFILE_HPP
#define GREGFILE_HPP

#include <fcntl.h>
#include <cstdint>
#include <iostream>
#include <unistd.h>
#include <cinttypes>

#define BUFF_SIZE 8192

namespace gtd {
    class file {
        int _fd = -1;
        const char *_fname{};
        bool _dont_close = false;
        /* char buff[BUFF_SIZE];
        char *bpos = &buff[0];
        char *wpos = bpos; */
    public:
        constexpr static size_t read_all = (size_t) -1;
        file() = default;
        file(const char *path, int flags = O_CREAT | O_RDWR, mode_t mode = S_IRUSR | S_IWUSR) noexcept
            : _fd{::open(path, flags, mode)}, _fname{path} {
            //if (this->_fd == -1)
                // throw std::ios_base::failure{"Error: could not open file.\n"};
        }
        file(int fd, bool dont_close = true) : _fd{fd}, _dont_close{dont_close} {}
        bool good() const noexcept {
            return this->_fd > -1;
        }
        bool bad() const noexcept {
            return this->_fd == -1;
        }
        bool write(const char *data, uint64_t size) {
            if (!data || !size || this->_fd == -1)
                return false;
            uint64_t rem = size;
            uint64_t written;
            do {
                if ((written = ::write(this->_fd, data, rem)) == -1)
                    return false;
                rem -= written;
                data += written;
            } while (rem);
            return true;
        }
        ssize_t read(char *out, size_t size = read_all) {
            if (!out || !size || this->_fd == -1)
                return 0;
            ssize_t tot = 0;
            ssize_t readd;
            if (size == read_all) {
                while (1) {
                    if ((readd = ::read(this->_fd, out, read_all)) == -1)
                        return -1;
                    // std::cout << "between" << std::endl;
                    if (!readd) // means EOF was reached
                        return tot;
                    out += readd;
                    tot += readd;
                }
            }
            size_t rem = size;
            do {
                if ((readd = ::read(this->_fd, out, rem)) == -1)
                    return -1;
                if (!readd)
                    return tot;
                rem -= readd;
                out += readd;
                tot += readd;
            } while (rem);
            return tot;
        }
        void close() noexcept { // ignores "_dont_close"
            if (this->_fd > -1) {
                ::close(this->_fd);
                this->_fd = -1;
            }
        }
        ~file() {
            if (this->_fd > -1 && !this->_dont_close)
                ::close(this->_fd);
        }
    };
}

#undef BUFF_SIZE
#endif
