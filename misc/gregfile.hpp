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
        int fd = -1;
        const char *fname{};
        /* char buff[BUFF_SIZE];
        char *bpos = &buff[0];
        char *wpos = bpos; */
    public:
        constexpr static uint64_t read_all = (uint64_t) -1;
        file() = default;
        file(const char *path, int flags = O_CREAT | O_RDWR, mode_t mode = S_IRUSR | S_IWUSR) noexcept
            : fd{::open(path, flags, mode)}, fname{path} {
            //if (this->fd == -1)
                // throw std::ios_base::failure{"Error: could not open file.\n"};
        }
        bool good() const noexcept {
            return this->fd > -1;
        }
        bool bad() const noexcept {
            return this->fd == -1;
        }
        bool write(const char *data, uint64_t size) {
            if (!data || !size || this->fd == -1)
                return false;
            uint64_t rem = size;
            uint64_t written;
            do {
                if ((written = ::write(this->fd, data, rem)) == -1)
                    return false;
                rem -= written;
                data += written;
            } while (rem);
            return true;
        }
        bool read(char *out, uint64_t size = read_all) {
            if (!out || !size || this->fd == -1)
                return false;
            uint64_t readd;
            if (size == read_all) {
                while (1) {
                    if ((readd = ::read(this->fd, out, 8192)) == -1)
                        return false;
                    // std::cout << "between" << std::endl;
                    if (!readd) // means EOF was reached
                        return true;
                    out += readd;
                }
            }
            uint64_t rem = size;
            do {
                if ((readd = ::read(this->fd, out, rem)) == -1)
                    return false;
                rem -= readd;
                out += readd;
            } while (rem);
            return true;
        }
        void close() noexcept {
            if (this->fd > -1) {
                ::close(this->fd);
                this->fd = -1;
            }
        }
        ~file() {
            if (this->fd > -1)
                ::close(this->fd);
        }
    };
}

#undef BUFF_SIZE
#endif
