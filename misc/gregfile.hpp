#ifndef GREGFILE_HPP
#define GREGFILE_HPP

#include <fcntl.h>
#include <cstdint>
#include <unistd.h>
#include <cinttypes>

#define BUFF_SIZE 8192

namespace gtd {
    class file {
        int fd = -1;
        const char *fname{};
        char buff[BUFF_SIZE];
        char *bpos = &buff[0];
        char *wpos = bpos;
    public:
        file() = default;
        file(const char *path, int flags = O_CREAT, mode_t mode = S_IRUSR | S_IWUSR) : fd{open(path, flags, mode)}, fname{path} {
            //if (this->fd == -1)
                // throw std::ios_base::failure{"Error: could not open file.\n"};
        }
        bool good() const noexcept {
            return this->fd > -1;
        }
        bool bad() const noexcept {
            return this->fd == -1;
        }
        bool write_data(const char *data, uint64_t size) {
            if (!data || !size || this->fd == -1)
                return false;
            uint64_t rem = size;
            uint64_t written;
            do {
                if ((written = write(this->fd, data, rem)) == -1)
                    return false;
                rem -= written;
                data += written;
            } while (rem);
            return true;
        }
        void close_file() {
            if (this->fd > -1) {
                close(this->fd);
                this->fd = -1;
            }
        }
        ~file() {
            if (this->fd > -1)
                close(this->fd);
        }
    };


}

#undef BUFF_SIZE

#endif
