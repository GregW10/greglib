#ifndef GREGMMAPPER_HPP
#define GREGMMAPPER_HPP

#if !defined(__unix__) && !defined(__APPLE__)
#error "Can only compiled on UNIX systems.\n"
#endif

#include <stdexcept>
#include <unistd.h>
#include <sys/mman.h>

namespace gtd {
    class unsuccessful_mmap : public std::runtime_error {
    public:
        unsuccessful_mmap() : std::runtime_error{"Error: mmap() call unsuccessful.\n"} {}
        explicit unsuccessful_mmap(const char *msg) : std::runtime_error{msg} {}
    };
    class unsuccessful_munmap : public std::runtime_error {
    public:
        unsuccessful_munmap() : std::runtime_error{"Error: munmap() call unsuccessful.\n"} {}
        explicit unsuccessful_munmap(const char *msg) : std::runtime_error{msg} {}
    };
    class mmapper {
        /* Simple class used for allocating and automatically deallocating memory using `mmap`. */
        void *_mem{};
        size_t _size{};
        const long _psize = sysconf(_SC_PAGESIZE);
        size_t _npages{}; // total number of mem. pages allocated
        size_t _esize{}; // effective number of bytes allocated
    public:
        mmapper() {
            if (_psize == -1)
                throw std::runtime_error{"Error: could not obtain page size.\n"};
        }
        explicit mmapper(size_t _s) : _size{_s}, _npages{_size/_psize + (_size % _psize != 0)}, _esize{_psize*_npages} {
            if (_psize == -1)
                throw std::runtime_error{"Error: could not obtain page size.\n"};
            if (!_s)
                return;
            this->_mem = mmap(nullptr, this->_esize, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
            if (this->_mem == MAP_FAILED)
                throw unsuccessful_mmap{};
        }
        size_t zero() noexcept {
            if (!this->_esize)
                return 0;
            char *ptr = static_cast<char*>(this->_mem);
            size_t counter = this->_esize;
            while (counter --> 0)
                *ptr++ = 0;
            return this->_esize;
        }
        size_t capacity() const noexcept {
            return this->_esize;
        }
        size_t num_pages() const noexcept {
            return this->_npages;
        }
        void *get() const noexcept {
            return this->_mem;
        }
        void *release() noexcept { // releases pointer to mapped memory without unmapping it
            this->_size = 0;
            this->_esize = 0;
            this->_npages = 0;
            void *_ptr = this->_mem;
            this->_mem = nullptr;
            return _ptr;
        }
        void *reset(size_t _s = 0) {
            if (!_s) {
                // this->~mmapper();
                if (this->_mem)
                    if (munmap(this->_mem, this->_esize) == -1)
                        throw unsuccessful_munmap{};
                this->_size = 0;
                this->_npages = 0;
                this->_esize = 0;
                return (this->_mem = nullptr);
            }
            size_t _prev_npages = this->_npages;
            this->_size = _s;
            this->_npages = _size/_psize + (_size % _psize != 0);
            if (this->_npages == _prev_npages) {
                // this->_esize = _psize*_npages;
                return this->_mem; // no need to unmap and re-map if the number of pages would be identical
            }
            // this->~mmapper();
            if (this->_mem)
                if (munmap(this->_mem, this->_esize) == -1)
                    throw unsuccessful_munmap{};
            this->_esize = _psize*_npages;
            this->_mem = mmap(nullptr, this->_esize, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
            if (this->_mem == MAP_FAILED)
                throw unsuccessful_mmap{};
            return this->_mem;
        }
        ~mmapper() {
            if (this->_mem)
                // if (munmap(this->_mem, this->_esize) == -1)
                //     throw unsuccessful_munmap{};
                munmap(this->_mem, this->_esize);
            /* I have decided not to throw an exception here if the deallocation fails, since destructors should not
             * throw, and, more importantly, given the setup of my class, none of the three reasons for `munmap` failing
             * (as listed in the man pages) should occur. */
            // this->_mem = nullptr;
        }
        operator bool() const noexcept {
            return this->_mem;
        }
    };
}
#endif
