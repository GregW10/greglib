#ifndef GREGRAND_HPP
#define GREGRAND_HPP

#include "gregfile.hpp"
#include "gregmisc.hpp"
#include "gregmmapper.hpp"

namespace gtd {
    class invalid_buffer_size : std::logic_error {
    public:
        invalid_buffer_size() : std::logic_error{"Error: invalid buffer size.\n"} {}
        invalid_buffer_size(const char *msg) : std::logic_error{msg} {}
    };
    template <std::integral T = unsigned char>
    class random {
        const char *_fname;
        gtd::file _f;
        uint64_t _bnels = dbuffs; // number of elements of T in buffer, not buffer size
        uint64_t _bs = sizeof(T)*_bnels;
        mmapper mapper{_bs};
        T *_buff = (T*) mapper.get();
        T *_end = _buff + _bnels;
        T *_next = _end;
    public:
        constexpr static uint64_t dbuffs = 16384;
        random(uint64_t buffs = dbuffs, bool urandom = false) : _fname{urandom ? "/dev/urandom" : "/dev/random"},
                                                                _f{_fname, O_RDONLY},
                                                                _bnels{buffs} {
            if (!buffs) {
                throw invalid_buffer_size{"Error: buffer size cannot be zero! Duhhhh....\n"};
            }
            if (_f.bad()) {
                char errbuff[40];
                strcpy_c(errbuff, "Error: could not open \"");
                strcpy_c(errbuff, _fname);
                strcpy_c(errbuff, "\".\n");
                throw std::ios_base::failure{errbuff};
            }
        }
        T next() {
            if (_next == _end) {
                if (!_f.read((char*) _buff, _bs)) {
                    char errbuff[40];
                    strcpy_c(errbuff, "Error: could not read from \"");
                    strcpy_c(errbuff, _fname);
                    strcpy_c(errbuff, "\".\n");
                    throw std::ios_base::failure{errbuff};
                }
                _next = _buff;
            }
            return *_next++;
        }
        bool fill(T *buff, uint64_t nels) noexcept {
            while (nels --> 0 && _next != _end)
                *buff++ = *_next++;
            if (!nels)
                return true;
            if (!_f.read((char*) buff, nels*sizeof(T)))
                return false;
            return true;
        }
    };
}

#endif
