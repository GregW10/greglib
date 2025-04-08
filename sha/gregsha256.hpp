#ifndef GREGSHA256_HPP
#define GREGSHA256_HPP

#include "gregshabase.hpp"

namespace gtd {
    class sha256 { // : public shasum {
        uint64_t _size = 0; // total number of bytes of data seen so far
        bool _done = false;
        uint32_t _H[8] = {
            0x6a09e667,
            0xbb67ae85,
            0x3c6ef372,
            0xa54ff53a,
            0x510e527f,
            0x9b05688c,
            0x1f83d9ab,
            0x5be0cd19,
        };
        uint32_t _M[16];
        uint32_t _W[64];
        char *_par[64]; // partially-filled block
        char *_nbyte = (char *) _par;
        char _npos = 0;
        void complete_block(const char *dat) {
            /* FUNCTION IS WRONG AT THE MOMENT! FIX! */
            char beg = _npos / 4;
            uint32_t *mptr = _M;
            char *pptr = _par;
            while (beg --> 0) {
                *mptr++ = (*pptr << 24) + (*(pptr + 1) << 16) + (*(pptr + 2) << 8) + *(pptr + 3);
                pptr += 4;
            }
            char rem = _npos % 4;
            if (rem) {
                *pptr = 0;
                beg = 24;
                while (rem --> 0) {
                    *mptr += (*pptr++ << beg);
                    beg -= 8;
                }
                while (beg >= 0) {
                    *mptr += (*dat++ << beg);
                    beg -= 8;
                }

            }
            else {

            }
        }
        void process_data(const char *data, uint64_t size) {
            _size += size;
            uint64_t bytes;
            if ((bytes = _npos + size) < 64) {
                memcpy(_partial, data, size);
                _nbyte += size;
                _npos += size;
                return;
            }
            uint64_t to_process = bytes / 64;
            bytes %= 64;

            memcpy(_nbyte, data, size);


            _nbyte = ((char *) _nblock) + bytes;
            _npos = bytes;
        }
    public:
        sha256() = default;
        sha256(const char *data, uint64_t size) {
            if (!size)
                return;
            process_data(data, size);
        }
    };
}

#endif

