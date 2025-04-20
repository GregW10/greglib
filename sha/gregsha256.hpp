#ifndef GREGSHA256_HPP
#define GREGSHA256_HPP

#include <cstring>
#include "gregshabase.hpp"
#include "gregsha-constants.hpp"

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
        char  _par[64]; // partially-filled block
        char *_nbyte = (char *) _par;
        unsigned char _npos = 0;
        void copy_block(const char *dat) {
            uint32_t *mptr = _M;
            uint32_t *wptr = _W;
            char counter = 16;
            while (counter --> 0) {
                *mptr = (*dat << 24) + (*(dat + 1) << 16) + (*(dat + 2) << 8) + *(dat + 3);
                *wptr++ = *mptr++; // remember to define a separate function for this if BIG-ENDIAN
                dat += 4;
            }
        }
        /* void complete_block(const char *dat) {
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
        } */
        void process_block(const char *data) {
            copy_block(data);
            char counter = 48;
            uint32_t *wptr = ((uint32_t *) _W) + 16;
            while (counter --> 0) {
                *wptr = sigma1_256(*(wptr - 2)) + *(wptr - 7) + sigma0_256(*(wptr - 15)) + *(wptr - 16);
                ++wptr;
            }
            uint32_t a = _H[0];
            uint32_t b = _H[1];
            uint32_t c = _H[2];
            uint32_t d = _H[3];
            uint32_t e = _H[4];
            uint32_t f = _H[5];
            uint32_t g = _H[6];
            uint32_t h = _H[7];
            counter = 64;
            const uint32_t *kptr = sha256_constants;
            wptr = (uint32_t *) _W;
            uint32_t T1;
            uint32_t T2;
            while (counter --> 0) {
                T1 = h + Sigma1_256(e) + ch(e, f, g) + *kptr++ + *wptr++;
                T2 = Sigma0_256(a) + maj(a, b, c);
                h = g;
                g = f;
                f = e;
                e = d + T1;
                d = c;
                c = b;
                b = a;
                a = T1 + T2;
                std::cout << "va" << std::endl;
            }
            _H[0] += a;
            _H[1] += b;
            _H[2] += c;
            _H[3] += d;
            _H[4] += e;
            _H[5] += f;
            _H[6] += g;
            _H[7] += h;
                for (int i = 0; i < 8; ++i)
                    std::cout << std::hex << _H[i] << std::endl;
        }
        void process_data(const char *data, uint64_t size) {
            _size += size;
            uint64_t bytes = _npos + size;
            if (bytes < 64) {
                memcpy(_nbyte, data, size);
                _nbyte += size;
                _npos += size;
                return;
            }
            uint64_t N;
            if (bytes % 64) {
                N = 64 - _npos; // N just acting as a temp. variable
                memcpy(_nbyte, data, N);
                data += N;
                process_block(_par);
                _npos += (size % 64);
                _npos = _npos % 64; // WATCH OUT
                _nbyte = ((char *) _par) + _npos;
                N = (bytes - N)/64;
            } else
                N = bytes / 64;
            while (N --> 0) {
                process_block(data);
                data += 64;
            }
            /* _size += size;
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
            _npos = bytes; */
        }
        void finish() {
            if (_npos < 55) {
                *_nbyte++ = 0b10000000;
                ++_npos;
                bye:
                while (_npos++ < 56)
                    *_nbyte++ = 0;
                _par[56] = (_size >> 56);
                _par[57] = (_size >> 48);
                _par[58] = (_size >> 40);
                _par[59] = (_size >> 32);
                _par[60] = (_size >> 24);
                _par[61] = (_size >> 16);
                _par[62] = (_size >> 8);
                _par[63] =  _size;
                process_block(_par);
                _done = true;
                return;
            }
            *_nbyte++ = 0b10000000;
            ++_npos;
            while (_npos < 64) {
                *_nbyte++ = 0;
                ++_npos;
            }
            process_block(_par);
            _nbyte = _par;
            _npos = 0;
            goto bye;
        }
    public:
        sha256() = default;
        sha256(const char *data, uint64_t size) {
            if (!data || !size)
                return;
            process_data(data, size);
        }
        void add_data(const char *data, uint64_t size) noexcept {
            if (!data || !size || _done)
                return;
            process_data(data, size);
        }
        const char *digest() noexcept {
            static unsigned char _digest[32];
            if (_done)
                return (char *) _digest;
            finish();
            unsigned char *dptr = _digest;
            unsigned char *hptr = (unsigned char *) _H;
            char counter = 8;
            while (counter --> 0) {
                *dptr++ = (*hptr >> 24);
                *dptr++ = (*hptr >> 16);
                *dptr++ = (*hptr >>  8);
                *dptr++ =  *hptr++;
            }
            return (char *) _digest;
        }
    };
}

#endif

