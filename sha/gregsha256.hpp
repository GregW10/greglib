#ifndef GREGSHA256_HPP
#define GREGSHA256_HPP

#include <cstring>
#include <iostream>
#include "gregshabase.hpp"
#include "gregsha-constants.hpp"

namespace gtd {
    class sha256 { // : public shasum {
        uint64_t _size = 0; // total number of bytes of data seen so far
        bool     _done = false;
        uint32_t _H[8] = { // initial hash values
            sha256_init_hash_vals[0],
            sha256_init_hash_vals[1],
            sha256_init_hash_vals[2],
            sha256_init_hash_vals[3],
            sha256_init_hash_vals[4],
            sha256_init_hash_vals[5],
            sha256_init_hash_vals[6],
            sha256_init_hash_vals[7]
        };
        uint32_t _M[16]; // message block
        uint32_t _W[64]; // message schedule
        unsigned char _par[64]{}; // partially-filled block
        unsigned char *_nbyte = (unsigned char *) _par; // pointer to the next char to assign in the partially-filled block
        unsigned char _npos = 0; // position of the next char to assign in the partially-filled block
        unsigned char _digest[32]; // final digest as raw char array -> only set if "_done" is true
        void _copy_block(const unsigned char *dat) {
            uint32_t *mptr = _M;
            uint32_t *wptr = _W;
            char counter = 16;
            while (counter --> 0) {
                *mptr = (*dat << 24) + (*(dat + 1) << 16) + (*(dat + 2) << 8) + *(dat + 3);
                *wptr++ = *mptr++; // remember to define a separate function for this if BIG-ENDIAN
                dat += 4;
            }
        }
        void _process_block(const unsigned char *data) {
            // print_hex(data, 64);
            if constexpr (little_endian())
                _copy_block(data);
            else {
                memcpy(_M, data, 64);
                memcpy(_W, data, 64);
            }
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
            const uint32_t *kptr;
            // if constexpr (little_endian())
                // kptr = sha256_constants_le;
            // else
            kptr = sha256_constants;//_be;
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
            }
            _H[0] += a;
            _H[1] += b;
            _H[2] += c;
            _H[3] += d;
            _H[4] += e;
            _H[5] += f;
            _H[6] += g;
            _H[7] += h;
        }
        void _process_data(const unsigned char *data, uint64_t size) {
            _size += size;
            uint64_t bytes = _npos + size;
            if (bytes < 64) {
                //std::cout << "small\n";
                memcpy(_nbyte, data, size);
                _nbyte += size;
                _npos += size;
                return;
            }
                //std::cout << "big\n";
            uint64_t N;
            if (bytes % 64) {
                //std::cout << "here\n";
                N = 64 - _npos; // N just acting as a temp. variable
                memcpy(_nbyte, data, N);
                data += N;
                _process_block(_par);
                _npos += (size % 64);
                _npos = _npos % 64; // WATCH OUT
                _nbyte = ((unsigned char *) _par) + _npos;
                N = (bytes - N)/64;
                while (N --> 0) {
                    _process_block(data);
                    data += 64;
                }
                memcpy(_par, data, _npos);
            } else {
                if (_npos) {
                    N = 64 - _npos;
                    memcpy(_nbyte, data, N);
                    data += N;
                    _npos = 0;
                    _nbyte = (unsigned char *) _par;
                    N = bytes/64 - 1;
                    _process_block(_par);
                } else
                    N = bytes/64;
                while (N --> 0) {
                    _process_block(data);
                    data += 64;
                }
            }
        }
        void _finish() { // leaves object in an unusable state
            if (_npos <= 55) {
                *_nbyte++ = 0b10000000;
                ++_npos;
                bye:
                while (_npos < 56) {
                    *_nbyte++ = 0;
                    ++_npos;
                }
                _size *= 8;
                _par[56] = (_size >> 56);
                _par[57] = (_size >> 48);
                _par[58] = (_size >> 40);
                _par[59] = (_size >> 32);
                _par[60] = (_size >> 24);
                _par[61] = (_size >> 16);
                _par[62] = (_size >> 8);
                _par[63] =  _size;
                // print_hex(_par, 64);
                this->_process_block(_par);
                this->_set_final_digest();
                _done = true;
                return;
            }
            *_nbyte++ = 0b10000000;
            ++_npos;
            while (_npos < 64) {
                *_nbyte++ = 0;
                ++_npos;
            }
            // print_hex(_par, 64);
            _process_block(_par);
            _nbyte = _par;
            _npos = 0;
            goto bye;
        }
        void _init_hash_vals() {
            _H[0] = sha256_init_hash_vals[0];
            _H[1] = sha256_init_hash_vals[1];
            _H[2] = sha256_init_hash_vals[2];
            _H[3] = sha256_init_hash_vals[3];
            _H[4] = sha256_init_hash_vals[4];
            _H[5] = sha256_init_hash_vals[5];
            _H[6] = sha256_init_hash_vals[6];
            _H[7] = sha256_init_hash_vals[7];
        }
        void _set_final_digest() {
            unsigned char *dptr = _digest;
            uint32_t *hptr = (uint32_t *) _H;
            char counter = 8;
            while (counter --> 0) {
                *dptr++ = (*hptr >> 24);
                *dptr++ = (*hptr >> 16);
                *dptr++ = (*hptr >>  8);
                *dptr++ =  *hptr++;
            }
        }
    public:
        sha256() = default;
        sha256(const void *data, uint64_t size, bool _hash = false) {
            if (!data || !size)
                return;
            _process_data((unsigned char *) data, size);
            if (_hash)
                this->_finish();
        }
        bool add_data(const void *data, uint64_t size) noexcept {
            if (!data || !size || _done)
                return false;
            _process_data((unsigned char *) data, size);
            return true;
        }
        void hash() noexcept {
            if (!this->_done)
                this->_finish();
        }
        const unsigned char *digest() noexcept {
            if (!_done)
                this->_finish();
            return _digest;
        }
        void reset() noexcept {
            _done = false;
            _npos = 0;
            _nbyte = (unsigned char *) _par;
            _size = 0;
            this->_init_hash_vals();
        }
        sha256 &operator=(const sha256 &other) noexcept {
            if (&other == this)
                return *this;
            this->_size = other._size;
            this->_npos = other._npos;
            this->_nbyte = other._nbyte;
            this->_done = other._done;
            memcpy(this->_H, other._H, 32);
            if (this->_npos)
                memcpy(this->_par, other._par, this->_npos);
            if (this->_done)
                memcpy(this->_digest, other._digest, 32);
            return *this;
        }
        friend std::ostream &operator<<(std::ostream &os, const sha256 &sha) {
            static const char *charset = "0123456789abcdef";
            if (!sha._done)
                return os << "[unfinished]";
            char buff[65];
            buff[64] = 0;
            const unsigned char *hptr = sha._digest;
            char *bptr = buff;
            unsigned char counter = 32;
            while (counter --> 0) {
                *bptr++ = *(charset + (*hptr >> 4));
                *bptr++ = *(charset + (*hptr++ & 0x0f));
            }
            return os << buff;
        }
        friend bool operator<(const sha256 &one, const sha256 &two) {
            const uint32_t *optr = (uint32_t *) one._H;
            const uint32_t *tptr = (uint32_t *) two._H;
            char counter = 8;
            while (counter --> 0) {
                if (*optr < *tptr)
                    return true;
                if (*optr++ > *tptr++)
                    return false;
            }
            return false;
        }
        friend bool operator<=(const sha256 &one, const sha256 &two) {
            const uint32_t *optr = (uint32_t *) one._H;
            const uint32_t *tptr = (uint32_t *) two._H;
            char counter = 8;
            while (counter --> 0) {
                if (*optr < *tptr)
                    return true;
                if (*optr++ > *tptr++)
                    return false;
            }
            return true;
        }
        friend bool operator==(const sha256 &one, const sha256 &two) {
            const uint32_t *optr = (uint32_t *) one._H;
            const uint32_t *tptr = (uint32_t *) two._H;
            char counter = 8;
            while (counter --> 0)
                if (*optr++ != *tptr++)
                    return false;
            return true;
        }
        friend bool operator!=(const sha256 &one, const sha256 &two) {
            const uint32_t *optr = (uint32_t *) one._H;
            const uint32_t *tptr = (uint32_t *) two._H;
            char counter = 8;
            while (counter --> 0)
                if (*optr++ != *tptr++)
                    return true;
            return false;
        }
        friend bool operator>(const sha256 &one, const sha256 &two) {
            return two < one;
        }
        friend bool operator>=(const sha256 &one, const sha256 &two) {
            return two <= one;
        }
    };
}
#endif

