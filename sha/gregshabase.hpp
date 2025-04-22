#ifndef GREGSHABASE_HPP
#define GREGSHABASE_HPP

#include <bit>
#include <cstdint>
#include <concepts>
#include <cinttypes>

namespace gtd {
    template <std::integral T>
    T rotr(T v, uint64_t n) {
        static constexpr size_t w = sizeof(T)*8;
        n = n % w;
        return (v >> n) | (v << (w - n));
    }
    template <std::integral T>
    T rotl(T v, uint64_t n) {
        static constexpr size_t w = sizeof(T)*8;
        n = n % w;
        return (v << n) | (v >> (w - n));
    }
    template <std::integral T> requires (sizeof(T) == 4 || sizeof(T) == 8)
    T ch(T x, T y, T z) {
        return (x & y) ^ (~x & z);
    }
    template <std::integral T> requires (sizeof(T) == 4 || sizeof(T) == 8)
    T maj(T x, T y, T z) {
        return (x & y) ^ (x & z) ^ (y & z);
    }
    uint32_t Sigma0_256(uint32_t x) {
        return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22);
    }
    uint32_t Sigma1_256(uint32_t x) {
        return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25);
    }
    uint32_t sigma0_256(uint32_t x) {
        return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3);
    }
    uint32_t sigma1_256(uint32_t x) {
        return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10);
    }
    uint64_t Sigma0_512(uint64_t x) {
        return rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39);
    }
    uint64_t Sigma1_512(uint64_t x) {
        return rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41);
    }
    uint64_t sigma0_512(uint64_t x) {
        return rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7);
    }
    uint64_t sigma1_512(uint64_t x) {
        return rotr(x, 19) ^ rotr(x, 61) ^ (x >> 6);
    }
    bool sha256_memcpy(uint32_t *dst, const char *src, uint64_t nwords) {
        if (!dst || !src || !nwords)
            return false;
        while (nwords --> 0) {
            *dst++ = (*src << 24) + (*(src + 1) << 16) + (*(src + 2) << 8) + *(src + 3);
            src += 4;
        }
        return true;
    }
    consteval bool little_endian() {
        static_assert(std::endian::native == std::endian::little ||
                      std::endian::native == std::endian::big,
                      "Unrecognised endianness!\n"
        );
        return std::endian::native == std::endian::little;
    }
    void print_hex(const void *data, uint64_t size) {
        if (!data || !size)
            return;
        const unsigned char *ptr = (const unsigned char *) data;
        while (size --> 0)
            printf("%02x", *ptr++);
        putchar('\n');
    }
    class shasum {
    protected:
        uint64_t _size = 0;
        bool _done = false;
    public:
        shasum() = default;
        shasum(uint64_t size) : _size{size} {}
        virtual bool update(const char *data, uint64_t size) noexcept = 0;
        virtual char *digest() const noexcept = 0;
        virtual char *hexdigest() const noexcept = 0;
        virtual void reset() noexcept = 0;
    };
    template <uint16_t N>
    class sha_const {

    };
}
#endif
