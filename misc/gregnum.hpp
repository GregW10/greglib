#ifndef GREGNUM_H
#define GREGNUM_H

#include <vector>
#include <cinttypes>
#include <cstdint>
#include <iostream>
#include <iomanip>

namespace gtd {
    class big_integer {
        bool neg = false;
        // bool zero = true;
        std::vector<uint64_t> data;
        constexpr static uint64_t half = 0x0fffffffffffffff;
        constexpr static uint64_t big  = 0xffffffffffffffff;
        static void add_to_existing(std::vector<uint64_t> &dat, const uint64_t *to_add, uint64_t size) {
            // std::cout << dat.size() << ", " << size << std::endl;
            if (size > dat.size()) {
                if (*(to_add + size - 1) > half)
                    dat.resize(size + 1);
                else
                    dat.resize(size);
                //std::cout << "NOOOOO" << std::endl;
            }
            // std::cout << std::dec << dat.back() << ", " << half << std::endl;
            if (dat.back() > half) {
                dat.resize(size + 1);
                // std::cout << "NOOOOO" << std::endl;
            }
            uint64_t *dptr = dat.data();
            uint64_t res;
            bool carry = 0;
            while (size --> 0) {
                res = *dptr + *to_add + carry;
                carry = res < *dptr || res < *to_add++;
                *dptr++ = res;
            }
            if (carry)
                *dptr = 1;
        }
        static void sub_from_existing(std::vector<uint64_t> &dat, const uint64_t *to_sub, uint64_t size);
    public:
        big_integer() : data{std::initializer_list<uint64_t>{0}} {}
        big_integer(uint64_t num, bool negative = false) : neg{negative}, data{std::initializer_list<uint64_t>{num}} {}
        big_integer &operator+=(uint64_t num) {
            if (this->neg) {
                this->sub_from_existing(this->data, &num, 1);
                return *this;
            }
            add_to_existing(this->data, &num, 1);
            return *this;
        }
        big_integer &operator+=(const big_integer &other) {
            if (&other == this) {
                // multiply by 2
            }
            add_to_existing(this->data, other.data.data(), other.data.size());
            return *this;
        }
        operator bool() const noexcept {
            for (const uint64_t &v : this->data)
                if (v)
                    return true;
        }
        friend std::ostream &operator<<(std::ostream &os, const big_integer &bi) {
            // if (!bi)
            //     return os << '0';
            if (bi.neg)
                os << '-';
            os << std::hex;
            uint64_t counter = bi.data.size();
            const uint64_t *dptr = bi.data.data() + counter - 1;
            while (counter > 0) {
                if (*dptr)
                    break;
                --dptr;
                --counter;
            }
            if (counter == big)
                return os << '0';
            while (counter --> 0)
                os << *dptr--;
            return os;
        }
    };
}

#endif
