#ifndef GREGGEN_HPP
#define GREGGEN_HPP

#include <type_traits>

#define MILLION 1'000'000.0l
#define BILLION 1'000'000'000.0l
#define TRILLION 1'000'000'000'000.0l

namespace gml {
    template <typename T>
    concept Numeric = requires (T value) {
        T{1};
        // T{1.0l};
        T{T{}};
        {value + value} -> std::convertible_to<T>;
        {1 + value} -> std::convertible_to<T>;
        {value + 1} -> std::convertible_to<T>;
        {1.0l + value} -> std::convertible_to<T>;
        {value + 1.0l} -> std::convertible_to<T>;
        {1 - value} -> std::convertible_to<T>;
        {value - 1} -> std::convertible_to<T>;
        {1.0l - value} -> std::convertible_to<T>;
        {value - 1.0l} -> std::convertible_to<T>;
        {value*value} -> std::convertible_to<T>;
        {1*value} -> std::convertible_to<T>;
        {value*1} -> std::convertible_to<T>;
        {1.0l*value} -> std::convertible_to<T>;
        {value*1.0l} -> std::convertible_to<T>;
        {1/value} -> std::convertible_to<T>;
        {value/1} -> std::convertible_to<T>;
        {1.0l/value} -> std::convertible_to<T>;
        {value/1.0l} -> std::convertible_to<T>;
        {value += value} -> std::convertible_to<T>;
        {value += 1} -> std::convertible_to<T>;
        {value += 1.0l} -> std::convertible_to<T>;
        {value -= 1} -> std::convertible_to<T>;
        {value -= 1.0l} -> std::convertible_to<T>;
        {value *= value} -> std::convertible_to<T>;
        {value *= 1} -> std::convertible_to<T>;
        {value *= 1.0l} -> std::convertible_to<T>;
        {value /= 1} -> std::convertible_to<T>;
        {value /= 1.0l} -> std::convertible_to<T>;
    };
    template <typename U, typename V>
    using addComT = decltype(std::declval<U>() + std::declval<V>());
    template <typename U, typename V>
    using subComT = decltype(std::declval<U>() - std::declval<V>());
    template <typename U, typename V>
    using mulComT = decltype(std::declval<U>()*std::declval<V>());
    template <typename U, typename V>
    using divComT = decltype(std::declval<U>()/std::declval<V>());
    typedef uint64_t ull_t;
    namespace gen {
        bool memcopy(void *dst, const void *src, size_t elem_size, size_t num_elems) {
            if (!dst || !src || !elem_size || !num_elems)
                return false;
            char *dest = static_cast<char*>(dst);
            const char *source = static_cast<const char*>(src);
            unsigned long long tot_bytes = elem_size*num_elems;
            while (tot_bytes --> 0)
                *dest++ = *source++;
            return true;
        }
        bool memswap(void *obj1, void *obj2, size_t obj_size) {
            if (!obj1 || !obj2 || !obj_size)
                return false;
            char *o1 = static_cast<char*>(obj1);
            char *o2 = static_cast<char*>(obj2);
            char temp;
            while (obj_size --> 0) {
                temp = *o1;
                *o1++ = *o2;
                *o2++ = temp;
            }
            return true;
        }
        template <typename T>
        bool memswap(T *obj1, T *obj2) {
            if (!obj1 || !obj2)
                return false;
            size_t size = sizeof(T);
            if (!(size % sizeof(uint64_t))) {
                uint64_t *o1 = (uint64_t*) (obj1);
                uint64_t *o2 = (uint64_t*) (obj2);
                size /= sizeof(uint64_t);
                uint64_t temp;
                while (size --> 0) {
                    temp = *o1;
                    *o1++ = *o2;
                    *o2++ = temp;
                }
                return true;
            }
            if (!(size % sizeof(uint32_t))) {
                uint32_t *o1 = (uint32_t*) (obj1);
                uint32_t *o2 = (uint32_t*) (obj2);
                size /= sizeof(uint32_t);
                uint32_t temp;
                while (size --> 0) {
                    temp = *o1;
                    *o1++ = *o2;
                    *o2++ = temp;
                }
                return true;
            }
            if (!(size % sizeof(uint16_t))) {
                uint16_t *o1 = (uint16_t*) (obj1);
                uint16_t *o2 = (uint16_t*) (obj2);
                size /= sizeof(uint16_t);
                uint16_t temp;
                while (size --> 0) {
                    temp = *o1;
                    *o1++ = *o2;
                    *o2++ = temp;
                }
                return true;
            }
            char *o1 = (char*) (obj1);
            char *o2 = (char*) (obj2);
            char temp;
            while (size --> 0) {
                temp = *o1;
                *o1++ = *o2;
                *o2++ = temp;
            }
            return true;
        }
        template <typename to, typename from> requires (std::is_convertible<from, to>::value)
        bool copy(to *dst, const from *src, uint64_t _num) {
            if (!dst || !src)
                return false;
            while (_num --> 0)
                *dst++ = *src++;
            return true;
        }
        uint64_t strlen_c(const char *str) {
            if (!str)
                return -1;
            uint64_t count = 0;
            while (*str++) ++count;
            return count;
        }
        char *strcpy_c(char *dst, const char *src) {
            if (!dst || !src)
                return nullptr;
            char *org = dst;
            while (*src)
                *dst++ = *src++;
            *dst = 0;
            return org;
        }
        bool endswith(const char *str, const char *with) {
            if (!str || !with || !*str)
                return false;
            uint64_t slen = 1;
            while (*++str) ++slen;
            uint64_t wlen = 0;
            while (*with++) ++wlen;
            --with;
            if (wlen > slen)
                return false;
            while (wlen --> 0)
                if (*--str != *--with)
                    return false;
            return true;
        }
        bool equals(const char *s1, const char *s2) {
            if (!s1 || !s2)
                return false;
            while (*s1 || *s2)
                if (*s1++ != *s2++)
                    return false;
            return true;
        }
        template <typename InputIt> requires requires (InputIt it, decltype(*std::declval<InputIt>()) val) {
            *it;
            it++;
            {val += val} -> std::same_as<decltype(val)>;
        }
        auto sum(InputIt _begin, InputIt _end) {
            typename std::remove_reference<decltype(*_begin)>::type _tot{};
            while (_begin != _end)
                _tot += *_begin++;
            return _tot;
        }
        template <typename InputIt> requires requires (InputIt it, decltype(*std::declval<InputIt>()) val) {
            *it;
            ++it;
            {val *= val} -> std::same_as<decltype(val)>;
        }
        auto product(InputIt _begin, InputIt _end) {
            if (_begin == _end)
                return typename std::remove_reference<decltype(*std::declval<InputIt>())>::type{};
            auto _prod = *_begin;
            while (++_begin != _end)
                _prod *= *_begin;
            return _prod;
        }
    }
}
#endif
