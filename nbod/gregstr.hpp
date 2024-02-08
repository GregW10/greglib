// Copyright (c) 2023 Gregor Anton Randall Hartl Watters
// This software is protected under the MIT license. Please see the LICENSE file for more information.

#ifndef GREGSTR_HPP
#define GREGSTR_HPP

//
// Created by Gregor Hartl Watters on 04/04/2022.
//

#ifndef __cplusplus
#error "The gregstr.hpp header file is a C++ header file only."
#endif

#include <cstdlib>
#include <string>
#include <iostream>
#include <iterator>
#include <cstddef>
#include <limits>
#include <vector>
#include <sys/stat.h>
#include <functional>
#include <ctime>

#ifndef _WIN32
#include <unistd.h>
#include <pwd.h>
#define FILE_SEP '/'
#define RESET_TXT_FLAGS "\033[0m"
#define BLACK_TXT_START "\033[30m"
#define RED_TXT_START "\033[31m"
#define GREEN_TXT_START "\033[32m"
#define YELLOW_TXT_START "\033[33m"
#define BLUE_TXT_START "\033[34m"
#define MAGENTA_TXT_START "\033[35m"
#define CYAN_TXT_START "\033[36m"
#define WHITE_TXT_START "\033[37m"
#define BOLD_TXT_START "\033[1m"
#define UNDERLINED_TXT_START "\033[4m"
#define BLACK_TXT(str) BLACK_TXT_START str RESET_TXT_FLAGS
#define RED_TXT(str) RED_TXT_START str RESET_TXT_FLAGS
#define GREEN_TXT(str) GREEN_TXT_START str RESET_TXT_FLAGS
#define YELLOW_TXT(str) YELLOW_TXT_START str RESET_TXT_FLAGS
#define BLUE_TXT(str) BLUE_TXT_START str RESET_TXT_FLAGS
#define MAGENTA_TXT(str) MAGENTA_TXT_START str RESET_TXT_FLAGS
#define CYAN_TXT(str) CYAN_TXT_START str RESET_TXT_FLAGS
#define WHITE_TXT(str) WHITE_TXT_START str RESET_TXT_FLAGS
#define BOLD_TXT(str) BOLD_TXT_START str RESET_TXT_FLAGS
#define UNDERLINED_TXT(str) UNDERLINED_TXT_START str RESET_TXT_FLAGS
#else
#include <shlobj.h>
#define FILE_SEP '\\'
#define RESET_TXT_FLAGS ""
#define BLACK_TXT_START ""
#define RED_TXT_START ""
#define GREEN_TXT_START ""
#define YELLOW_TXT_START ""
#define BLUE_TXT_START ""
#define MAGENTA_TXT_START ""
#define CYAN_TXT_START ""
#define WHITE_TXT_START ""
#define BOLD_TXT_START ""
#define UNDERLINED_TXT_START ""
#define BLACK_TXT(str) str
#define RED_TXT(str) str
#define GREEN_TXT(str) str
#define YELLOW_TXT(str) str
#define BLUE_TXT(str) str
#define MAGENTA_TXT(str) str
#define CYAN_TXT(str) str
#define WHITE_TXT(str) str
#define BOLD_TXT(str) str
#define UNDERLINED_TXT(str) str
#endif

namespace gtd {
    unsigned char to_upper(unsigned char ch);
    unsigned char to_lower(unsigned char ch);
    size_t strlen_c(const char *str);
    char *memset_c(char *str, char ch, size_t n_chars);
    char *strcpy_c(char *dest, const char *source);
    inline char *setchr_c(char *str, char ch, size_t pos);
    int strcmp_c(const char *str1, const char *str2);
    int strncmp_c(const char *str1, const char *str2, size_t n);
    bool isdigit_c(char ch);
    class NullPointerError : public std::exception {
    private:
        char *message;
    public:
        NullPointerError() {
            const char def_msg[] = "A null pointer cannot be passed as a string.";
            size_t length = strlen_c(def_msg);
            message = new char[length + 1];
            memset_c(message, '\0', length + 1);
            strcpy_c(message, def_msg);
        }
        explicit NullPointerError(const char *msg) {
            if (msg == nullptr)
                throw std::invalid_argument("nullptr passed as message.");
            size_t len_p1 = strlen_c(msg) + 1;
            message = new char[len_p1];
            memset_c(message, '\0', len_p1);
            strcpy_c(message, msg);
        }
        [[nodiscard]] const char *what() const noexcept override {
            return message;
        }
        ~NullPointerError() override {
            delete[] message;
        }
    };
    class EmptyStringError : public std::exception {
    private:
        static constexpr const char *default_msg = "This operation cannot be performed on an emtpy string.";
        char *msg = nullptr;
    public:
        EmptyStringError() = default;
        explicit EmptyStringError(const char *message) {
            if (message == nullptr)
                throw NullPointerError();
            size_t length = strlen_c(message);
            msg = new char[length + 1];
            strcpy_c(msg, message);
            setchr_c(msg, 0, length);
        }
        [[nodiscard]] const char *what() const noexcept override {
            if (msg == nullptr)
                return default_msg;
            return msg;
        }
        ~EmptyStringError() override {
            delete[] msg; // safe to delete nullptr in C++
        }
    };
    template <typename T>
    concept charIt = requires (T it) {
        {*it} -> std::convertible_to<char>;
        {++it} -> std::same_as<T&>;
        {it++} -> std::same_as<T>;
        {--it} -> std::same_as<T&>;
        {it--} -> std::same_as<T>;
    };
    class String {
    private:
        char *data = nullptr;
        char *string = nullptr;
        size_t length_w_null = 0;
        size_t size = 0;
        size_t space_front = 0;
        size_t space_back = 0;
        bool is_empty = true;
        bool shrunk = false; // only true after calling shrink_to_fit() and before memory is resized again
        bool start_left = false; // only true if push_left is set to true when a String object is constructed
        static const size_t MIN_SIZE = 32;
        static inline bool pow_of_2(size_t num) {
            return (num != 0) && ((num & (num - 1)) == 0); // any power of two will be a 1 followed by zeros in binary,
        } // so, e.g., 10000 & 01111 = 00000 (which is zero) - bitwise AND of n and n-1 for n power of 2 is always zero
        void set_size(bool def_double = true) {
            if (shrunk && !pow_of_2(size)) {
                if (size < MIN_SIZE) {
                    size = MIN_SIZE;
                }
                else {
                    unsigned char n = 0;
                    while ((size >> (++n)) != 0); // find position from RHS of the first non-zero bit (n)
                    size = (~size >> (sizeof(size_t)*8 - n)) + 1; // bitwise negation makes all leading zero bits 1,
                } // then this is bit-shifted by n, which yields 2^n - 1, so +1 after
            }
            if (def_double) {
                size *= 2;
            }
            while (size < length_w_null) {
                size *= 2;
            }
        }
        void constructor(const char *str, bool after_empty = false) {
            shrunk = false;
            size = MIN_SIZE;
            length_w_null = strlen_c(str) + 1;
            set_size(false);
            if (after_empty) { // this is to avoid re-assigning the memory if not necessary
                if (size != MIN_SIZE) {
                    delete[] data;
                    data = new char[size];
                }
            }
            else {
                data = new char[size];
            }
            memset_c(data, '\0', size);
            string = data + get_first_pos();
            // setchr_c(string, '\0', length_w_null - 1);
            strcpy_c(string, str);
            space_front = get_first_pos();
            space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front : space_front + 1);
            is_empty = false;
        }
        void empty_constructor() {
            shrunk = false;
            size = MIN_SIZE;
            data = new char[size];
            memset_c(data, '\0', size);
            string = nullptr;
            length_w_null = 0;
            if (start_left) {
                space_front = 0;
                space_back = size - length_w_null;
            }
            else {
                space_front = size / 2;
                space_back = size / 2;
            }
            is_empty = true;
        }
        template <charIt It>
        void it_constructor(It beg, It ending) {
            shrunk = false;
            length_w_null = ending - beg + 1;
            size = MIN_SIZE;
            set_size(false);
            data = new char[size];
            memset_c(data, '\0', size);
            string = data + get_first_pos();
            char *ptr = string;
            while (beg != ending)
                *ptr++ = *beg++;
            setchr_c(string, '\0', length_w_null - 1);
            space_front = get_first_pos();
            space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front : space_front + 1);
            is_empty = false;
        }
        inline static bool is_even(size_t num) {
            return num % 2 == 0;
        }
        inline static bool is_punc(const char &ch) { // reference because it is a private func. only used with l-values
            return (ch > 0 && ch < 48) || (ch > 57 && ch < 65) || (ch > 90 && ch < 97) || ch > 122;
        }
        [[nodiscard]] inline size_t get_first_pos() const noexcept {
            return start_left ? 0 : (is_even(length_w_null) ? size/2 - length_w_null/2 : size/2 - (length_w_null+1)/2);
        }
    public:
        static const size_t nopos = -1;
        class RevIterator;
        class Iterator {
        public:
            Iterator() : ptr(nullptr) {}
            Iterator(char *pointer) : ptr(pointer) {}
            char &operator*() const {
                return *ptr;
            }
            char *operator->() const {
                return ptr;
            }
            virtual Iterator &operator++() {
                ++ptr;
                return *this;
            }
            Iterator operator++(int) {
                Iterator tmp = *this;
                ++ptr;
                return tmp;
            }
            virtual Iterator &operator--() {
                --ptr;
                return *this;
            }
            Iterator operator--(int) {
                Iterator tmp = *this;
                --ptr;
                return tmp;
            }
            friend bool operator==(const Iterator &A, const Iterator &B);
            friend bool operator!=(const Iterator &A, const Iterator &B);
            friend bool operator<(const Iterator &A, const Iterator &B);
            friend bool operator>(const Iterator &A, const Iterator &B);
            friend bool operator<=(const Iterator &A, const Iterator &B);
            friend bool operator>=(const Iterator &A, const Iterator &B);
            Iterator &operator=(const Iterator &other) {
                if (*this == other) {
                    return *this;
                }
                this->ptr = other.ptr;
                return *this;
            }
            Iterator &operator=(char *pointer) {
                if (this->ptr == pointer) {
                    return *this;
                }
                this->ptr = pointer;
                return *this;
            }
            Iterator operator+(std::ptrdiff_t offset) {
                return {ptr + offset};
            }

            Iterator operator-(std::ptrdiff_t offset) {
                return {ptr - offset};
            }
            friend std::ostream &operator<<(std::ostream &out, const Iterator &A);
            friend std::ptrdiff_t operator-(const Iterator &A, const Iterator &B);
            friend std::ptrdiff_t operator-(const String::RevIterator &A, const String::RevIterator &B);
            ~Iterator() = default;
        protected:
            char *ptr;
        };
        class RevIterator : virtual public Iterator {
        public:
            RevIterator() = default;
            RevIterator(char *pointer) : Iterator(pointer) {}
            RevIterator &operator++() override {
                --ptr;
                return *this;
            }
            RevIterator operator++(int) { // not pretty, but necessary
                RevIterator tmp = *this;
                --ptr;
                return tmp;
            }
            RevIterator &operator--() override {
                ++ptr;
                return *this;
            }
            RevIterator operator--(int) { // same
                RevIterator tmp = *this;
                ++ptr;
                return tmp;
            }
            friend bool operator<(const RevIterator &A, const RevIterator &B);
            friend bool operator>(const RevIterator &A, const RevIterator &B);
            friend bool operator<=(const RevIterator &A, const RevIterator &B);
            friend bool operator>=(const RevIterator &A, const RevIterator &B);
            RevIterator operator+(std::ptrdiff_t offset) {
                return {ptr - offset};
            }
            RevIterator operator-(std::ptrdiff_t offset) {
                return {ptr + offset};
            }
        };
        class ConstIterator : virtual public Iterator {
        private:
            const char *cptr = nullptr;
        public:
            ConstIterator() = default;
            ConstIterator(char *pointer) : Iterator(pointer), cptr(pointer){}
            virtual const char *operator->() {
                cptr = ptr;
                return cptr;
            }
            virtual const char &operator*() {
                cptr = ptr;
                return *cptr;
            }
        };
        class ConstRevIterator : public RevIterator, public ConstIterator {
        public:
            ConstRevIterator() = default;
            ConstRevIterator(char *pointer) : Iterator(pointer), RevIterator(pointer), ConstIterator(pointer) {}
            const char *operator->() override {
                return ConstIterator::operator->();
            }
            const char &operator*() override {
                return ConstIterator::operator*();
            }
        };
        String() {
            start_left = false;
            empty_constructor();
        }
        explicit String(bool push_left) {
            start_left = push_left;
            empty_constructor();
        }
        String(char ch, bool push_left = false) {
            start_left = push_left;
            if (!ch) {
                empty_constructor();
                return;
            }
            const char str[2]{ch, '\0'};
            constructor(str);
        }
        String(const char *str, bool push_left = false) {
            start_left = push_left;
            if (str == nullptr) {
                throw NullPointerError();
            }
            if (!strlen_c(str)) {
                empty_constructor();
                return;
            }
            constructor(str);
        }
        String(const std::string &str, bool push_left = false) {
            start_left = push_left;
            if (str.empty()) {
                empty_constructor();
                return;
            }
            constructor(str.c_str());
        }
        String(const String &str, bool push_left = false) {
            start_left = push_left;
            if (str.empty()) {
                empty_constructor();
                return;
            }
            constructor(str.c_str());
        }
        String(String &&str) noexcept { // move constructor
            this->operator=(std::move(str));
        }
        template <charIt ForwardIterator>
        String(ForwardIterator beg, ForwardIterator ending, bool push_left = false) {
            start_left = push_left;
            if (beg >= ending) {
                empty_constructor();
                return;
            }
            it_constructor(beg, ending);
        }
        String &append_back(const char *str) {
            if (str == nullptr) {
                throw NullPointerError();
            }
            if (is_empty) {
                constructor(str, true);
                return *this;
            }
            size_t l = strlen_c(str);
            size_t old_l = length_w_null - 1;
            length_w_null += l;
            if (l > space_back) {
                if (shrunk) {
                    set_size(false);
                    shrunk = false;
                }
                else {
                    set_size();
                }
                char *old_data = data;
                char *old_str = string;
                data = new char[size];
                memset_c(data, '\0', size);
                string = data + get_first_pos();
                strcpy_c(string, old_str);
                delete[] old_data;
                space_front = get_first_pos();
                space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front : space_front+1);
                strcpy_c(string + old_l, str);
                return *this;
            }
            strcpy_c(string + old_l, str);
            space_back -= l;
            return *this;
        }
        String &append_back(const std::string &str) {
            if (str.empty()) {
                return *this;
            }
            return append_back(str.c_str());
        }
        String &append_back(const String &str) {
            if (str.is_empty) {
                return *this;
            }
            return append_back(str.c_str());
        }
    private:
        template <typename T> requires (std::is_integral<T>::value) // only here for completeness
        unsigned char process_integral_num(T value) {
            if (value < 0) {
                this->push_back('-');
                value = -value;
            }
            size_t val = value;
            unsigned char n = 1;
            size_t power = 1;
            while (val > 0) {
                val = value;
                power *= 10;
                val /= power;
                ++n;
            }
            while (true) {
                power /= power == 1 ? 1 : 10;
                this->push_back((char) ((value / power) + 48)); // 48 in ASCII = '0'
                value = value % power;
                if (power == 1) {
                    break;
                }
            }
            return n;
        }
    public:
        template <typename T> requires (std::is_fundamental<T>::value)
        String &append_back(T value, size_t num_dec_places = 15) {
            if constexpr (std::is_integral<T>::value) {
                process_integral_num(value);
                return *this;
            }
            unsigned char num;
            size_t org_dec_p = num_dec_places;
            if (value == 0) {
                this->push_back('0');
                if (num_dec_places > 0) {
                    this->push_back('.');
                    while (num_dec_places --> 0)
                        this->push_back('0');
                }
            }
            bool negative = false;
            if (value < 0) {
                value = -value;
                this->push_back('-');
                negative = true;
            }
            T org_org_partial_val = value - ((size_t) value);
            size_t org_org_int_val = (size_t) value;
            T org_val = value;
            size_t n = 1;
            long double floatPart;
            bool reduced = false;
            if (value == 1) {
                this->push_back('1');
                if (num_dec_places > 0) {
                    this->push_back('.');
                    while (num_dec_places --> 0)
                        this->push_back('0');
                }
                return *this;
            }
            if (num_dec_places == 0) {
                floatPart = value - ((size_t) value);
                if (value > 1) {
                    if (floatPart > 0.5)
                        process_integral_num(((size_t) value) + 1);
                    else
                        process_integral_num((size_t) value);
                    return *this;
                }
                if (floatPart > 0.5)
                    this->push_back('1');
                else
                    this->push_back('0');
                return *this;
            }
            if (value > 1) {
                num = process_integral_num((size_t) (value));
                value -= (size_t) value;
                org_val = value;
                reduced = true;
            }
            floatPart = value;
            if (reduced)
                this->push_back('.');
            else
                this->append_back("0.");
            unsigned char count = 0;
            while (value < 1 && count < num_dec_places) {
                n *= 10;
                value = org_val*n; // done to minimise accumulation of floating-point error
                this->push_back('0');
                ++count;
            }
            floatPart = value;
            this->pop_back();
            num_dec_places -= count;
            n = 1;
            size_t max_n = 1000000000000000000;
            max_n *= 10; // to avoid compiler warnings on too-long-to-be-signed integer literals
            while (num_dec_places --> 0 && n < (size_t) max_n)
                n *= 10;
            floatPart *= n;
            size_t intFloatPart = (size_t) floatPart;
            unsigned char remainder = ((size_t) (floatPart*10)) % 10;
            if (10 - remainder > 5) {
                process_integral_num(intFloatPart);
                return *this;
            }
            if (intFloatPart == 0) {
                this->push_back('1');
                return *this;
            }
            if ((intFloatPart + 1) / (n*10) == 1) {
                if (reduced && org_org_partial_val > 0.1) {
                    this->pop_back();
                    while (num --> 0)
                        this->pop_back();
                    if (negative)
                        this->push_back('-');
                    process_integral_num(org_org_int_val + 1);
                    this->push_back('.');
                    while (org_dec_p --> 0)
                        this->push_back('0');
                    return *this;
                }
                if (*(string + length_w_null - 2) != '.')
                    this->pop_back();
                else {
                    this->pop_back();
                    this->pop_back();
                    while (num --> 0) {
                        this->pop_back();
                    }
                    process_integral_num(org_org_int_val + 1);
                    this->push_back('.');
                    while (org_dec_p --> 0) {
                        this->push_back('0');
                    }
                    return *this;
                }
            }
            process_integral_num(intFloatPart + 1);
            return *this;
        }
        String &append_front(const char *str) {
            if (str == nullptr)
                throw NullPointerError();
            if (is_empty) {
                constructor(str, true);
                return *this;
            }
            size_t l = strlen_c(str);
            length_w_null += l;
            char gone = *string;
            if (l > space_front) {
                if (shrunk) {
                    set_size(false);
                    shrunk = false;
                }
                else {
                    set_size();
                }
                char *old_data = data;
                char *old_str = string;
                data = new char[size];
                memset_c(data, '\0', size);
                string = data + get_first_pos();
                strcpy_c(string + l, old_str);
                delete[] old_data;
                strcpy_c(string, str);
                space_front = get_first_pos();
                space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front : space_front+1);
            } else {
                strcpy_c(string - l, str);
                string -= l;
                space_front -= l;
            }
            setchr_c(string, gone, l);
            return *this;
        }
        size_t print_all_chars() const noexcept {
            const char *ptr = data;
            const char *end = data + size;
            for (; ptr < end; ++ptr) {
                if (*ptr == 0) {
                    putchar('0');
                    continue;
                }
                putchar(*ptr);
            }
            putchar(10);
            return size;
        }
        void both() const noexcept { // debugging method, will eventually get rid of it
            space();
            print_all_chars();
        }
        String &append_front(const std::string &str) {
            if (str.empty())
                return *this;
            return append_front(str.c_str());
        }
        String &append_front(const String &str) {
            if (str.is_empty)
                return *this;
            return append_front(str.c_str());
        }
        String &push_back(char ch) {
            if (!ch)
                return *this;
            const char str[2] = {ch, 0};
            return append_back(str);
        }
        String &push_front(char ch) {
            if (!ch)
                return *this;
            const char str[2] = {ch, 0};
            return append_front(str);
        }
        String &pop_back() noexcept {
            if (is_empty)
                return *this;
            if (length_w_null == 2 || length_w_null == 1) {
                this->clear();
                return *this;
            }
            *(string + length_w_null - 2) = 0;
            --length_w_null;
            ++space_back;
            return *this;
        }
        String &pop_front() noexcept {
            if (is_empty)
                return *this;
            if (length_w_null == 2 || length_w_null == 1) {
                this->clear();
                return *this;
            }
            *string++ = 0;
            --length_w_null;
            ++space_front;
            return *this;
        }
        String &reverse() noexcept {
            if (is_empty || length_w_null == 1 || length_w_null == 2)
                return *this;
            char *start = string;
            char *final = start + length_w_null - 2;
            char temp;
            do {
                temp = *start;
                *start = *final;
                *final = temp;
            } while (--final > ++start);
            return *this;
        }
        String &fill(char ch) noexcept {
            if (is_empty || !ch)
                return *this;
            char *ptr = string + length_w_null - 2;
            while (ptr >= string)
                *ptr-- = ch;
            return *this;
        }
        String &fill(size_t starting_index, size_t end_index, char ch) noexcept { // including end_index
            if (starting_index >= end_index || starting_index < 0 || end_index < 0 || !ch || is_empty ||
            starting_index >= length_w_null - 1 || end_index > length_w_null - 1)
                return *this;
            char *start = string + starting_index;
            const char *final = end_index > length_w_null - 1 ? string + length_w_null - 1 : string + end_index;
            while (start <= final)
                *start++ = ch;
            return *this;
        }
        String &fill_data(char ch) noexcept {
            if (!ch)
                return *this;
            char *start = data;
            const char *final = data + size;
            for (; start != final; ++start)
                *start = ch;
            *(--start) = 0;
            string = data;
            length_w_null = size;
            space_front = 0;
            space_back = 0;
            is_empty = false;
            return *this;
        }
        template <charIt ForwardIterator>
        String &assign(ForwardIterator start, ForwardIterator final) {
            if (start >= final)
                return *this;
            delete[] data;
            it_constructor(start, final);
            return *this;
        }
        String substr(size_t starting_index, size_t end_index) const {
            if (starting_index < 0 || end_index < 0 || starting_index >= end_index || starting_index > length_w_null -
                2 || is_empty)
                return {};
            if (end_index > length_w_null - 1)
                end_index = length_w_null - 1;
            // Iterator start = begin() + (long int) starting_index;
            // Iterator final = begin() + (long int) end_index;
            // return {start, final};
            return {string + starting_index, string + end_index};
        }
        String substr(char ch) const {
            if (!ch || is_empty)
                return {};
            char *start = string;
            const char *final = string + length_w_null - 1;
            for (; start != final; ++start)
                if (*start == ch)
                    return {start};
            return {};
        }
        String substr(const char *str) const {
            size_t len = strlen_c(str);
            if (str == nullptr || !len || len > length_w_null - 1 || is_empty)
                return {};
            if (len == length_w_null - 1) {
                if (!strcmp_c(string, str))
                    return {string};
                return {};
            }
            const char *start_m = string;
            const char *start;
            const char *final_m = string + len;
            const char *final_fixed = string + length_w_null - 1;
            size_t count;
            const char *ptr;
            while (final_fixed >= final_m) {
                ptr = str;
                count = 0;
                for (start = start_m; start != final_m; ++count)
                    if (*start++ != *ptr++)
                        break;
                if (count == len)
                    return {start_m};
                ++start_m;
                ++final_m;
            }
            return {};
        }
        String rsubstr(char ch) const {
            if (!ch || is_empty)
                return {};
            const char *start = string;
            const char *final = start + length_w_null - 1;
            --final;
            for (; start != final; --final)
                if (*final == ch)
                    return {string + (final - start)};
            if (*start == ch)
                return {string};
            return {};
        }
        String rsubstr(const char *str) const noexcept {
            size_t len = strlen_c(str);
            if (str == nullptr || !len || len > length_w_null - 1 || is_empty)
                return {};
            if (len == length_w_null - 1) {
                if (!strcmp_c(string, str))
                    return {string};
                return {};
            }
            const char *start_fixed = string;
            const char *start_m = string + length_w_null - 1 - len;
            const char *start;
            const char *final_m = string + length_w_null - 1;
            size_t count = 0;
            while (start_fixed <= start_m) {
                for (start = start_m; start != final_m; ++start) {
                    if (*start != *(str + count))
                        break;
                    ++count;
                }
                if (count == len)
                    return {start_m};
                count = 0;
                --start_m;
                --final_m;
            }
            return {};
        }
        String &insert(char ch, size_t pos = nopos) {// moves the chars in whichever direction has the least chars to
            if (is_empty) {                         // move, or in the direction which does not cause a resize in memory
                const char str[2]{ch, 0};
                constructor(str, true);
                return *this;
            }
            if (pos > length_w_null - 2) {
                if (space_back != 0 || (space_front == 0 && space_back == 0)) {
                    push_back(ch);
                    return *this;
                }
                pos = length_w_null - 1;
            }
            if (!pos) {
                push_front(ch);
                return *this;
            }
            char *to_place = string + pos;
            if (space_front != 0 || space_back != 0) {
                if (pos >= length_w_null / 2 || (pos < length_w_null / 2 && space_front == 0)) {
                    char *ptr = string + length_w_null - 2;
                    char *end = string + length_w_null - 1;
                    while (ptr >= to_place)
                        *end-- = *ptr--;
                    --space_back;
                    *++to_place = ch;
                } else {
                    char *beg = string - 1;
                    char *ptr = string;
                    while (ptr <= to_place)
                        *beg++ = *ptr++;
                    --string;
                    --space_front;
                    *--to_place = ch;
                }
                ++length_w_null;
                return *this;
            }
            return *this;
        }
        String &to_upper() noexcept {
            if (is_empty)
                return *this;
            char *ptr = string + length_w_null - 1;
            while (--ptr >= string)
                if (*ptr >= 97 && *ptr <= 122)
                    *ptr -= 32;
            return *this;
        }
        String &to_lower() noexcept {
            if (is_empty)
                return *this;
            char *ptr = string + length_w_null - 1;
            while (--ptr >= string)
                if (*ptr >= 65 && *ptr <= 90)
                    *ptr += 32;
            return *this;
        }
        bool isnumeric() const noexcept {
            if (is_empty)
                return false;
            const char *ptr = string;
            const char *final = ptr + length_w_null - 1;
            bool has_dot = false;
            while (ptr != final) {
                if (*ptr == '.') {
                    if (has_dot)
                        return false;
                    has_dot = true;
                    ++ptr;
                    continue;
                }
                if (!isdigit_c(*ptr++))
                    return false;
            }
            return true;
        }
        bool contains(char ch) const noexcept {
            if (is_empty)
                return false;
            char *start = string;
            const char *final = start + length_w_null - 1;
            while (start != final)
                if (*start++ == ch)
                    return true;
            return false;
        }
        bool contains(const char *str) const noexcept {
            size_t length = strlen_c(str);
            if (is_empty || str == nullptr || !*str || length > length_w_null - 1)
                return false;
            if (length == length_w_null - 1)
                return *this == str;
            const char *sptr;
            const char *ptr = string;
            const char *end_ptr = string + length_w_null - length;
            const char *vptr;
            size_t i;
            size_t count;
            for (; ptr <= end_ptr; ++ptr) {
                for (count = 0, i = 0, sptr = str, vptr = ptr; i < length; ++i, ++sptr, ++vptr) {
                    if (*sptr != *vptr)
                        break;
                    ++count;
                }
                if (count == length)
                    return true;
            }
            return false;
        }
        size_t strip(char ch = ' ') noexcept { // returns the number of characters removed
            if (is_empty || !contains(ch) || !ch)
                return 0;
            char *ptr = string;
            const char *from;
            size_t count = 0;
            while (*ptr) {
                if (*ptr == ch) {
                    for (from = ptr; *ptr;)
                        *ptr++ = *++from;
                    --length_w_null;
                    --space_back;
                    ptr = string;
                    ++count;
                    continue;
                }
                ++ptr;
            }
            if (length_w_null < 2)
                this->clear();
            return count;
        }
        size_t strip(const char *characters) noexcept {
            if (is_empty)
                return 0;
            size_t count = 0;
            while (*characters)
                count += strip(*characters++);
            return count;
        }
        size_t remove(char ch) noexcept {
            if (is_empty || !contains(ch) || !ch)
                return nopos;
            char *ptr = string;
            size_t pos = 0;
            while (true) {
                if (*ptr == ch) {
                    for (; *ptr; ++ptr)
                        *ptr = *(ptr + 1);
                    --space_back;
                    --length_w_null;
                    if (length_w_null == 1)
                        this->clear();
                    return pos;
                }
                ++ptr;
                ++pos;
            }
        }
        size_t remove(const char *str) {
            size_t len = strlen_c(str);
            if (is_empty || str == nullptr || !*str || length_w_null - 1 < len)
                return nopos;
            if (length_w_null - 1 == len) {
                if (!strcmp_c(string, str)) {
                    this->clear();
                    return 0;
                }
                return nopos;
            }
            char *optr = string;
            char *ptr = string;
            char *sub = new char[len + 1];
            char *optr_sub;
            size_t count;
            for (size_t pos = 0; pos < length_w_null - len; ++pos) {
                count = 0;
                optr_sub = optr;
                while (count < len) {
                    *sub = *optr_sub;
                    ++count; ++sub; ++optr_sub;
                }
                *sub = 0;
                sub -= len;
                if (!strcmp_c(str, sub)) {
                    for (size_t i = 0; i < len; ++i) {
                        ptr = optr;
                        for (; *ptr; ++ptr) {
                            *ptr = *(ptr + 1);
                        }
                        --space_back;
                        --length_w_null;
                    }
                    if (length_w_null == 1)
                        this->clear();
                    delete[] sub;
                    return pos;
                }
                ++optr;
                ++ptr;
            }
            delete[] sub;
            return nopos;
        }
        size_t r_remove(char ch) noexcept {
            if (is_empty || !contains(ch) || !ch)
                return nopos;
            char *ptr = string + length_w_null - 2;
            size_t pos = length_w_null - 2;
            while (true) {
                if (*ptr == ch) {
                    for (; *ptr; ++ptr)
                        *ptr = *(ptr + 1);
                    --space_back;
                    --length_w_null;
                    if (length_w_null == 1)
                        this->clear();
                    return pos;
                }
                --ptr;
                --pos;
            }
        }
        std::vector<String> split(char delim = 32) const {
            if (is_empty)
                return {};
            if (!contains(delim))
                return {string};
            std::vector<String> retvec;
            char *ptr = string;
            gtd::String part{true};
            while (*ptr) {
                while (*ptr != delim && *ptr) {
                    part.push_back(*ptr);
                    ++ptr;
                }
                if (!part.empty())
                    retvec.push_back(part);
                part.clear();
                ++ptr;
            }
            return retvec;
        }
        bool startswith(const char *beg) const noexcept {
            if (!contains(beg))
                return false;
            return !strncmp_c(string, beg, strlen_c(beg));
        }
        bool endswith(const char *end) const noexcept {
            if (!contains(end))
                return false;
            return !strcmp_c(string + length_w_null - strlen_c(end) - 1, end);
        }
        bool adopt_text(const char *path) noexcept {
            if (path == nullptr || !*path)
                return false;
            struct stat buffer{};
            FILE *fp;
            if (stat(path, &buffer) == -1 || S_ISDIR(buffer.st_mode) || !buffer.st_size ||
                (fp = fopen(path, "r")) == nullptr || fgetc(fp) == EOF)
                return false;
            shrunk = false;
            if (!is_empty)
                this->clear();
            fseek(fp, 0, SEEK_END);
            size_t count = ftell(fp);
            // while (fgetc(fp) != EOF) // instead of using struct stat, because st_size would include EOF & others
            //     ++count;
            fseek(fp, 0, SEEK_SET);
            length_w_null = count + 1;
            set_size(false);
            data = new char[size];
            string = data + get_first_pos();
            space_front = get_first_pos();
            space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front+1 : space_front+2);
            is_empty = false;
            memset_c(data, 0, size);
            fread(string, sizeof(char), count, fp);
            fclose(fp);
            return true;
        }
        size_t word_count() const noexcept { // NEEDS WORK
            if (is_empty)
                return 0;
            size_t word_count = 0;
            size_t word_len = 0;
            bool has_backspace = false;
            bool has_carriage_return = false;
            const char *ptr = string;
            while (*ptr) {
                if (*ptr >= 33 && *ptr <= 126) {
                    word_len = 1;
                    ++word_count;
                    while (*(++ptr) >= 33 && *ptr <= 126)
                        ++word_len;
                    if (!*ptr)
                        return word_count;
                    if (*ptr == '\b' && (*(ptr + 1) >= 33 && *(ptr + 1) <= 126))
                        --word_count;
                }
                ++ptr;
            }
        }
        size_t shift_center() noexcept {
            if (is_empty || space_back == space_front || space_back == space_front + 1)
                return 0;
            size_t move_by = space_back > space_front ? (is_even(length_w_null) ? (space_back - space_front) / 2 :
                    (space_back - space_front - 1) / 2) : (is_even(length_w_null) ? (space_front - space_back) / 2 :
                            (space_front - space_back + 1) / 2);
            char *first_ptr = string;
            char *fixed_first = first_ptr;
            char *last_ptr = string + length_w_null - 1;
            char *less_ptr = nullptr;
            char *fixed_last = last_ptr;
            if (space_back > space_front) {
                for (size_t i = 0; i < move_by; ++i) {
                    last_ptr = fixed_last;
                    less_ptr = last_ptr - 1;
                    for (; last_ptr > first_ptr; --last_ptr, --less_ptr) {
                        *last_ptr = *less_ptr;
                        *less_ptr = 0;
                    }
                    ++first_ptr;
                    ++fixed_last;
                }
                string += move_by;
                space_front += move_by;
                space_back -= move_by;
                return move_by;
            }
            for (size_t i = 0; i < move_by; ++i) {
                first_ptr = fixed_first;
                less_ptr = first_ptr - 1;
                for (; first_ptr <= last_ptr; ++first_ptr, ++less_ptr)
                    *less_ptr = *first_ptr;
                --fixed_first;
                --last_ptr;
            }
            string -= move_by;
            space_front -= move_by;
            space_back += move_by;
            return move_by;
        }
        size_t shift_left() noexcept {
            if (is_empty || !space_front)
                return 0;
            char *less_ptr = string - 1;
            char *fixed_ptr = string;
            char *ptr = string;
            const char *end = string + length_w_null;
            while (fixed_ptr > data) {
                while (ptr < end) {
                    *less_ptr = *ptr;
                    ++ptr; ++less_ptr;
                }
                --end;
                --fixed_ptr;
                ptr = fixed_ptr;
                less_ptr = fixed_ptr - 1;
            }
            size_t retval = space_front;
            space_front = 0;
            space_back = size - length_w_null;
            string = data;
            return retval;
        }
        size_t shift_right() noexcept {
            if (is_empty || !space_back)
                return 0;
            char *fixed_beg = string;
            char *ptr = string + length_w_null - 1;
            char *less_ptr = ptr - 1;
            char *fixed_end = ptr;
            const char *end = data + size - 1;
            while (fixed_end < end) {
                while (less_ptr >= fixed_beg) {
                    *ptr = *less_ptr;
                    --ptr; --less_ptr;
                }
                *fixed_beg = 0;
                ++fixed_beg; ++fixed_end;
                ptr = fixed_end; less_ptr = ptr - 1;
            }
            string += space_back;
            size_t retval = space_back;
            space_back = 0;
            space_front = size - length_w_null;
            return retval;
        }
        size_t transform_chars(const std::function<void(char&)> &func) { // returns the number of characters altered
            if (is_empty)
                return 0;
            char c;
            size_t transformed_cnt = 0;
            for (char &ch : *this) {
                c = ch;
                func(ch);
                if (!ch) {
                    ch = c;
                    continue;
                }
                if (c != ch)
                    ++transformed_cnt;
            }
            return transformed_cnt;
        }
        size_t transform_chars(const std::function<void(char&)> &&func) { // to allow lambda expressions
            return transform_chars(func);
        } // NEEDS LOTS OF WORK (below)
        size_t transform_words(const std::function<const char *(char*)> &func, bool check_punctuation = true) { // returns the number of words altered
            if (is_empty || !contains(32))
                return 0;
            char *start_ptr = string;
            char *ptr = string;
            char *start_of_word = nullptr;
            char *word = nullptr;
            const char *word_copy = nullptr;
            const char *new_word = nullptr;
            size_t word_len = 0;
            size_t words_altered = 0;
            bool end_punc = false;
            bool tripunc = false;
            bool altered = false;
            size_t num_puncs = 0;
            while (*ptr) {
                if (*ptr == 32) {
                    ++ptr;
                    continue;
                }
                end_punc = false;
                tripunc = false;
                altered = false;
                num_puncs = 0;
                start_ptr = ptr;
                start_of_word = ptr;
                word_len = 0;
                while (*ptr && *ptr != 32) {
                    if (is_punc(*ptr))
                        ++num_puncs;
                    ++word_len;
                    ++ptr;
                }
                if (check_punctuation) {
                    if (num_puncs > 0) {
                        if (num_puncs == 1) {
                            if (is_punc(*(ptr - 1)) && word_len > 1) {
                                end_punc = true;
                                --ptr;
                                --word_len;
                            }
                            else if (is_punc(*start_ptr) && word_len > 1) {
                                ++start_ptr; ++start_of_word;
                                --word_len;
                            }
                        }
                        else if (num_puncs == 2) {
                            if (is_punc(*(ptr - 1)) && is_punc(*start_ptr) && word_len > 2) {
                                end_punc = true;
                                --word_len; --word_len;
                                --ptr; ++start_ptr; ++start_of_word;
                            }
                        }
                        else if (num_puncs == 3) {
                            if (*(ptr - 1) == '.' && (*(ptr - 2) == '\'' || *(ptr - 2) == '"') &&
                                (*start_ptr == '\'' || *start_ptr == '"')) {
                                tripunc = true;
                                --word_len; --word_len; --word_len;
                                --ptr; --ptr; ++start_ptr; ++start_of_word;
                            }
                        }
                    }
                }
                word = new char[word_len + 1];
                word_copy = word;
                for (size_t i = 0; i < word_len; ++i, ++start_of_word)
                    *(word + i) = *start_of_word;
                *(word + word_len) = 0;
                new_word = func(word);
                if (!strcmp_c(new_word, word)) { // same word, so same length
                    if (new_word != word) {
                        delete[] word;
                        return nopos;
                    }
                    for (size_t i = 0; start_ptr < ptr; ++start_ptr, ++i) {
                        if (!altered && *start_ptr != *(word + i)) {
                            altered = true;
                            ++words_altered;
                        }
                        if (!*(word + i))
                            *(word + i) = 32;
                        *start_ptr = *(word + i);
                    }
                    if (end_punc)
                        ++ptr;
                    else if (tripunc) {
                        ++ptr;
                        ++ptr;
                    }
                }
                else {

                    size_t new_len = strlen_c(new_word);
                }
                delete[] word;
            }
            return words_altered;
        }
        size_t shrink_to_fit() {
            if (is_empty || (!space_front && !space_back))
                return 0;
            shrunk = true;
            const char *old_ptr = data;
            data = new char[length_w_null];
            strcpy_c(data, string);
            delete[] old_ptr;
            string = data;
            size_t retval = size - length_w_null;
            size = length_w_null;
            space_front = space_back = 0;
            return retval;
        }
        void space() const noexcept { // debugging method - will eventually be removed
            printf("Length w null: %zu, size: %zu, space front: %zu, space back: %zu, string - data: %ld\n",
                   length_w_null, size, space_front, space_back, string - data);
        }
        void clear() {
            delete[] data;
            empty_constructor();
        }
        //void erase() noexcept {
        //
        //}
        // erase_chars() shifts the string in whichever direction saves computation
        String &erase_chars(size_t start_pos = 0, size_t end_pos = nopos) noexcept { // erases excluding end_pos
            if (is_empty || start_pos >= length_w_null - 1 || start_pos >= end_pos)
                return *this;
            if (end_pos > length_w_null - 1)
                end_pos = length_w_null - 1;
            char *str = string + start_pos;
            char *end = string + end_pos - 1;
            char *end_end = string + length_w_null - 1;
            size_t diff = end_pos - start_pos;
            if (!start_pos) { // case for entire beginning-of-string being erased
                while (string <= end)
                    *string++ = 0;
                if (end_pos == length_w_null - 1) { // case for entire string being erased
                    is_empty = true;
                    string = nullptr;
                    if (start_left) {
                        space_front = 0;
                        space_back = size;
                    } else
                        space_front = space_back = size / 2;
                    length_w_null = 0;
                    return *this;
                }
                space_front += diff;
                length_w_null -= diff;
                return *this;
            }
            if (end_pos == length_w_null - 1) { // case for entire end-of-string being erased
                while (end >= str)
                    *end-- = 0;
                space_back += diff;
                length_w_null -= diff;
                return *this;
            }
            length_w_null -= diff;
            if (str - string >= end_end - end) { // case for end-of-string being shifted left
                ++end;
                while (end <= end_end)
                    *str++ = *end++;
                while (str <= end_end)
                    *str++ = 0;
                space_back += diff;
                return *this;
            } // below: case for end-of-string being shifted right
            --str;
            while (str >= string)
                *end-- = *str--;
            while (end >= string)
                *end-- = 0;
            string += diff;
            space_front += diff;
            return *this;
        }
        size_t find(char c) const noexcept {
            if (is_empty)
                return nopos;
            const char *str = string;
            const char *end = string + length_w_null - 1;
            size_t count = 0;
            while (str < end) {
                if (*str++ == c)
                    return count;
                ++count;
            }
            return nopos;
        }
        size_t r_find(char c) const noexcept {
            if (is_empty)
                return nopos;
            const char *str = string + length_w_null - 2;
            const char *end = string;
            size_t count = length_w_null - 2;
            while (str >= end) {
                if (*str-- == c)
                    return count;
                --count;
            }
            return nopos;
        }
        Iterator begin() noexcept {
            if (is_empty) return {};
            return {string};
        }
        Iterator end() noexcept {
            if (is_empty) return {};
            return {string + length_w_null - 1};
        }
        RevIterator rbegin() noexcept {
            if (is_empty) return {};
            return {string + length_w_null - 2};
        }
        RevIterator rend() noexcept {
            if (is_empty) return {};
            return {string - 1};
        }
        ConstIterator cbegin() const noexcept {
            if (is_empty) return {};
            return {string};
        }
        ConstIterator cend() const noexcept {
            if (is_empty) return {};
            return {string + length_w_null - 1};
        }
        ConstRevIterator crbegin() const noexcept {
            if (is_empty) return {};
            return {string + length_w_null - 2};
        }
        ConstRevIterator crend() const noexcept {
            if (is_empty) return {};
            return {string - 1};
        }
        char &front() const {
            if (is_empty)
                throw EmptyStringError("front() cannot be called on an empty string.\n");
            return *string;
        }
        char &back() const {
            if (is_empty)
                throw EmptyStringError("back() cannot be called on an empty string.\n");
            return *(string + length_w_null - 2);
        }
        size_t get_size() const noexcept {
            return size;
        }
        size_t get_length() const noexcept {
            return !length_w_null ? 0 : length_w_null - 1;
        }
        bool empty() const noexcept {
            return is_empty;
        }
        const char *c_str() const noexcept {
            if (is_empty)
                return nullptr; // perhaps change this
            return string;
        }
        std::string str() const {
            return {string};
        }
        ~String() {
            delete[] data;
        }
        char &operator[](size_t index) const {
            if (is_empty)
                throw std::out_of_range("No characters to be accessed in an empty string.\n");
            if (index > length_w_null - 2)
                throw std::out_of_range("You are indexing a character that is out of range.\n");
            return *(string + index);
        }
        String operator+(const char *str) {
            return String{*this}.append_back(str);
        }
        String operator+(const String &str) {
            return String{*this}.append_back(str);
        }
        String operator+(const std::string &str) {
            return String{*this}.append_back(str);
        }
        String &operator+=(const char *str) {
            return this->append_back(str);
        }
        String &operator+=(const std::string &str) {
            return this->append_back(str);
        }
        String &operator+=(const String &str) {
            return this->append_back(str);
        }
        String &operator=(const char *str) {
            if (str == nullptr)
                throw NullPointerError();
            if (*this == str)
                return *this;
            if (is_empty) {
                constructor(str, true);
                return *this;
            }
            this->clear();
            if (!strlen_c(str))
                return *this;
            this->constructor(str);
            return *this;
        }
        String &operator=(const std::string &str) {
            this->clear();
            if (str.empty())
                return *this;
            return *this = str.c_str();
        }
        String &operator=(const String &str) {
            if (this == &str)
                return *this;
            return *this = str.c_str();
        }
        String &operator=(String &&str) noexcept {
            if (this == &str)
                return *this;
            this->string = str.string;
            str.string = nullptr;
            this->data = str.data;
            str.data = nullptr;
            this->length_w_null = str.length_w_null;
            str.length_w_null = 0;
            this->is_empty = str.is_empty;
            str.is_empty = true;
            this->size = str.size;
            str.size = 0;
            this->space_front = str.space_front;
            str.space_front = 0;
            this->space_back = str.space_back;
            str.space_back = 0;
            this->start_left = str.start_left;
            str.start_left = false;
            this->shrunk = str.shrunk;
            str.shrunk = false;
            return *this;
        }
        String &operator++() {
            if (is_empty)
                return *this;
            for (char &ch : *this)
                if (ch < 127)
                    ++ch;
            return *this;
        }
        String operator++(int) {
            if (is_empty)
                return {};
            String new_str = *this;
            for (char &ch : *this)
                if (ch < 127)
                    ++ch;
            return new_str;
        }
        String &operator--() {
            if (is_empty)
                return *this;
            for (char &ch : *this)
                if (ch > 1)
                    --ch;
            return *this;
        }
        String operator--(int) {
            if (is_empty)
                return {};
            String new_str = *this;
            for (char &ch : *this)
                if (ch > 1)
                    --ch;
            return new_str;
        }
        explicit operator bool() const noexcept {
            return !this->is_empty;
        }
        friend bool operator==(const String&, const char*) noexcept;
        friend std::ostream &operator<<(std::ostream&, const String&);
        friend std::istream &operator>>(std::istream&, String&);
        friend std::istream &getline(std::istream&, String&, char);
    };
    bool operator==(const String &g_str, const char *str) noexcept {
        if (str == nullptr)
            return false;
        if (g_str.is_empty)
            return !*str;
        const char *cpy_str = g_str.string;
        size_t count = 0;
        while (*cpy_str && *str) {
            if (*cpy_str++ != *str++)
                return false;
            ++count;
        }
        return count == g_str.length_w_null - 1 && !*str;
    }
    bool operator==(const String &str1, const String &str2) {
        return str1 == str2.c_str();
    }
    bool operator==(const String &g_str, const std::string &str) {
        return g_str == str.c_str();
    }
    std::ostream &operator<<(std::ostream &os, const String &str) {
        return os << str.string;
    }
    std::istream &operator>>(std::istream &is, String &str) {
        if (!is.good())
            return is;
        if (is.peek() == EOF) {
            is.setstate(std::ios::failbit);
            return is;
        }
        bool org_s_left = str.start_left;
        if (!org_s_left)
            str.start_left = true;
        if (!str.is_empty)
            str.clear();
        int c;
        while ((c = is.get()) != 32 && c != '\n' && is.good())
            str.push_back((char) c);
        is.clear();
        str.start_left = org_s_left;
        return is;
    }
    std::istream &getline(std::istream &is, String &str, char delim = '\n') {
        if (!is.good() || is.peek() == EOF)
            return is;
        bool org_s_left = str.start_left;
        if (!org_s_left)
            str.start_left = true;
        if (!str.is_empty)
            str.clear();
        int ch;
        while ((ch = is.get()) != delim && is.good())
            str.push_back((char) ch);
        str.start_left = org_s_left;
        return is;
    }
    // called for Iterator, RevIterator, ConstIterator and ConstRevIterator
    bool operator==(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr == B.ptr;
    }
    bool operator!=(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr != B.ptr;
    }
    // called for Iterator and ConstIterator
    bool operator<(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr < B.ptr;
    }
    bool operator>(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr > B.ptr;
    }
    bool operator<=(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr <= B.ptr;
    }
    bool operator>=(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr >= B.ptr;
    }
    std::ptrdiff_t operator-(const String::Iterator &A, const String::Iterator &B) {
        return A.ptr - B.ptr;
    }
    std::ostream &operator<<(std::ostream &out, const String::Iterator &A) {
        return out << static_cast<void *>(A.ptr);
    }
    // std::ostream &operator>>(const String &str, std::ostream &os) {
    //     if (str.is_empty)
    //         return os;
    //     return os << str.string;
    // }
    // called for RevIterator and ConstRevIterator
    bool operator<(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr > B.ptr;
    }
    bool operator>(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr < B.ptr;
    }
    bool operator<=(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr >= B.ptr;
    }
    bool operator>=(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return false;
        return A.ptr <= B.ptr;
    }
    std::ptrdiff_t operator-(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr)
            return 0;
        return B.ptr - A.ptr;
    }
    bool operator==(const char *str, const String &string) {
        return string == str;
    }
    bool operator==(const std::string &str, const String &string) {
        return string == str;
    }
    String operator+(const char *str, const String &string) {
        return String{string}.append_front(str);
    }
    String operator+(const std::string &str, const String &string) {
        return String{string}.append_front(str);
    }
    size_t strlen_c(const char *str) {
        if (str == nullptr)
            return 0;
        size_t length_c = 0;
        while (*str++) ++length_c;
        return length_c;
    }
    char *memset_c(char *str, char ch, size_t n_chars) {
        if (str == nullptr)
            return nullptr;
        for (size_t i = 0; i < n_chars; i++)
            *(str + i) = ch;
        return str;
    }
    char *strcpy_c(char *dest, const char *source) {
        if (dest == nullptr || source == nullptr)
            return nullptr;
        char *ptr = dest;
        while ((*dest++ = *source++));
        return ptr;
    }
    char *strcat_c(char *_dest, const char *_src) {
        if (_dest == nullptr || _src == nullptr)
            return nullptr;
        char *ptr = _dest;
        if (!*ptr) goto after;
        while (*++ptr);
        after:
        while ((*ptr++ = *_src++));
        return _dest;
    }
    int strcmp_c(const char *str1, const char *str2) {
        if (str1 == nullptr || str2 == nullptr)
            return -128;
        while (*str1 && *str2) {
            if (*str1 != *str2)
                return *str1 - *str2;
            ++str1;
            ++str2;
        }
        return *str1 - *str2;
    }
    int strncmp_c(const char *str1, const char *str2, size_t n) {
        if (str1 == nullptr || str2 == nullptr)
            return -128;
        size_t count = 0;
        while (*str1 && *str2 && count++ < n) {
            if (*str1 != *str2)
                return *str1 - *str2;
            ++str1;
            ++str2;
            // ++count;
        }
        return 0;
    }
    inline char *setchr_c(char *str, const char ch, size_t pos) {
        if (str == nullptr)
            return nullptr;
        *(str + pos) = ch;
        return str;
    }
    inline unsigned char to_upper(unsigned char ch) {
        return ch <= 122 && ch >= 97 ? ch - 32 : ch;
    }
    inline unsigned char to_lower(unsigned char ch) {
        return ch <= 90 && ch >= 65 ? ch + 32 : ch;
    }
    std::string to_upper(const char *str) {
        size_t length = strlen_c(str);
        std::string ret_string;
        for (size_t i = 0; i < length; i++)
            ret_string.push_back((char) to_upper(*(str + i)));
        return ret_string;
    }
    std::string to_upper(const std::string &str) {
        std::string ret_string;
        for (const char &ch: str)
            ret_string.push_back((char) to_upper(ch));
        return ret_string;
    }
    void char_upper(char &c) {
        if (c >= 97 && c <= 122)
            c -= 32;
    }
    void char_lower(char &c) {
        if (c >= 65 && c <= 90)
            c += 32;
    }
    void string_upper(std::string &str) {
        std::for_each(str.begin(), str.end(), char_upper);
    }
    void string_upper(char *str) {
        if (str == nullptr)
            return;
        while (*str) {
            *str = (char) to_upper(*str);
            ++str;
        }
    }
    void string_lower(std::string &str) {
        std::for_each(str.begin(), str.end(), char_lower);
    }
    void string_lower(char *str) {
        if (str == nullptr)
            return;
        for (; *str; ++str)
            *str = (char) to_lower(*str);
    }
    inline bool isdigit_c(char ch) {
        return ch >= 48 && ch <= 57;
    }
    template <typename T = uint64_t> requires (std::integral<T>)
    T to_int(const char *str) {
        if (!str || !*str)
            throw std::invalid_argument{"Error: null or empty string encountered.\n"};
        T retval{};
        int mul_by;
        if (*str == '-') {
            mul_by = -1;
            if (!*++str)
                throw std::invalid_argument{"Error: a negative sign on its own is not a valid number.\n"};
        } else
            mul_by = 1;
        goto start_loop;
        while (*str) {
            retval *= 10;
            start_loop:
            if (!isdigit_c(*str))
                throw std::invalid_argument{"Error: non-numeric character encountered.\n"};
            retval += *str++ - 48;
        }
        return mul_by*retval;
    }
    bool is_numeric(const char *str) {
        if (str == nullptr || !*str)
            return false;
        while (*str)
            if (!isdigit_c(*str++))
                return false;
        return true;
    }
    bool is_numeric(const std::string &str) {
        if (str.empty())
            return false;
        for (const char &ch: str)
            if (!isdigit_c(ch))
                return false;
        return true;
    }
    bool contains(const char *str, char ch) {
        if (str == nullptr || !*str || !ch)
            return false;
        while (*str)
            if (*str++ == ch)
                return true;
        return false;
    }
    size_t count_char(const char *str, char ch = 32) {
        if (str == nullptr || !*str)
            return 0;
        size_t num = 0;
        while (*str)
            if (*str++ == ch)
                ++num;
        return num;
    }
    char **strsplit(const char *str, char delim = 32) {
        if (str == nullptr || !*str)
            return nullptr;
        char **retptr;
        if (!delim || !contains(str, delim)) {
            retptr = new char*[2*sizeof(char *)];
            *(retptr + 1) = nullptr;
            *retptr = new char[strlen_c(str) + 1];
            strcpy_c(*retptr, str);
            return retptr;
        }
        size_t str_len = 0;
        size_t str_count = 0;
        if (*str == delim)
            ++str;
        const char *beg = str;
        size_t num_strings = count_char(str, delim) + 1;
        retptr = new char*[num_strings + 1];
        while (str_count != num_strings) {
            if (*str == delim || *str == 0) {
                *(retptr + str_count) = new char[str_len + 1];
                for (size_t i = 0; i < str_len; ++i, ++beg)
                    *(*(retptr + str_count) + i) = *beg;
                *(*(retptr + str_count) + str_len) = 0;
                ++beg;
                ++str_count;
                str_len = 0;
                ++str;
                continue;
            }
            ++str;
            ++str_len;
        }
        *(retptr + num_strings) = nullptr;
        return retptr;
    }
    void free_split_str(const char *const *array) {
        if (array == nullptr || *array == nullptr)
            return;
        const char *const *ptr = array;
        while (*ptr != nullptr) {
            delete[] *ptr;
            ++ptr;
        }
        delete[] *ptr;
        delete[] array;
    }
    inline void clear_cin() {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // discards characters until newline found
    }
    template<typename charT = char*>
    charT get_home_path() = delete;
    template<> char *get_home_path<char *>() {
#ifdef _WIN32
        static char home_path_c[MAX_PATH];
        if (SHGetFolderPathA(nullptr, CSIDL_PROFILE, nullptr, SHGFP_TYPE_CURRENT, home_path_c) != S_OK)
            return nullptr;
        return home_path_c;
#else
        struct passwd *pwd = getpwuid(getuid());
        if (pwd == nullptr)
            return nullptr;
        return pwd->pw_dir; // the passwd struct is statically allocated by getpwuid(), so this is safe
#endif
    }
    template<> std::string get_home_path<std::string>() { // returns empty string in case of error
#ifdef _WIN32
        std::string path{MAX_PATH, 0};
        if (SHGetFolderPathA(nullptr, CSIDL_PROFILE, nullptr, SHGFP_TYPE_CURRENT, path.data()) != S_OK)
            return path.erase();
        return path;
#else
        struct passwd *pwd = getpwuid(getuid());
        if (pwd == nullptr)
            return {};
        return {pwd->pw_dir};
#endif
    }
    template<> String get_home_path<String>() {
#ifdef _WIN32
        static char home_path_c[MAX_PATH];
        if (SHGetFolderPathA(nullptr, CSIDL_PROFILE, nullptr, SHGFP_TYPE_CURRENT, home_path_c) != S_OK)
            return {};
        return {home_path_c};
#else
        struct passwd *pwd = getpwuid(getuid());
        if (pwd == nullptr)
            return {};
        return {pwd->pw_dir};
#endif
    }
    consteval char file_sep() {
#ifdef _WIN32
        return '\\';
#else
        return '/';
#endif
    }
    const char *get_date_and_time() { // returns date & time with '_' instead of spaces and 'h', 'm', 's' inst. of colon
        time_t t = time(nullptr);
        char *ptr = ctime(&t);
        char *end = ptr + 24; // ctime() always returns a string with 26 chars (including '\0'), end points to '\n'
        for (unsigned char i = 0; i < 4; ++i) {
            *end = *(end - 1);
            --end;
        }
        *end-- = '_';
        *end = 's';
        *(end -= 3) = 'm'; // easy to avoid a loop here
        *(end -= 3) = 'h';
        *(end -= 3) = '_';
        *(end -= 3) = '_';
        *(end - 4) = '_';
        if (!isdigit_c(*++end)) // case for first 9 days of month
            *end = 48;
        return ptr;
    }
    template<typename whatever>
    void print(whatever w) {
        std::cout << w << std::endl;
    }
    void print(bool w) {
        std::cout << std::boolalpha << w << std::endl;
    }
    template <typename T, typename... types>
    void print(T arg, types... args) {
        print(arg);
        if constexpr (sizeof...(args))
            print(args...);
    }
    namespace word_transforms {
        typedef const char *(*WT)(char *);
        const char *reverse(char *word) {
            size_t length = strlen_c(word);
            if (word == nullptr || !length || length == 1)
                return word;
            char *beg = word;
            char *end = word + strlen_c(word) - 1;
            char temp;
            do {
                temp = *beg;
                *beg = *end;
                *end = temp;
            } while (++beg < --end);
            return word;
        }
        const char *cap_first_letter(char *word) {
            if (word == nullptr || !*word)
                return word;
            *word = *word <= 122 && *word >= 97 ? *word - 32 : *word;
            return word;
        }
    }
}
#endif
