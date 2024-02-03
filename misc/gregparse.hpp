#ifndef GREGPARSE_HPP
#define GREGPARSE_HPP

#include <vector>
#include <algorithm>
#include <list>
#include <sstream>
#include <string>
#include <regex>

namespace gtd {
    bool str_eq(const char *s1, const char *s2) {
        if (!s1 || !s2)
            return false;
        while (*s1 || *s2)
            if (*s1++ != *s2++)
                return false;
        return true;
    }
#ifndef GREGALG_HPP
    bool memcopy(void *dst, const void *src, size_t bytes) {
        if (!bytes || !dst || !src)
            return false;
        char *cdst = (char*) dst;
        const char *csrc = (char*) src;
        while (bytes --> 0)
            *cdst++ = *csrc++;
        return true;
    }
#endif
#ifndef GREGSTR_HPP
    size_t strlen_c(const char *str) {
        if (!str || !*str)
            return -1;
        size_t len = 0;
        while (*str++) ++len;
        return len;
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
#endif
    class duplicate_error : public std::invalid_argument {
    public:
        duplicate_error() : std::invalid_argument{"Error: duplicate flags/arguments found.\n"} {}
        duplicate_error(const char *msg) : std::invalid_argument{msg} {}
    };
    class missing_arg_error : public std::invalid_argument {
    public:
        missing_arg_error() : std::invalid_argument{"Error: insufficient arguments for flags.\n"} {}
        missing_arg_error(const char *msg) : std::invalid_argument{msg} {}
    };
    class arg_conversion_error : public std::invalid_argument {
        std::string _arg;
    public:
        arg_conversion_error(const std::string &arg) :
        std::invalid_argument{"Error: the argument \"" + arg + "\" could not be converted to the required type.\n"},
        _arg{arg} {}
        const std::string &get_arg() const noexcept {
            return this->_arg; // this allows the user to access the problem argument within a catch block, if desired
        }
    };
    class parser {
        int argc;
        const char *const *argv;
        std::list<std::string> cmdl{};
    public:
        enum dup_pol {
            THROW = 1,
            USE_FIRST = 2,
            USE_LAST = 4,
        };
        parser(int _argc, const char *const *_argv, bool skip_first = true) {
            if (_argc < 0) {
                err:
                throw std::invalid_argument{"Error: number of arguments cannot be negative.\n"};
            }
            if (skip_first) {
                if (_argc == 0)
                    goto err;
                --_argc;
                ++_argv;
            }
            this->argc = _argc;
            this->argv = _argv;
            this->cmdl.assign(_argv, _argv + _argc);
        }
    private:
        std::pair<std::list<std::string>::iterator, bool> process_flag(const std::string &flag,
                                                                       dup_pol duplication_policy = USE_FIRST) {
            size_t flen = flag.length();
            if (!flen)
                throw std::invalid_argument{"Error: an empty string cannot be passed as a flag.\n"};
            if (flen == 1)
                throw std::invalid_argument{"Error: a 1-character string cannot be passed as a flag.\n"};
            if (flag[0] != '-')
                throw std::invalid_argument{"Error: a flag MUST begin with a hyphen.\n"};
            if (flag[1] == '-') {
                if (flen == 2)
                    throw std::invalid_argument{"Error: a double-hyphen flag must have at least 3 characters.\n"};
            }
            else {
                if (flen > 2)
                    throw std::invalid_argument{"Error: single-hyphen flag can only be two characters long.\n"};
            }
            if (duplication_policy != THROW && duplication_policy != USE_FIRST && duplication_policy != USE_LAST)
                throw std::invalid_argument{"Error: invalid duplication policy.\n"};
            auto it = this->cmdl.begin();
            auto eit = this->cmdl.end();
            bool have = false;
            decltype(it) arg_it;
            if (duplication_policy == THROW) {
                while (it != eit) {
                    if (*it == flag) {
                        if (have)
                            throw duplicate_error{};
                        it = this->cmdl.erase(it);
                        if (it == eit)
                            throw missing_arg_error{};
                        else
                            arg_it = it;
                        have = true;
                    }
                    ++it;
                }
            }
            else if (duplication_policy == USE_FIRST) {
                while (it != eit) {
                    if (*it == flag) {
                        it = this->cmdl.erase(it);
                        if (it == eit)
                            throw missing_arg_error{};
                        else if (have)
                            it = this->cmdl.erase(it);
                        else
                            arg_it = it;
                        have = true;
                    }
                    ++it;
                }
            }
            else {
                while (it != eit) {
                    if (*it == flag) {
                        it = this->cmdl.erase(it);
                        if (it == eit)
                            throw missing_arg_error{};
                        else if (have)
                            this->cmdl.erase(arg_it);
                        arg_it = it;
                        have = true;
                    }
                    ++it;
                }
            }
            return {arg_it, have};
        }
    public:
        template <typename T> requires requires (T arg, std::istringstream stream) {
            stream >> arg;
        }
        T get_arg(const std::string &flag, const T &def_val = T{}, dup_pol duplication_policy = USE_FIRST) {
            auto [arg_it, have] = process_flag(flag, duplication_policy);
            if constexpr (std::same_as<T, std::string>) {
                if (have) {
                    std::string retstr{std::move(*arg_it)};
                    this->cmdl.erase(arg_it);
                    return retstr;
                }
            }
            if (have) {
                T retval;
                if (!(std::istringstream{*arg_it} >> retval)) // test for conversion error
                    throw arg_conversion_error{*arg_it};
                this->cmdl.erase(arg_it);
                return retval;
            }
            return def_val;
        }
        template <>
        bool get_arg<bool>(const std::string &flag, const bool &def_val, dup_pol duplication_policy) {
            size_t flen = flag.length();
            if (!flen)
                throw std::invalid_argument{"Error: an empty string cannot be passed as a flag.\n"};
            if (flen == 1)
                throw std::invalid_argument{"Error: a 1-character string cannot be passed as a flag.\n"};
            if (flag[0] != '-')
                throw std::invalid_argument{"Error: a flag MUST begin with a hyphen.\n"};
            bool dub = false;
            if (flag[1] == '-') {
                if (flen == 2) {
                    throw std::invalid_argument{"Error: a double-hyphen flag must have at least 3 characters.\n"};
                }
                dub = true;
            } else {
                if (flen > 2)
                    throw std::invalid_argument{"Error: single-hyphen flag can only be two characters long.\n"};
            }
            if (duplication_policy != THROW && duplication_policy != USE_FIRST && duplication_policy != USE_LAST)
                throw std::invalid_argument{"Error: invalid duplication policy.\n"};
            auto it = this->cmdl.begin();
            auto eit = this->cmdl.end();
            bool have = false;
            if (dub) {
                if (duplication_policy == THROW) {
                    while (it != eit) {
                        if (*it == flag) {
                            if (have)
                                throw duplicate_error{};
                            it = this->cmdl.erase(it);
                            have = true;
                            continue;
                        }
                        ++it;
                    }
                }
                else {
                    while (it != eit) {
                        if (*it == flag) {
                            it = this->cmdl.erase(it);
                            have = true;
                            continue;
                        }
                        ++it;
                    }
                }
            }
            else {
                std::string::size_type pos;
                if (duplication_policy == THROW) {
                    while (it != eit) {
                        if (it->operator[](0) == '-') {
                            if (it->length() == 2 && it->operator[](1) == flag[1]) {
                                if (have)
                                    throw duplicate_error{};
                                it = this->cmdl.erase(it);
                                have = true;
                                continue;
                            }
                            if ((pos = it->find(flag[1])) != std::string::npos) {
                                if (have)
                                    throw duplicate_error{};
                                it->erase(pos, 1);
                                have = true;
                                if (it->find(flag[1]) != std::string::npos) // case in which there is still an
                                    throw duplicate_error{}; // identical flag character remaining in the argument
                            }
                        }
                        ++it;
                    }
                }
                else {
                    while (it != eit) {
                        if (it->operator[](0) == '-') {
                            if (it->length() == 2 && it->operator[](1) == flag[1]) {
                                it = this->cmdl.erase(it);
                                have = true;
                                continue;
                            }
                            if ((pos = it->find(flag[1])) != std::string::npos) {
                                it->erase(pos, 1);
                                have = true;
                                while ((pos = it->find(flag[1])) != std::string::npos)
                                    it->erase(pos, 1);
                                if (it->length() == 1) {
                                    it = this->cmdl.erase(it);
                                    continue;
                                }
                            }
                        }
                        ++it;
                    }
                }
            }
            return have == !def_val;
        }
        template <typename T>
        T get_arg(const std::string &flag, T (*func)(const char*), const T &def_val,
                  dup_pol duplication_policy = USE_FIRST) { // convenience method
            char *arg = get_arg(flag, duplication_policy);
            T retval = arg ? func(arg) : def_val;
            delete [] arg;
            return retval;
        }
        char *get_arg(const std::string &flag, dup_pol duplication_policy = USE_FIRST) {
            auto [arg_it, have] = process_flag(flag, duplication_policy);
            char *ptr = strcpy_c(new char[arg_it->length() + 1], arg_it->data());
            return have ? this->cmdl.erase(arg_it), ptr : nullptr;
        }
        char *get_arg(const std::regex &rgx, dup_pol duplication_policy = USE_FIRST) {
            bool have = false;
            if (duplication_policy == THROW) {
                auto it = this->cmdl.begin();
                auto eit = this->cmdl.end();
                decltype(it) arg_it;
                while (it != eit) {
                    if (std::regex_match(*it, rgx)) {
                        if (have)
                            throw duplicate_error{};
                        arg_it = it;
                        // it = this->cmdl.erase(it);
                        have = true;
                        // continue;
                    }
                    ++it;
                }
                char *ptr = strcpy_c(new char[arg_it->length() + 1], arg_it->data());
                return have ? this->cmdl.erase(arg_it), ptr : nullptr;
            }
            if (duplication_policy == USE_FIRST) {
                auto it = this->cmdl.begin();
                auto eit = this->cmdl.end();
                decltype(it) arg_it;
                while (it != eit) {
                    if (std::regex_match(*it, rgx)) {
                        char *ptr = strcpy_c(new char[it->length() + 1], it->data());
                        this->cmdl.erase(it);
                        return ptr;
                    }
                }
                return nullptr;
            }
            if (duplication_policy == USE_LAST) {
                auto rit = this->cmdl.rbegin();
                auto reit = this->cmdl.rend();
                while (rit != reit) {
                    if (std::regex_match(*rit, rgx)) {
                        char *ptr = strcpy_c(new char[rit->length() + 1], rit->data());
                        this->cmdl.erase(--rit.base());
                        return ptr;
                    }
                }
            }
            return nullptr;
        }
        bool empty() const noexcept {
            return this->argc;
        }
        auto begin() noexcept {
            return this->cmdl.begin();
        }
        auto end() noexcept {
            return this->cmdl.end();
        }
        auto cbegin() const noexcept {
            return this->cmdl.cbegin();
        }
        auto cend() const noexcept {
            return this->cmdl.cend();
        }
        ~parser() {
            // std::for_each(this->flags.begin(), this->flags.end(), [](const char *flag){delete [] flag;});
        }
    };
}
#endif
