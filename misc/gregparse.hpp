#ifndef GREGPARSE_HPP
#define GREGPARSE_HPP

#include <vector>
#include <algorithm>
#include <list>
#include <sstream>
#include <string>
#include <regex>
#include "gregmisc.hpp"

namespace gtd {
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
        int argc; // mark for removal?
        const char *const *argv; // mark for removal? Maybe use to return pointer to arg, rather than reallocated arg
        std::list<std::pair<int, std::string>> cmdl{};
    public:
        enum dup_pol {
            THROW = 1,
            USE_FIRST = 2,
            USE_LAST = 4,
        };
        parser(int _argc, const char *const *_argv, int to_skip = 1) {
            if (_argc < 0) {
                err:
                throw std::invalid_argument{"Error: number of arguments cannot be negative.\n"};
            }
            if (to_skip > _argc)
                throw std::invalid_argument{"Error: skip number too large - cannot skip more arguments than "
                                            "those present.\n"};
            while (to_skip --> 0) {
                --_argc;
                ++_argv;
            }
            this->argc = _argc;
            this->argv = _argv;
            _argc = 0; // use as counter
            while (_argc < this->argc)
                this->cmdl.emplace_back(_argc++, *_argv++);
            // this->cmdl.assign(_argv, _argv + _argc);
        }
    private:
        std::pair<std::list<std::pair<int, std::string>>::iterator, bool> process_flag(const std::string &flag,
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
                    if (it->second == flag) {
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
                    if (it->second == flag) {
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
                    if (it->second == flag) {
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
                    std::string retstr{std::move(arg_it->second)};
                    this->cmdl.erase(arg_it);
                    return retstr;
                }
            }
            else {
                if (have) {
                    T retval;
                    if (!(std::istringstream{arg_it->second} >> retval)) // test for conversion error
                        throw arg_conversion_error{arg_it->second};
                    this->cmdl.erase(arg_it);
                    return retval;
                }
            }
            return def_val;
        }
        // template <>
        // bool get_arg<bool>(const std::string &flag, const bool &def_val, dup_pol duplication_policy);
        template <typename T>
        T get_arg(const std::string &flag, T (*func)(const char*), const T &def_val,
                  dup_pol duplication_policy = USE_FIRST) { // convenience method
            const char *arg = get_arg(flag, duplication_policy);
            return arg ? func(arg) : def_val;
        }
        const char *get_arg(const std::string &flag, dup_pol duplication_policy = USE_FIRST) {
            auto [arg_it, have] = process_flag(flag, duplication_policy);
            if (!have)
                return nullptr;
            int index = arg_it->first;
            this->cmdl.erase(arg_it);
            return *(this->argv + index);
        }
        const char *get_arg(const std::regex &rgx, dup_pol duplication_policy = USE_FIRST) {
            bool have = false;
            if (duplication_policy == THROW) {
                auto it = this->cmdl.begin();
                auto eit = this->cmdl.end();
                decltype(it) arg_it;
                while (it != eit) {
                    if (std::regex_match(it->second, rgx)) {
                        if (have)
                            throw duplicate_error{};
                        arg_it = it;
                        // it = this->cmdl.erase(it);
                        have = true;
                        // continue;
                    }
                    ++it;
                }
                if (have) {
                    int index = arg_it->first;
                    this->cmdl.erase(arg_it);
                    return *(this->argv + index);
                }
                return nullptr;
            }
            if (duplication_policy == USE_FIRST) {
                auto it = this->cmdl.begin();
                auto eit = this->cmdl.end();
                // decltype(it) arg_it;
                while (it != eit) {
                    if (std::regex_match(it->second, rgx)) {
                        int index = it->first;
                        this->cmdl.erase(it);
                        return *(this->argv + index);
                    }
                    ++it;
                }
                return nullptr;
            }
            if (duplication_policy == USE_LAST) {
                auto rit = this->cmdl.rbegin();
                auto reit = this->cmdl.rend();
                while (rit != reit) {
                    if (std::regex_match(rit->second, rgx)) {
                        int index = rit->first;
                        this->cmdl.erase(--rit.base());
                        return *(this->argv + index);
                    }
                    ++rit;
                }
            }
            return nullptr;
        }
        std::list<std::pair<int, std::string>>::size_type remaining() const noexcept {
            return this->cmdl.size();
        }
        bool empty() const noexcept {
            return this->cmdl.empty();
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
    template <>
    bool parser::get_arg<bool>(const std::string &flag, const bool &def_val, dup_pol duplication_policy) {
        size_t flen = flag.length();
        if (!flen)
            throw std::invalid_argument{"Error: an empty string cannot be passed as a flag.\n"};
        if (flen == 1)
            throw std::invalid_argument{"Error: a 1-character string cannot be passed as a flag.\n"};
        if (flag[0] != '-')
            throw std::invalid_argument{"Error: a flag MUST begin with a hyphen.\n"};
        bool dub = false;
        if (flag[1] == '-') {
            if (flen == 2)
                throw std::invalid_argument{"Error: a double-hyphen flag must have at least 3 characters.\n"};
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
                    if (it->second == flag) {
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
                    if (it->second == flag) {
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
                    if (it->second.operator[](0) == '-' && it->second.operator[](1) != '-') {
                        if (it->second.length() == 2 && it->second.operator[](1) == flag[1]) {
                            if (have)
                                throw duplicate_error{};
                            it = this->cmdl.erase(it);
                            have = true;
                            continue;
                        }
                        if ((pos = it->second.find(flag[1])) != std::string::npos) {
                            if (have)
                                throw duplicate_error{};
                            it->second.erase(pos, 1);
                            have = true;
                            if (it->second.find(flag[1]) != std::string::npos) // case in which there is still an
                                throw duplicate_error{}; // identical flag character remaining in the argument
                        }
                    }
                    ++it;
                }
            }
            else {
                while (it != eit) {
                    if (it->second.operator[](0) == '-' && it->second.operator[](1) != '-') {
                        if (it->second.length() == 2 && it->second.operator[](1) == flag[1]) {
                            it = this->cmdl.erase(it);
                            have = true;
                            continue;
                        }
                        if ((pos = it->second.find(flag[1])) != std::string::npos) {
                            it->second.erase(pos, 1);
                            have = true;
                            while ((pos = it->second.find(flag[1])) != std::string::npos)
                                it->second.erase(pos, 1);
                            if (it->second.length() == 1) {
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
}
#endif
