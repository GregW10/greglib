#ifndef GREGCSV_HPP
#define GREGCSV_HPP

#include "ml/gregvct.hpp"

namespace gtd {
    class csv_format_error : public std::ios_base::failure {
    public:
        explicit csv_format_error(const char *msg) : std::ios_base::failure{msg} {}
        csv_format_error() : std::ios_base::failure{"Error: invalid CSV format.\n"} {}
    };
    template <gml::Numeric IntType, gml::Numeric FloatType>
    class csv_col {
        char *_name{};
        std::vector<char*> *_s{};
        std::vector<IntType> *_i{};
        std::vector<FloatType> *_f{};
        csv_col() = default;
    public:
        explicit csv_col(const char *name, uint64_t len = 0) {
            if (!name || !*name)
                throw std::invalid_argument{"Error: a name must be passed.\n"};
            if (len) {
                _name = new char[len + 1];
                char *ptr = _name;
                while (len --> 0)
                    *ptr++ = *name++;
                *ptr = 0;
                return;
            }
            _name = new char[gml::gen::strlen_c(name) + 1];
            gml::gen::strcpy_c(_name, name);
        }
        ~csv_col() {
            delete [] _name;
        }
        template <gml::Numeric, gml::Numeric>
        friend class csv;
    };
    template <gml::Numeric IntType, gml::Numeric FloatType>
    class csv {
        std::vector<csv_col<IntType, FloatType>> _cols;
        void load_csv(const char *path) {
            if (!path || !*path)
                throw std::invalid_argument{"Error: nullptr or empty path.\n"};
            std::ifstream file{path, std::ios_base::in};
            if (!file.good())
                throw std::ios_base::failure{"Error: could not open file.\n"};
            std::string line;
            std::getline(file, line);
            if (line.empty())
                throw csv_format_error{"Error: empty CSV file -> nothing to load.\n"};
            // if (line.front() == '\n')
            //     throw csv_format_error{"Error: BLAH.\n"};
            const char *sptr = line.data();
            const char *eptr = line.data();
            uint64_t len = 0;
            uint64_t nempty = 0;
            uint64_t nreps;
            std::vector<std::string> cnames;
            std::string cname;
            std::vector<std::string>::size_type pos = 0;
            std::vector<std::string>::iterator it;
            std::vector<std::string>::iterator end;
            bool oot = false;
            do {
                if (*eptr != ',') {
                    if (*eptr == '\n') {
                        oot = true;
                        goto _rest;
                    }
                    ++eptr;
                    ++len;
                    continue;
                }
                _rest:
                if (!len) { // case for empty column
                    pos = 0;
                    cname = "NONAME_";
                    _start:
                    cname += std::to_string(++nempty);
                    it = cnames.begin() + pos;
                    end = cnames.end();
                    while (it != end) {
                        ++pos;
                        if (cname == *it++) {
                            cname.erase(cname.begin() + 7, cname.end());
                            goto _start;
                        }
                    }
                    cnames.push_back(cname);
                } else {
                    cnames.emplace_back(sptr, len);
                    len = 0;
                }
                if (oot)
                    break;
                ++eptr;
                sptr = eptr;
            } while (1);
            while (std::getline(file, line)) {

            }
            file.close();
        }
    public:
        explicit csv(const char *path) {}
    };
}

#endif
