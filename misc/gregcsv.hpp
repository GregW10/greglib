#ifndef GREGCSV_HPP
#define GREGCSV_HPP

#include "ml/gregvct.hpp"
#include "gregmisc.hpp"
#include "gregmmapper.hpp"
#include <unistd.h>
#include <fcntl.h>
#include <cinttypes>

#define BUFFER_SIZE 8'192

namespace gtd {
    class csv_format_error : public std::ios_base::failure {
    public:
        explicit csv_format_error(const char *msg) : std::ios_base::failure{msg} {}
        csv_format_error() : std::ios_base::failure{"Error: invalid CSV format.\n"} {}
    };
    template <gml::Numeric IntType, gml::Numeric FloatType>
    class csv_col {
        using I = IntType;
        using F = FloatType;
        using VC = std::vector<char*>;
        using VI = std::vector<I>;
        using VF = std::vector<F>;
        char *_name{};
        VC *_s = new VC{};
        VI *_i{};
        VF *_f{};
        bool alloced = false;
    public:
        csv_col() = default;
        explicit csv_col(uint64_t index) {
            uint64_t num_digs = ((uint64_t) log10l((long double) index)) + 1;
            _name = new char[num_digs + 1];
            alloced = true;
            snprintf(_name, num_digs + 1, "%" PRIu64, index);
        }
        explicit csv_col(const char *name, uint64_t len = 0) {
            if (!name || !*name)
                throw std::invalid_argument{"Error: a name must be passed.\n"};
            if (len) {
                _name = new char[len + 1];
                char *ptr = _name;
                while (len --> 0)
                    *ptr++ = *name++;
                *ptr = 0;
            } else {
                _name = new char[gml::gen::strlen_c(name) + 1];
                gml::gen::strcpy_c(_name, name);
            }
            alloced = true;
            this->_s = new std::vector<char*>{};
        }
        ~csv_col() {
            if (alloced)
                delete [] _name;
            if (this->_s)
                delete this->_s;
            if (this->_i)
                delete this->_i;
            if (this->_f)
                delete this->_f;
        }
        template <gml::Numeric, gml::Numeric>
        friend class csv;
    };
    template <gml::Numeric IntType, gml::Numeric FloatType>
    class csv {
        using I = IntType;
        using F = FloatType;
        std::vector<csv_col<I, F>> _cols;
        gtd::mmapper _smap{};
        char *_sdata{};
        uint64_t lf = 0; // longest field
        uint64_t nf = 0; // number of fields
        uint64_t nl = 0; // number of lines (including header)
        void load_csv(const char *path, bool header) {
            if (!path || !*path)
                throw std::invalid_argument{"Error: nullptr or empty path.\n"};
            struct stat buff{};
            if (stat(path, &buff) == -1)
                throw std::ios_base::failure{"Error: could not obtain file information.\n"};
            if (!buff.st_size)
                throw csv_format_error{"Error: invalid CSV format - file size is zero.\n"};
            int fd = open(path, O_RDONLY);
            if (fd == -1)
                throw std::ios_base::failure{"Error: couldn't open CSV file.\n"};
            _sdata = (char*) _smap.reset(buff.st_size);
            if (gtd::read_all(fd, _sdata, buff.st_size) != buff.st_size)
                throw std::ios_base::failure{"Error: could not load CSV data.\n"};
            close(fd);
            char *sbeg = _sdata; // pointer to start of field string
            char *send = _sdata; // pointer to end of field string
            char *fend = _sdata + buff.st_size - 1; // pointer to last char in buffer
            if (*fend != '\n')
                throw csv_format_error{"Error: CSV file does not end in newline.\n"};
            uint64_t flen; // field length
            if (header) {
                do {
                    std::cout << "top" << std::endl;
                    if (*send == ',' || *send == '\n') {
                        if ((flen = (uint64_t) (send - sbeg)) > this->lf)
                            this->lf = flen;
                        this->_cols.emplace_back();
                        this->_cols.back()._name = sbeg;
                        sbeg = send + 1;
                        ++this->nf;
                        if (*send == '\n') {
                            *send++ = 0;
                            break;
                        }
                        *send = 0;
                    }
                    std::cout << "bottom" << std::endl;
                } while (++send != fend);
                ++this->nl;
            } else {
                do {
                    if (*send == ',' || *send == '\n') {
                        if ((flen = (uint64_t) (send - sbeg)) > this->lf)
                            this->lf = flen;
                        this->_cols.emplace_back(this->nf);
                        this->_cols.back()._s->emplace_back(sbeg);
                        sbeg = send + 1;
                        ++this->nf;
                        if (*send == '\n') {
                            *send++ = 0;
                            break;
                        }
                        *send = 0;
                    }
                    std::cout << "bottom" << std::endl;
                } while (++send != fend);
                ++this->nl; // keep this like this or add an extra line for inserted indices header?
            }
        }
    public:
        explicit csv(const char *path, bool header = true) {
            this->load_csv(path, header);
            for (const auto &c : _cols) {
                std::cout << c._name << std::endl;
                for (const auto &s : *c._s)
                    std::cout << s << std::endl;
            }
        }
    };
}

#endif
