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
    class csv;
    template <gml::Numeric IntType, gml::Numeric FloatType>
    class csv_col {
        using I = IntType;
        using F = FloatType;
        using VC = std::vector<char*>;
        using VI = std::vector<I>;
        using VF = std::vector<F>;
        char *_name{};
        VC *_s{};
        VI *_i{};
        VF *_f{};
        bool alloced = false;
        uint64_t _lf = 0; // longest field
    public:
        csv_col() : _s{new VC{}} {}
        explicit csv_col(uint64_t index, uint64_t *out_digs = nullptr) : _s{new VC{}} {
            uint64_t num_digs = ((uint64_t) log10l((long double) index)) + 1;
            _name = new char[num_digs + 1];
            alloced = true;
            snprintf(_name, num_digs + 1, "%" PRIu64, index);
            if (out_digs)
                *out_digs = num_digs;
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
            this->_s = new VC{};
        }
        csv_col(const csv_col<I, F> &other) : _name{new char[gtd::strlen_c(other._name) + 1]},
                                              _s{other._s ? new VC{*other._s} : nullptr},
                                              _i{other._i ? new VI{*other._i} : nullptr},
                                              _f{other._f ? new VF{*other._f} : nullptr},
                                              alloced{true} {
            gtd::strcpy_c(this->_name, other._name);
        }
        csv_col(csv_col<I, F> &&other) : _name{other._name},
                                         _s{other._s},
                                         _i{other._i},
                                         _f{other._f},
                                         alloced{other.alloced} {
            other._name = nullptr;
            other._s = nullptr;
            other._i = nullptr;
            other._f = nullptr;
            other.alloced = false;
        }
        csv_col<I, F> &operator=(const csv_col<I, F> &other) {
            _name = new char[gtd::strlen_c(other._name) + 1];
            _s = other._s ? new VC{*other._s} : nullptr;
            _i = other._i ? new VI{*other._i} : nullptr;
            _f = other._f ? new VF{*other._f} : nullptr;
            alloced = true;
        }
        csv_col<I, F> &operator=(csv_col<I, F> &&other) {
            this->_name = other._name;
            this->_s = other._s;
            this->_i = other._i;
            this->_f = other._f;
            this->alloced = other.alloced;
            other._name = nullptr;
            other._s = nullptr;
            other._i = nullptr;
            other._f = nullptr;
            other.alloced = false;
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
        template <gml::Numeric II, gml::Numeric FF>
        friend std::ostream &operator<<(std::ostream &os, const csv<II, FF> &f);
    };
    template <gml::Numeric IntType = uint64_t, gml::Numeric FloatType = long double>
    class csv {
        using I = IntType;
        using F = FloatType;
        std::vector<csv_col<I, F>> _cols;
        gtd::mmapper _smap{};
        char *_sdata{};
        uint64_t lf = 0; // longest field
        uint64_t nf = 0; // number of fields
        uint64_t nl = 0; // number of lines (including header)
        void load_csv(const char *path, int options) {
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
            char *fend = _sdata + buff.st_size; // pointer to past last char in buffer
            if (*(fend - 1) != '\n')
                throw csv_format_error{"Error: CSV file does not end in newline.\n"};
            uint64_t flen; // field length
            uint64_t tlen; // temp var.
            if (options & read_hdr) {
                do {
                    if (*send == ',' || *send == '\n') {
                        //if ((flen = (uint64_t) (send - sbeg)) > this->lf)
                            //this->lf = flen;
                        this->_cols.emplace_back();
                        this->_cols.back()._name = sbeg;
                        this->_cols.back()._lf = (uint64_t) (send - sbeg);
                        sbeg = send + 1;
                        ++this->nf;
                        if (*send == '\n') {
                            *send++ = 0;
                            break;
                        }
                        *send = 0;
                    }
                    ++send;
                } while (1);
                ++this->nl;
            } else {
                do {
                    if (*send == ',' || *send == '\n') {
                        //if ((flen = (uint64_t) (send - sbeg)) > this->lf)
                            //this->lf = flen;
                        this->_cols.emplace_back(this->nf, &flen);
                        this->_cols.back()._s->emplace_back(sbeg);
                        tlen = (uint64_t) (send - sbeg);
                        this->_cols.back()._lf = flen > tlen ? flen : tlen;
                        sbeg = send + 1;
                        ++this->nf;
                        if (*send == '\n') {
                            *send++ = 0;
                            break;
                        }
                        *send = 0;
                    }
                    ++send;
                } while (1);
                this->nl = 2; // keep this like this or add an extra line for inserted indices header?
            }
            uint64_t fc = 0; // field counter
            csv_col<I, F> *colbeg = this->_cols.data();
            csv_col<I, F> *colptr = colbeg;
            do {
                if (*send == ',' || *send == '\n') {
                    if (++fc == this->nf && *send == ',') {
                        char error[60];
                        snprintf(error, 60, "Error: too many fields found on line %" PRIu64 "\n.", this->nl + (options & read_hdr == read_hdr));
                        throw csv_format_error{error};
                    }
                    //if ((flen = (uint64_t) (send - sbeg)) > this->lf)
                        //this->lf = flen;
                    if ((flen = (uint64_t) (send - sbeg)) > colptr->_lf)
                        colptr->_lf = flen;
                    colptr->_s->emplace_back(sbeg);
                    sbeg = send + 1;
                    if (*send == '\n') {
                        ++this->nl;
                        if (fc < this->nf) {
                            char error[60];
                            snprintf(error, 60, "Error: too few fields found on line %" PRIu64 "\n.", this->nl + (options & read_hdr == read_hdr));
                            throw csv_format_error{error};
                        }
                        *send = 0;
                        fc = 0;
                        colptr = colbeg;
                        continue;
                    }
                    *send = 0;
                    ++colptr;
                }
            } while (++send != fend);
        }
    public:
        static constexpr int read_hdr = 1; // if set, treat the first row of the CSV as the header
        static constexpr int conv_int = 2; // if set, attempt to convert values to integers
        static constexpr int conv_flt = 4; // if set, attempt to convert values to floats
        static constexpr int conv_all = 8; // if set, do not "give up" on columns whose values are not all numeric (store non-numeric ones as NaN) - requires at least one of "conv_int" or "conv_flt"
        explicit csv(const char *path, int options = read_hdr) {
            this->load_csv(path, options);
        }
        template <gml::Numeric II, gml::Numeric FF>
        friend std::ostream &operator<<(std::ostream &os, const csv<II, FF> &f) {
            static char buff[256] = {0};
            uint64_t counter;
            if (!buff[0]) {
                counter = 256;
                char *ptr = buff;
                while (counter --> 0)
                    *ptr++ = '-';
            }
            const csv_col<I, F> *colbeg = f._cols.data();
            const csv_col<I, F> *colptr = colbeg;
            uint64_t totw = f.nf*3 + 1; // total summed width of all maximum field lengths of columns
            counter = f.nf;
            while (counter --> 0)
                totw += colptr++->_lf;
#define HYPHENS \
            counter = totw; \
            while (1) { \
                if (counter <= 256) { \
                    os.write(buff, counter); \
                    break; \
                } \
                os.write(buff, 256); \
                counter -= 256; \
            } \
            os << '\n';
            HYPHENS
            counter = f.nf;
            colptr = colbeg;
            while (counter --> 0) {
                os << "| " << std::left << std::setw(colptr->_lf) << colptr++->_name << ' ';
            }
            os << "|\n";
            HYPHENS
            counter = 0; // starting at 1 since header has been printed
            uint64_t ic; // inner counter
            uint64_t tt = f.nl - 1;
            while (counter < tt) {
                ic = f.nf;
                colptr = colbeg;
                while (ic --> 0) {
                    os << "| " << std::left << std::setw(colptr->_lf) << (colptr++->_s->operator[](counter)) << ' ';
                }
                os << "|\n";
                ++counter;
            }
            HYPHENS
#undef HYPHENS
            return os;
        }
    };
}
#endif
