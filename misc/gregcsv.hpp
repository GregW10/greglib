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
    template <std::floating_point>
    class csv;
    template <std::floating_point T>
    class csv_col {
        using VC = std::vector<char*>;
        using VT = std::vector<T>;
        char *_name{};
        VC *_s{};
        VT *_v{};
        bool alloced = false;
        uint64_t _lf = 0; // longest field
        VC _to_free{};
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
        csv_col(const csv_col<T> &other) : _name{new char[gtd::strlen_c(other._name) + 1]},
                                              _s{other._s ? new VC{*other._s} : nullptr},
                                              _v{other._v ? new VT{*other._v} : nullptr},
                                              alloced{true} {
            gtd::strcpy_c(this->_name, other._name);
        }
        csv_col(csv_col<T> &&other) noexcept : _name{other._name},
                                                  _s{other._s},
                                                  _v{other._v},
                                             alloced{other.alloced} {
            other._name = nullptr;
            other._s = nullptr;
            other._v = nullptr;
            other.alloced = false;
        }
        const char *name() const noexcept {
            return this->_name;
        }
        uint64_t size() const noexcept {
            return this->_s->size();
        }
        bool has_nums() const noexcept {
            return this->_v;
        }
        void push(const char *str, bool dup_str = true) {
            /* Appends string to string array.
             * If float array set, attempts to convert to float. If fails, adds NaN. */
            if (dup_str) {
                char *nstr = new char[gtd::strlen_c(str) + 1];
                gtd::strcpy_c(nstr, str);
                this->_to_free.emplace_back(nstr);
                this->_s->emplace_back(nstr);
            } else
                this->_s->emplace_back(str);
            if (this->_v) {
                T val;
                if (gtd::to_float(str, &val))
                    this->_v->emplace_back(val);
                else
                    this->_v->emplace_back(std::numeric_limits<T>::quiet_NaN());
            }
        }
        void push(const T &_num) {
            if (this->_v)
                this->_v->emplace_back(_num);
            if (_num != _num) {
                this->_s->emplace_back("nan");
                return;
            }
            char *buffer = new char[128];
            this->_to_free.emplace_back(buffer);
            int chars_written;
            if constexpr (std::same_as<T, float>)
                chars_written = snprintf(buffer, 128, "%f", _num);
            if constexpr (std::same_as<T, double>)
                chars_written = snprintf(buffer, 128, "%lf", _num);
            if constexpr (std::same_as<T, long double>)
                chars_written = snprintf(buffer, 128, "%Lf", _num);
            if (chars_written > this->_lf)
                this->_lf = chars_written;
            this->_s->emplace_back(buffer);
        }
        char *operator()(uint64_t sindex) { // no bounds-checking
            return this->_s->operator[](sindex);
        }
        T &operator[](uint64_t nindex) { // no bounds-checking
            return this->_v->operator[](nindex);
        }
        csv_col<T> &operator=(const csv_col<T> &other) {
            _name = new char[gtd::strlen_c(other._name) + 1];
            _s = other._s ? new VC{*other._s} : nullptr;
            _v = other._v ? new VT{*other._v} : nullptr;
            alloced = true;
        }
        csv_col<T> &operator=(csv_col<T> &&other) noexcept {
            this->_name = other._name;
            this->_s = other._s;
            this->_v = other._v;
            this->alloced = other.alloced;
            other._name = nullptr;
            other._s = nullptr;
            other._v = nullptr;
            other.alloced = false;
        }
        ~csv_col() {
            if (alloced)
                delete [] _name;
            if (this->_s)
                delete this->_s;
            if (this->_v)
                delete this->_v;
            for (char *& _str : this->_to_free)
                delete [] _str;
        }
        template <std::floating_point>
        friend class csv;
        // template <std::floating_point U>
        friend std::ostream &operator<<(std::ostream&, const csv<T>&);
    };
    template <std::floating_point T = long double>
    class csv {
        std::vector<csv_col<T>> _cols;
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
            csv_col<T> *colbeg = this->_cols.data();
            csv_col<T> *colptr = colbeg;
            do {
                if (*send == ',' || *send == '\n') {
                    if (++fc == this->nf && *send == ',') {
                        char error[60];
                        snprintf(error, 60, "Error: too many fields found on line %" PRIu64 "\n.", this->nl + (options & read_hdr == read_hdr));
                        throw csv_format_error{error};
                    }
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
            if (options & conv_num) {
                std::istringstream iss; // remove the use of a string stream and instead use your float-conversion function
                T val;
                uint64_t nl_m1 = this->nl - 1;
                char **cptr;
                bool at_least_one;
                if (options & conv_all) {
                    if (options & prelimfp) {
                        std::vector<bool> bvec;
                        bvec.reserve(nl_m1);
                        bool bval;
                        for (csv_col<T> &_col : this->_cols) {
                            at_least_one = false;
                            for (const char* const& _p : *(_col._s)) {
                                iss.str(_p);
                                bvec.push_back((bval = (bool) (iss >> val)));
                                at_least_one |= bval;
                            }
                            if (!at_least_one) {
                                bvec.clear();
                                continue;
                            }
                            _col._v = new std::vector<T>{};
                            _col._v->reserve(nl_m1);
                            cptr = _col._s->data();
                            for (const auto &_b : bvec) {
                                if (_b) {
                                    iss.str(*cptr++);
                                    iss >> val;
                                    _col._v->emplace_back(val);
                                } else
                                    _col._v->emplace_back(std::numeric_limits<T>::quiet_NaN());
                            }
                        }
                    } else {
                        for (csv_col<T> &_col : this->_cols) {
                            _col._v = new std::vector<T>{};
                            _col._v->reserve(nl_m1);
                            at_least_one = false;
                            for (const char* const& _p : *(_col._s)) {
                                iss.str(_p);
                                if (iss >> val) {
                                    _col._v->emplace_back(val);
                                    at_least_one = true;
                                }
                                else
                                    _col._v->emplace_back(std::numeric_limits<T>::quiet_NaN());
                            }
                            if (!at_least_one) {
                                delete _col._v;
                                _col._v = nullptr;
                            }
                        }
                    }
                    return;
                }
                if (options & prelimfp) {
                    bool all_good;
                    for (csv_col<T> &_col : this->_cols) {
                        all_good = true;
                        for (const char* const& _p : *(_col._s)) {
                            iss.str(_p);
                            if (!(iss >> val)) {
                                all_good = false;
                                break;
                            }
                        }
                        if (!all_good)
                            continue;
                        _col._v = new std::vector<T>{};
                        _col._v->reserve(nl_m1);
                        for (const char* const& _p : *(_col._s)) {
                            iss.str(_p);
                            iss >> val;
                            _col._v->emplace_back(val);
                        }
                    }
                    return;
                }
                for (csv_col<T> &_col : this->_cols) {
                    _col._v = new std::vector<T>{};
                    for (const char* const& _p : *(_col._s)) {
                        iss.str(_p);
                        if (iss >> val)
                            _col._v->emplace_back(val);
                        else {
                            delete _col._v;
                            _col._v = nullptr;
                            break;
                        }
                    }
                }
            }
        }
    public:
        static constexpr int read_hdr = 1; // if set, treat the first row of the CSV as the header
        static constexpr int conv_num = 2; // if set, attempt to convert values to type T
        static constexpr int conv_all = 4; // if set, do not "give up" on columns whose values are not all numeric (store non-numeric ones as NaN) - requires "conv_num"
        static constexpr int prelimfp = 8; // if set, perform a forward-pass over column first to check if all values can be converted, before allocating memory - requires "conv_num"
        explicit csv(const char *path, int options = read_hdr | conv_num) {
            this->load_csv(path, options);
        }
        csv_col<T> &operator[](const char *field) {
            if (!field)
                throw std::invalid_argument{"Error: field name cannot be nullptr.\n"};
            for (csv_col<T> &_col : this->_cols)
                if (gtd::str_eq(field, _col._name))
                    return _col;
            throw std::invalid_argument{"Error: field passed doesn't exist.\n"};
        }
        // template <std::floating_point U>
        friend std::ostream &operator<<(std::ostream &os, const csv<T> &f) {
            static char buff[256] = {0};
            uint64_t counter;
            if (!buff[0]) {
                counter = 256;
                char *ptr = buff;
                while (counter --> 0)
                    *ptr++ = '-';
            }
            const csv_col<T> *colbeg = f._cols.data();
            const csv_col<T> *colptr = colbeg;
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
