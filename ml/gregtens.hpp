#ifndef GREGTENS_HPP
#define GREGTENS_HPP

/* Some acronyms:
 *     LSD: Least significant dimension
 *     MSD: Most significant dimension
 * */

#include <iostream>
#include <fstream>
#include <vector>
#include <concepts>
#include <initializer_list>
#include <sstream>
#include <limits>
#include <sys/stat.h>
#include "greggen.hpp"

namespace gml {
    namespace exceptions {
        class tensor_error : public std::logic_error {
        public:
            tensor_error() : std::logic_error{"Error: an invalid operation has been performed on a tensor.\n"} {}
            explicit tensor_error(const char *msg) : std::logic_error{msg} {}
        };
        class invalid_indexing_error : public tensor_error {
        public:
            invalid_indexing_error() : tensor_error{"Error: an invalid attempt has been made to index a tensor.\n"} {}
            explicit invalid_indexing_error(const char *msg) : tensor_error{msg} {}
        };
        class indexing_dimension_error : public invalid_indexing_error {
            char *message{};
            static constexpr size_t max_buff_size = 126 + std::numeric_limits<unsigned int>::digits10 +
                    std::numeric_limits<size_t>::digits10; // max. possible length of string
        public:
            indexing_dimension_error() : invalid_indexing_error{"Error: an attempt has been made to index a tensor "
                                                                "with an invalid number of dimensions.\n"} {}
            explicit indexing_dimension_error(const char *msg) : invalid_indexing_error{msg} {}
            indexing_dimension_error(unsigned int rank, size_t num_indices) {
                message = new char[max_buff_size]{};
                snprintf(message, max_buff_size, "Error: it is impossible to index a tensor of rank (dimension) %u "
                                                 "using %zu indices (the rank must equal the number of indices).\n",
                                                 rank, num_indices);
            }
            const char *what() const noexcept override {
                return message ? message : invalid_indexing_error::what();
            }
            ~indexing_dimension_error() override {
                delete [] message;
            }
        };
        class out_of_bounds_error : public invalid_indexing_error {
        public:
            out_of_bounds_error() : invalid_indexing_error{"Error: index provided is out of bounds for "
                                                           "the given tensor.\n"} {}
            explicit out_of_bounds_error(const char *msg) : invalid_indexing_error{msg} {}
        };
        class invalid_tsr_format : public tensor_error {
        public:
            invalid_tsr_format() : tensor_error{"Error: the file format of the .tsr file is invalid.\n"} {}
            explicit invalid_tsr_format(const char *msg) : tensor_error{msg} {}
        };
        class invalid_mtsr_format : public invalid_tsr_format {
        public:
            invalid_mtsr_format() : invalid_tsr_format{"Error: the file format of the .mtsr file is invalid.\n"} {}
            explicit invalid_mtsr_format(const char *msg) : invalid_tsr_format{msg} {}
        };
    }
    template <Numeric T>
    class tensor;
    class tensor_shape {
    private:
        unsigned int _r{}; // rank of the tensor
        ull_t *_s{}; // its shape
        ull_t *sizes = _r ? new ull_t[_r - 1] : nullptr; // holds the num. of elements in previous layers
        virtual ull_t offset_at(const std::initializer_list<ull_t> &indices) const {
            unsigned int num_indices = (unsigned int) indices.size();
            if (num_indices != _r)
                throw exceptions::indexing_dimension_error(_r, num_indices);
            if (!_r)
                return 0; // in case of a 0D tensor (a scalar - pretty dumb to use a tensor for this)
            const ull_t *ptr = indices.begin();
            const ull_t *dim = this->_s;
            ull_t offset = *(ptr + _r - 1);
            ull_t *size = sizes;
            unsigned int counter = 1;
            while (counter++ < _r) {
                if (*ptr >= *dim++)
                    throw exceptions::out_of_bounds_error{"Index out of bounds.\n"};
                offset += (*ptr++)*(*size++);
            }
            return offset;
        }
        template <typename ...types> requires (std::is_integral<types>::value && ...)
        ull_t offset_at(types... indices) const {
            if (sizeof...(indices) != this->_r)
                throw exceptions::indexing_dimension_error(_r, sizeof...(indices));
            if (((indices < 0) || ...))
                throw std::invalid_argument{"Error: an index cannot be negative.\n"};
            if (!this->_r)
                return 0; // in case of a 0D tensor (a scalar - pretty dumb to use a tensor for this)
            ull_t idx[] = {static_cast<ull_t>(indices)...};
            // the algorithm I've commented out below allows me to remove the pre-calculated sizes (computed in
            // calc_sizes), but is more expensive if the dimensions of a tensor are not changed during many
            // repeated `offset_at` calls
            /* const ull_t *ptr = ((ull_t *) idx) + this->_r - 1; // point to last index (of LSD)
            const ull_t *dim = this->_s + this->_r - 1; // points to LSD
            unsigned int counter = this->_r;
            ull_t prod = 1;
            ull_t offset = 0;
            while (counter-- > 0) {
                if (*ptr >= *dim)
                    throw exceptions::out_of_bounds_error{"Index out of bounds.\n"};
                offset += prod*(*ptr--);
                prod *= *dim--;
            } */
            const ull_t *ptr = idx;
            const ull_t *dim = this->_s;
            ull_t offset = *(ptr + _r - 1);
            ull_t *size = sizes;
            unsigned int counter = 1;
            while (counter++ < _r) {
                if (*ptr >= *dim++)
                    throw exceptions::out_of_bounds_error{"Index out of bounds.\n"};
                offset += (*ptr++)*(*size++);
            }
            return offset;
        }
        void calc_sizes() noexcept {
            if (!_r || _r == 1)
                return;
            unsigned int counter = 1; // first dimension NOT required for index calculations
            ull_t *dim = _s + _r - 1; // points to LSD (least significant dimension)
            ull_t prod = 1;
            ull_t *size = sizes + _r - 2; // if _r == 1, size will never be de-referenced (so no SIGSEGV)
            while (counter++ < _r) {
                *size = prod*(*dim--);
                prod = *size--;
            }
        }
        explicit tensor_shape(bool alloc = true) {
            std::cout << "\n-----------\nDefault tensor_shape ctor being called.\n-------------" << std::endl;
            if (alloc)
                this->_s = new ull_t[0];
        }
    public:
        // shape() : _r{0}, _s{new ull_t[0]} {}
        explicit tensor_shape(unsigned int rank) : _r{rank}, _s{new ull_t[rank]} {calc_sizes();} // ERROR - NO SHAPE
        tensor_shape(unsigned int rank, ull_t *shape) : _r{rank}, _s{new ull_t[rank]} {
            gen::memcopy(_s, shape, sizeof(ull_t), rank);
            calc_sizes();
        }
        tensor_shape(const std::initializer_list<ull_t> &_shape) : _r{(unsigned int) _shape.size()},
                                                            _s{new ull_t[_r]} {
            gen::memcopy(_s, _shape.begin(), sizeof(ull_t), _r);
            calc_sizes();
        }
        tensor_shape(const tensor_shape &other) : _r{other._r}, _s{new ull_t[_r]} {
            gen::memcopy(_s, other._s, sizeof(ull_t), _r);
            if (_r)
                gen::memcopy(sizes, other.sizes, sizeof(ull_t), _r - 1);
        }
        tensor_shape(tensor_shape &&other) noexcept : _r{other._r}, _s{other._s}, sizes{other.sizes} {
            other._s = nullptr;
            other.sizes = nullptr;
        }
        tensor_shape &make_zero() {
            this->~tensor_shape();
            this->_r = 0;
            this->_s = new ull_t[0];
            this->sizes = nullptr;
            return *this;
        }
        tensor_shape &set(const std::initializer_list<ull_t> &new_shape) {
            this->_r = new_shape.size();
            this->~tensor_shape();
            this->_s = new ull_t[0];
            this->sizes = this->_r ? new ull_t[this->_r - 1] : nullptr;
            return *this;
        }
        ull_t volume() const noexcept {
            unsigned int counter = 0;
            ull_t vol = 1;
            ull_t *ptr = this->_s;
            while (counter++ < this->_r)
                vol *= *ptr++;
            return vol;
        }
        unsigned int rank() const noexcept {
            return _r;
        }
        tensor_shape copy() const {
            return {*this};
        }
        const ull_t *begin() const noexcept {
            return this->_s;
        }
        const ull_t *end() const noexcept {
            return this->_s + this->_r;
        }
        tensor_shape &operator=(const tensor_shape &other) {
            if (&other == this)
                return *this;
            if (other._r == this->_r) {
                gen::memcopy(_s, other._s, sizeof(ull_t), _r);
                if (_r)
                    gen::memcopy(sizes, other.sizes, sizeof(ull_t), _r - 1);
                return *this;
            }
            this->_r = other._r;
            this->~tensor_shape();
            this->_s = new ull_t[this->_r]{};
            if (this->_r) {
                this->sizes = new ull_t[this->_r - 1];
                gen::memcopy(this->_s, other._s, sizeof(ull_t), this->_r);
                gen::memcopy(this->sizes, other.sizes, sizeof(ull_t), this->_r - 1);
            }
            else
                this->sizes = nullptr;
            return *this;
        }
        tensor_shape &operator=(tensor_shape &&other) noexcept {
            if (&other == this)
                return *this;
            this->_r = other._r;
            this->~tensor_shape();
            this->_s = other._s;
            other._s = nullptr;
            this->sizes = other.sizes;
            other.sizes = nullptr;
            return *this;
        }
        ~tensor_shape() {
            delete [] _s;
            delete [] sizes;
        }
        friend std::ostream &operator<<(std::ostream&, const tensor_shape&);
        template <Numeric T>
        friend class tensor;
        template <Numeric T>
        friend class matrix;
        template <Numeric T>
        friend class vector;
        friend bool operator==(const tensor_shape&, const tensor_shape&);
        friend bool operator!=(const tensor_shape&, const tensor_shape&);
        template <Numeric U>
        friend std::ostream &operator<<(std::ostream&, const tensor<U>&);
    };
    std::ostream &operator<<(std::ostream &out, const tensor_shape &_shape) {
        if (!_shape._r)
            return out << "()";
        unsigned int counter = 0;
        ull_t *ptr = _shape._s;
        out << '(';
        goto start;
        while (counter < _shape._r) {
            out << ", ";
            start:
            out << *ptr++;
            ++counter;
        }
        return out << ')';
    }
    bool operator==(const tensor_shape &s1, const tensor_shape &s2) {
        if (s1._r != s2._r)
            return false;
        unsigned int counter = s1._r;
        ull_t *p1 = s1._s;
        ull_t *p2 = s2._s;
        while (counter-- > 0)
            if (*p1++ != *p2++)
                return false;
        return true;
    }
    bool operator!=(const tensor_shape &s1, const tensor_shape &s2) {
        return !operator==(s1, s2);
    }
    template <Numeric T>
    class tensor {
        /* A class used to represent a tensor object - a multidimensional array of numbers with certain mathematical
         * operations defined. */
        /* Notes:
         *     - There is a distinction between a rank/order 0 tensor and an empty tensor. A tensor of order 0 is,
         *       essentially, a scalar, and, thus, has only 1 element. Its shape is empty, (). An empty tensor is a
         *       special case in which no data is held at all. Its shape is also empty, but it is differentiated from a
         *       tensor of order 0 by the `data` pointer being `nullptr`.
         *     - The `begin` and `end` operations are still well-defined for empty `tensor` objects.
         * */
    protected:
        tensor_shape _shape{};
        ull_t vol{};
        T *data{}; // make sure initialisation order isn't causing any bugs
        static void print_array(std::ostream &os, std::streamsize width, const T *_data, unsigned int rank,
                                const ull_t *curr_dim, const ull_t *_sizes, unsigned int max_rank,
                                unsigned int spaces = 0) {
            if (!curr_dim || !*curr_dim) {
                os << "[]";
                return;
            }
            ull_t num_elems = *curr_dim;
            unsigned int counter = spaces;
            while (counter --> 0)
                os << ' ';
            os << '[';
            if (rank == 1) {
                goto start;
                while (--num_elems > 0) {
                    os << " ";
                    start:
                    os.width(width); // unfortunately, this is necessary, given that the width is reset
                    os << +(*_data++);
                }
            } else {
                unsigned int next_rank = rank - 1;
                const ull_t *next_dim = curr_dim + 1;
                const ull_t *next_size = next_rank == 1 ? nullptr : _sizes + 1; // for completeness
                print_array(os, width, _data, next_rank, next_dim, next_size, max_rank, 0);
                _data += *_sizes;
                --num_elems;
                while (num_elems-- > 0) {
                    os << '\n';
                    print_array(os, width, _data, next_rank, next_dim, next_size, max_rank, 1 + max_rank - rank);
                    _data += *_sizes;
                }
            }
            os << ']';
        }
        static std::streamsize max_width(const T *_data, ull_t size) {
            std::ostringstream oss;
            std::streamsize max_w = 0;
            std::streamsize width;
            while (size --> 0) {
                oss << *_data++;
                if ((width = oss.tellp()) > max_w)
                    max_w = width;
                oss.seekp(0);
            }
            return max_w;
        }
        static std::streamsize max_width(const T *_data, ull_t size, std::streamsize prec = precision)
        requires (std::floating_point<T>) {
            std::ostringstream oss;
            oss.precision(prec);
            std::streamsize max_w = 0;
            std::streamsize width;
            while (size --> 0) {
                oss << *_data++;
                if ((width = oss.tellp()) > max_w)
                    max_w = width;
                oss.seekp(0);
            }
            return max_w;
        }
        void load_tsr(std::ifstream &in, size_t st_size) {
            in.read((char *) &this->_shape._r, sizeof(uint32_t));
            if (st_size == (4 + sizeof(uint32_t))) { // case for empty tensor
                if (this->_shape._r)
                    throw exceptions::invalid_tsr_format{"Invalid .tsr format: file too small for non-empty tensor.\n"};
                this->_shape._s = new ull_t[0];
                return;
            }
            if (!this->_shape._r)
                throw exceptions::invalid_tsr_format{"Invalid .tsr format: file too large for empty tensor.\n"};
            uint64_t fsize;
            size_t arr_shape_size = this->_shape._r*sizeof(uint64_t);
            if (st_size < (fsize = 4 + sizeof(uint32_t) + arr_shape_size + sizeof(uint64_t)))
                throw exceptions::invalid_tsr_format{"Invalid .tsr format: file too small for rank read.\n"};
            this->_shape._s = new uint64_t[this->_shape._r];
            this->_shape.sizes = new uint64_t[this->_shape._r - 1];
            in.read((char *) this->_shape._s, arr_shape_size);
            uint64_t T_size;
            in.read((char *) &T_size, sizeof(uint64_t));
            if (T_size != sizeof(T))
                throw exceptions::invalid_tsr_format{"Invalid .tsr format: size of .tsr tensor data type does not "
                                                     "match the size of the templated type T of this instance.\n"};
            this->_shape.calc_sizes();
            this->vol = this->_shape.volume();
            size_t arr_size = sizeof(T)*this->vol;
            fsize += arr_size;
            if (st_size != fsize)
                throw exceptions::invalid_tsr_format{"Invalid .tsr format: invalid file size.\n"};
            this->data = new T[this->vol];
            in.read((char *) this->data, arr_size);
            in.close();
        }
        void process_tsr(const char *path) {
            if (!path || !*path)
                throw std::invalid_argument{"Error: null or empty path provided to .tsr file.\n"};
            struct stat buff{};
            if (stat(path, &buff) == -1)
                throw std::ios_base::failure{"Error: could not access file.\n"};
            if (buff.st_size < (4 + sizeof(uint32_t)))
                throw exceptions::invalid_tsr_format{"Invalid .tsr format: file size too small.\n"};
            std::ifstream in{path, std::ios_base::in | std::ios_base::binary};
            if (!in.good())
                throw std::ios_base::failure{"Error: could not open file.\n"};
            char header[4];
            in.read(header, 4*sizeof(char));
            if (header[0] || header[1] != 't' || header[2] != 's' || header[3] != 'r')
                throw exceptions::invalid_tsr_format{"Invalid .tsr format: invalid header.\n"};
            load_tsr(in, buff.st_size);
        }
        explicit tensor(bool alloc) : data{}, vol{0}, _shape(alloc) {}
    public:
        static inline std::streamsize precision = 6; // default precision of output streams, as per the C++ standard
        tensor() : data{}, vol{0} {} // constructs an empty tensor
        explicit tensor(const tensor_shape &_s):_shape{_s},vol{_s.volume()},data{new T[this->vol]}/* init. to zeros */{}
        explicit tensor(tensor_shape &&_s) :
        _shape{std::move(_s)}, vol{this->_shape.volume()}, data{new T[this->vol]} /* init. to zeros */ {}
        tensor(const std::vector<T> &_data, const tensor_shape &_s) : _shape{_s}, vol{_s.volume()} {
            if (this->vol != _data.size())
                throw std::invalid_argument{"Error: the number of elements in the data and requested shape of tensor "
                                            "do not match.\n"};
            this->data = new T[vol];
            gen::memcopy(this->data, _data.data(), sizeof(T), this->vol);
        }
        tensor(const T *_data, const tensor_shape &_s) : _shape{_s}, vol{_s.volume()} {
            this->data = new T[vol];
            gen::memcopy(this->data, _data, sizeof(T), this->vol);
        }
        tensor(const tensor<T> &other) : _shape{other.s}, vol{other.vol} {
            this->data = new T[vol];
            gen::memcopy(this->data, other.data, sizeof(T), other.vol);
        }
        tensor(tensor<T> &&other) noexcept : data{other.data}, _shape{std::move(other.s)}, vol{other.vol} {
            other.data = nullptr;
            other.s.make_zero();
            other.vol = 0;
        }
        explicit tensor(const char *tsr_path) : _shape{false} {
            process_tsr(tsr_path);
        }
        template <typename ...types> requires (std::is_integral<types>::value && ...)
        T &at(types... indices) {
            return *(this->data + this->_shape.offset_at(indices...));
        }
        template <typename ...types> requires (std::is_integral<types>::value && ...)
        const T &at(types... indices) const {
            return *(this->data + this->s.offset_at(indices...));
        }
        bool empty() const noexcept {
            return !this->data;
        }
        ull_t volume() const noexcept {
            if (this->data)
                return this->s.volume();
            return 0;
        }
        tensor<T> copy() const {
            return {*this};
        }
        tensor &reshape(const tensor_shape &new_shape) {
            this->_shape = new_shape;
            return *this;
        }
        const tensor_shape &shape() const noexcept {
            return this->_shape;
        }
        T *begin() noexcept {
            /* If `this->data` is `nullptr` (empty tensor), this operation is still defined. */
            return this->data;
        }
        T *end() noexcept {
            /* If `this->data` is `nullptr` (empty tensor), this operation is well-defined, as adding a zero offset to
             * a `nullptr` results in another `nullptr`. Thus, when performing iteration, the pointer returned by
             * `begin` will compare as equivalent to that returned by `end`, and terminate the loop. */
            return this->data + vol;
        }
        std::ofstream::pos_type to_tsr(const char *path) const { // writes the data of the tensor to a .tsr file
            if (!path || !*path)
                throw std::invalid_argument{"Error: nullptr or empty string passed as path.\n"};
            std::ofstream out{path, std::ios_base::out | std::ios_base::trunc};
            if (!out.good())
                throw std::ios_base::failure{"Error: could not open file.\n"};
            out.put(0); out.put('t'); out.put('s'); out.put('r');
            out.write((char *) &this->_shape._r, sizeof(uint32_t));
            if (!this->data) { // this means the tensor is empty, so nothing more is written (just the rank of 0)
                out.close();
                return 4 + sizeof(uint32_t); // == 8
            }
            out.write((char *) this->_shape._s, this->_shape._r*sizeof(ull_t));
            ull_t size = sizeof(T);
            out.write((char *) &size, sizeof(ull_t));
            out.write((char *) this->data, sizeof(T)*this->vol);
            std::ofstream::pos_type bytes = out.tellp();
            out.close();
            return bytes;
        }
        ~tensor() {
            delete [] data;
        }
        template <Numeric U>
        friend std::ostream &operator<<(std::ostream&, const tensor<U>&);
        template <Numeric U, Numeric V>
        friend bool operator==(const tensor<U>&, const tensor<V>&);
        template <Numeric U, Numeric V>
        friend bool operator!=(const tensor<U>&, const tensor<V>&);
    };
    template <Numeric T>
    std::ostream &operator<<(std::ostream &os, const tensor<T> &tens) {
        std::streamsize prev_width = os.width();
        std::streamsize max_w;
        if constexpr (std::floating_point<T>) {
            std::streamsize prev_prec = os.precision(tensor<T>::precision);
            max_w = tensor<T>::max_width(tens.data, tens.vol, tensor<T>::precision);
            tensor<T>::print_array(os, max_w, tens.data, tens._shape._r, tens._shape._s, tens._shape.sizes,
                                   tens._shape._r);
            os.precision(prev_prec);
        } else {
            max_w = tensor<T>::max_width(tens.data, tens.vol);//, tensor<T>::precision);
            tensor<T>::print_array(os, max_w, tens.data, tens._shape._r, tens._shape._s, tens._shape.sizes,
                                   tens._shape._r);
        }
        os.width(prev_width);
        return os;
    }
    template <Numeric U, Numeric V>
    bool operator==(const tensor<U> &t1, const tensor<V> &t2) {
        if (t1._shape != t2._shape)
            return false;
        U *p1 = t1.data;
        V *p2 = t2.data;
        ull_t counter = t1.vol;
        while (counter --> 0)
            if (*p1++ != *p2++)
                return false;
        return true;
    }
    template <Numeric U, Numeric V>
    bool operator!=(const tensor<U> &t1, const tensor<V> &t2) {
        return !operator==(t1, t2);
    }
    typedef tensor<long double> tens_ld;
    typedef tensor<double> tens_d;
    typedef tensor<float> tens_f;
    typedef tensor<unsigned long long> tens_ull;
    typedef tensor<long long> tens_ll;
    typedef tensor<unsigned long> tens_ul;
    typedef tensor<long> tens_l;
    typedef tensor<unsigned int> tens_ui;
    typedef tensor<int> tens_i;
}
#endif
