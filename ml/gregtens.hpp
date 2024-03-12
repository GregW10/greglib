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
        class dimension_mismatch_error : public tensor_error {
        public:
            dimension_mismatch_error() : tensor_error{"Error: there is a dimension mismatch between tensors.\n"} {}
            explicit dimension_mismatch_error(const char *msg) : tensor_error{msg} {}
        };
        class empty_tensor_error : public tensor_error {
        public:
            empty_tensor_error() :tensor_error{"Error: the given operation cannot be performed on an empty tensor.\n"}{}
            explicit empty_tensor_error(const char *msg) : tensor_error{msg} {}
        };
    }
    template <Numeric>
    class tensor;
    template <Numeric>
    class matrix;
    template <Numeric>
    class vector_base;
    template <Numeric>
    class vector;
    namespace tens_ops {
        template <Numeric U, Numeric V>
        tensor<mulComT<U, V>> hadamard(const tensor<U>&, const tensor<V>&);
        template <Numeric U, Numeric V>
        matrix<mulComT<U, V>> hadamard(const matrix<U>&, const matrix<V>&);
        template <Numeric U, Numeric V>
        vector<mulComT<U, V>> hadamard(const vector<U>&, const vector<V>&);
    }
    class tensor_shape {
        /* Class whose instances describe the shape of a given tensor. I chose to not make it a nested class inside
         * `tensor<T>`, as it should not and does not need to be templated.
         * Notes:
         *          - `_r` is the rank of the tensor. A 0D tensor is a scalar and has a rank of zero, with a volume of 1
         *            (i.e. it has one element). An empty tensor, on the other hand, is given a rank of -1.
         *          - `_s` is an array containing the dimensions of the tensor. It is always allocated for any non-empty
         *            tensor. A shape describing a 0D tensor simply allocates an array of zero size. An empty tensor, on
         *            the other hand, sets _s to be a nullptr.
         *          - `sizes` is an array of size `_r - 1`, containing the sizes of previous layers. It is used only for
         *            calculating the required offsets when indexing a tensor. If the tensor is 1D, then `sizes` is
         *            still allocated, but with zero size. If the tensor is 0D or empty, then `sizes` is `nullptr`. */
    private:
        int64_t _r{-1}; // rank of the tensor - is signed integral because has to be able to be -1 for empty tensor
        uint64_t *_s{}; // its shape
        uint64_t *sizes = _r > 0 ? new ull_t[_r - 1] : nullptr; // holds the num. of elements in previous layers
        virtual ull_t offset_at(const std::initializer_list<ull_t> &indices) const {
            if (this->_r == -1) // case for an empty tensor
                throw exceptions::empty_tensor_error{"Error: the offset for an empty tensor is undefined as no data is "
                                                     "held.\n"};
            int64_t num_indices = (int64_t) indices.size();
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
            if (this->_r == -1) // case for an empty tensor
                throw exceptions::empty_tensor_error{"Error: the offset for an empty tensor is undefined as no data is "
                                                     "held.\n"};
            if (sizeof...(indices) != this->_r)
                throw exceptions::indexing_dimension_error(_r, sizeof...(indices));
            if constexpr (sizeof...(indices)) {
                if (((indices < 0) || ...))
                    throw std::invalid_argument{"Error: an index cannot be negative.\n"};
            } else  // thus, if `offset_at` is ever called with zero indices, the compiler will generate a version of
                return 0; // this function which always returns zero
            // if (!this->_r)
            //     return 0; // in case of a 0D tensor (a scalar - pretty dumb to use a tensor for this)
            uint64_t idx[] = {static_cast<uint64_t>(indices)...};
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
            const uint64_t *ptr = idx;
            const uint64_t *dim = this->_s;
            uint64_t offset = *(ptr + _r - 1);
            uint64_t *size = sizes;
            unsigned int counter = 1;
            while (counter++ < _r) {
                if (*ptr >= *dim++)
                    throw exceptions::out_of_bounds_error{"Index out of bounds.\n"};
                offset += (*ptr++)*(*size++);
            }
            return offset;
        }
        void calc_sizes() noexcept {
            if (this->_r <= 1) // there are no previous layers if the tensor is 1D, 0D or empty
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
        // explicit tensor_shape(bool alloc) {
        //     if (alloc)
        //         this->_s = new ull_t[0];
        // }
        tensor_shape(int64_t rank, uint64_t *_shape, bool make_copy = true) : _r{rank} {
            /* Private ctor for efficient construction of a tensor_shape object (used internally by other classes). */
            /* This ctor is never called for creating a shape describing an empty tensor. */
            if (make_copy) {
                _s = new uint64_t[rank];
                gen::memcopy(_s, _shape, sizeof(uint64_t), rank);
            } else
                _s = _shape;
            this->calc_sizes();
        }
    public:
        tensor_shape() = default; // only way to construct an empty tensor shape
        explicit tensor_shape(uint32_t rank) : _r{rank}, _s{new uint64_t[rank]} {calc_sizes();}
        tensor_shape(uint32_t rank, uint64_t *_shape) : _r{rank}, _s{new uint64_t[rank]} {
            gen::memcopy(_s, _shape, sizeof(uint64_t), rank);
            calc_sizes();
        }
        tensor_shape(const std::initializer_list<uint64_t> &_shape) : _r{(unsigned int) _shape.size()},
                                                            _s{new uint64_t[_r]} {
            gen::memcopy(_s, _shape.begin(), sizeof(uint64_t), _r);
            calc_sizes();
        }
        tensor_shape(const tensor_shape &other) : _r{other._r} {
            if (this->_r == -1) {
                return;
            }
            this->_s = new uint64_t[_r];
            gen::memcopy(_s, other._s, sizeof(uint64_t), _r);
            if (_r)
                gen::memcopy(sizes, other.sizes, sizeof(uint64_t), _r - 1);
        }
        tensor_shape(tensor_shape &&other) noexcept : _r{other._r}, _s{other._s}, sizes{other.sizes} {
            other._r = -1; // the other shape is left as describing an empty tensor
            other._s = nullptr;
            other.sizes = nullptr;
        }
        tensor_shape &make_empty() { // makes the shape describe an empty tensor (not 0D, actually properly empty)
            this->~tensor_shape();
            this->_r = -1;
            this->_s = nullptr;
            this->sizes = nullptr;
            return *this;
        }
        tensor_shape &set(const std::initializer_list<uint64_t> &new_shape) {
            /* It is impossible to use this method to set the shape as describing an empty tensor. The `make_empty`
             * method should be used for that. */
            this->_r = (int64_t) new_shape.size();
            this->~tensor_shape();
            this->_s = new uint64_t[this->_r];
            this->sizes = this->_r ? new uint64_t[this->_r - 1] : nullptr;
            return *this;
        }
        uint64_t volume() const noexcept {
            if (this->_r == -1)
                return 0;
            int64_t counter = 0;
            uint64_t vol = 1;
            uint64_t *ptr = this->_s;
            while (counter++ < this->_r)
                vol *= *ptr++;
            return vol; // a 0D tensor (a scalar) shape will always return 1, as the tensor contains a single element
        }
        int64_t rank() const noexcept {
            return _r; // returns -1 for an empty tensor shape
        }
        tensor_shape copy() const {
            return {*this};
        }
        const ull_t *begin() const noexcept {
            return this->_s; // is `nullptr` for empty tensor, and non-dereference-able for a 0D tensor, though not null
        }
        const ull_t *end() const noexcept {
            if (this->_r <= 0)
                return this->_s; // thus, if the tensor shape is empty or 0D, `begin == end`
            return this->_s + this->_r;
        }
        tensor_shape &operator=(const tensor_shape &other) {
            if (&other == this)
                return *this;
            if (other._r == this->_r) {
                if (this->_r == -1) { // nothing to do or to copy if both shapes already describe an empty tensor
                    return *this;
                }
                gen::memcopy(_s, other._s, sizeof(uint64_t), _r);
                if (_r)
                    gen::memcopy(sizes, other.sizes, sizeof(uint64_t), _r - 1);
                return *this;
            }
            this->_r = other._r;
            delete [] this->_s;
            delete [] this->sizes;
            if (other._r == -1) { // case for *this becoming a shape describing an empty tensor
                this->_s = nullptr;
                this->sizes = nullptr;
                return *this;
            }
            this->_s = new uint64_t[this->_r]{}; // should I be default initialising if _r == 0 ?
            if (this->_r) {
                this->sizes = new uint64_t[this->_r - 1];
                gen::memcopy(this->_s, other._s, sizeof(uint64_t), this->_r);
                gen::memcopy(this->sizes, other.sizes, sizeof(uint64_t), this->_r - 1);
            }
            else
                this->sizes = nullptr;
            return *this;
        }
        tensor_shape &operator=(tensor_shape &&other) noexcept {
            if (&other == this)
                return *this;
            this->_r = other._r;
            other._r = -1;
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
        template <Numeric>
        friend class tensor;
        template <Numeric>
        friend class matrix;
        template <Numeric>
        friend class vector;
        friend bool operator==(const tensor_shape&, const tensor_shape&);
        friend bool operator!=(const tensor_shape&, const tensor_shape&);
        template <Numeric U>
        friend std::ostream &operator<<(std::ostream&, const tensor<U>&);
        template <Numeric U, Numeric V>
        friend tensor<mulComT<U, V>> tens_ops::hadamard(const tensor<U>&, const tensor<V>&);
        template <Numeric U, Numeric V>
        friend matrix<mulComT<U, V>> tens_ops::hadamard(const matrix<U>&, const matrix<V>&);
        template <Numeric U, Numeric V>
        friend vector<mulComT<U, V>> tens_ops::hadamard(const vector<U>&, const vector<V>&);
        template <Numeric U, Numeric V>
        friend tensor<addComT<U, V>> operator+(const tensor<U>&, const tensor<V>&);
        template <Numeric U, Numeric V>
        friend tensor<subComT<U, V>> operator-(const tensor<U>&, const tensor<V>&);
        template <Numeric U, Numeric V>
        friend tensor<mulComT<U, V>> operator*(const U&, const tensor<V>&);
        template <Numeric U, Numeric V>
        friend tensor<divComT<U, V>> operator/(const tensor<U>&, const V&);
        template <Numeric U, Numeric V>
        friend matrix<addComT<U, V>> operator+(const matrix<U>&, const matrix<V>&);
        template <Numeric U, Numeric V>
        friend matrix<subComT<U, V>> operator-(const matrix<U>&, const matrix<V>&);
        template <Numeric U, Numeric V>
        friend matrix<mulComT<U, V>> operator*(const U&, const matrix<V>&);
        template <Numeric U, Numeric V>
        friend matrix<divComT<U, V>> operator/(const matrix<U>&, const V&);
        template <Numeric U, Numeric V>
        friend vector<addComT<U, V>> operator+(const vector<U>&, const vector<V>&);
        template <Numeric U, Numeric V>
        friend vector<subComT<U, V>> operator-(const vector<U>&, const vector<V>&);
        template <Numeric U, Numeric V>
        friend vector<mulComT<U, V>> operator*(const U&, const vector<V>&);
        template <Numeric U, Numeric V>
        friend vector<divComT<U, V>> operator/(const vector<U>&, const V&);
        template <Numeric U, Numeric V>
        friend matrix<mulComT<U, V>> operator*(const matrix<U>&, const matrix<V>&);
        template <Numeric>
        friend class ffnn;
    };
    std::ostream &operator<<(std::ostream &out, const tensor_shape &_shape) {
        if (_shape._r == -1)
            return out << "(None)";
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
        if (s1._r == -1) // if the shapes describe an empty tensor, it is automatically known that they are equal
            return true;
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
         *       essentially, a scalar, and, thus, has only 1 element. Its shape is `()`. An empty tensor is a
         *       special case in which no data is held at all. Its shape is `(None)`, and it is differentiated from a
         *       tensor of order 0 by the `data` pointer being `nullptr` (and by fields in the `tensor_shape` object).
         *     - The `begin` and `end` operations are still well-defined for empty `tensor` objects.
         *     - There is another breed of empty tensor objects: these are tensor objects that have a positive rank, but
         *       that have at least one dimension that is zero. These are also "empty", as they cannot store any data
         *       (since even one dimension being zero leads to a total tensor volume of zero), but are not treated in
         *       exactly the same way as "truly" empty tensor objects with a rank of -1. Their rank is preserved as
         *       whatever it was stated to be, and the `data` pointer is allocated, just with a size of zero (so it is
         *       not `nullptr`, as in the case of a "truly" empty tensor with a rank of -1).
         * */
    protected:
        tensor_shape _shape{}; // empty shape by default
        uint64_t vol{};
        T *data{}; // make sure initialisation order isn't causing any bugs
        static void print_array(std::ostream &os, std::streamsize width, const T *_data, int64_t rank,
                                const uint64_t *curr_dim, const uint64_t *_sizes, unsigned int max_rank,
                                unsigned int spaces = 0) {
            if (!curr_dim || !*curr_dim) {
                os << "[]";
                return;
            }
            uint64_t num_elems = *curr_dim;
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
                int64_t next_rank = rank - 1;
                const uint64_t *next_dim = curr_dim + 1;
                const uint64_t *next_size = next_rank == 1 ? nullptr : _sizes + 1; // for completeness
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
        static std::streamsize max_width(const T *_data, uint64_t size) {
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
            in.read((char *) &this->_shape._r, sizeof(int64_t));
            if (this->_shape._r == -1) { // case for empty tensor
                if (st_size != (4 + sizeof(int64_t)))
                    throw exceptions::invalid_tsr_format{"Invalid .tsr format: invalid file size for empty tensor.\n"};
                // this->_shape._s = nullptr;
                return;
            }
            if (!this->_shape._r) { // case for 0D (scalar) tensor
                if (st_size != 4 + sizeof(int64_t) + sizeof(uint64_t) + sizeof(T))
                    throw exceptions::invalid_tsr_format{"Invalid .tsr format: invalid file size for a 0D tensor.\n"};
                uint64_t T_size;
                in.read((char *) &T_size, sizeof(uint64_t));
                if (T_size != sizeof(T))
                    throw exceptions::invalid_tsr_format{"Invalid .tsr format: sizeof(T) mismatch.\n"};
                this->vol = 1;
                this->data = new T[1];
                in.read((char *) this->data, sizeof(T));
                in.close();
                this->_shape._s = new uint64_t[0];
                // this->_shape.sizes = nullptr;
                return;
            }
            uint64_t fsize;
            size_t arr_shape_size = this->_shape._r*sizeof(uint64_t);
            if (st_size < (fsize = 4 + sizeof(int64_t) + arr_shape_size + sizeof(uint64_t)))
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
    protected:
        // ctor ONLY to be used within other methods that will manually initialise the shape object:
        // explicit tensor(bool alloc) : data{}, vol{0}, _shape(alloc) {}
        tensor(T *_data, int64_t rank, uint64_t *_s, bool copy_shape = true) :
        _shape{rank, _s, copy_shape}, vol{_shape.volume()}, data{_data} {}
        tensor(T *_data, int64_t rank, uint64_t *_s, uint64_t _volume, bool copy_shape = true) :
        _shape{rank, _s, copy_shape}, vol{_volume}, data{_data} {}
    public:
        static inline std::streamsize precision = 6; // default precision of output streams, as per the C++ standard
        tensor() = default; // constructs an empty tensor
        explicit tensor(const T &val) : _shape(0), vol{1}, data{new T[1]{val}} {} // constructs 0D tensor (one element)
        explicit tensor(const tensor_shape &_s) : _shape{_s}, vol{_s.volume()},
        data{this->_shape._r >= 0 ? new T[this->vol]{} : nullptr} /* init. to zeros if non-empty */ {}
        explicit tensor(tensor_shape &&_s) : _shape{std::move(_s)}, vol{this->_shape.volume()},
        data{this->_shape._r >= 0 ? new T[this->vol]{} : nullptr} {}
        tensor(const std::vector<T> &_data, const tensor_shape &_s) : _shape{_s}, vol{_s.volume()} {
            if (!this->vol) {
                if (!_data.empty())
                    throw exceptions::empty_tensor_error{"Error: an empty tensor cannot be initialised with a non-empty"
                                                         " std::vector.\n"};
                if (this->_shape._r >= 1) // case for an empty positive-rank tensor: `data` is allocated with zero size
                    this->data = new T[this->vol];
                return; // `data` only left as `nullptr` if `_r == -1`
            }
            if (this->vol != _data.size())
                throw std::invalid_argument{"Error: the number of elements in the data and requested shape of tensor "
                                            "do not match.\n"};
            this->data = new T[this->vol];
            gen::memcopy(this->data, _data.data(), sizeof(T), this->vol);
        }
        tensor(const T *_data, const tensor_shape &_s) : _shape{_s}, vol{_s.volume()} {
            if (!this->vol) {
                if (this->_shape._r >= 1) // same reasoning as above ctor
                    this->data = new T[this->vol];
                return;
            }
            if (!_data)
                throw std::invalid_argument{"Error: `nullptr` supplied as pointer to data, but non-empty shape was "
                                            "passed.\n"};
            this->data = new T[this->vol];
            gen::memcopy(this->data, _data, sizeof(T), this->vol);
        }
        tensor(const tensor<T> &other) : _shape{other._shape}, vol{other.vol} {
            if (!this->vol) {
                if (this->_shape._r >= 1) // same reasoning again
                    this->data = new T[this->vol];
                return;
            }
            this->data = new T[this->vol];
            gen::memcopy(this->data, other.data, sizeof(T), other.vol); // watch out, bypasses constructor of T if non-t
        }
        tensor(tensor<T> &&other) noexcept : _shape{std::move(other._shape)}, vol{other.vol}, data{other.data} {
            // other._shape.make_empty(); // redundant because the move already leaves the other shape empty
            other.vol = 0;
            other.data = nullptr; // the other tensor object is left empty (not 0D, but properly empty)
        }
        template <Numeric U> requires (std::is_convertible<U, T>::value)
        tensor(const tensor<U> &other) : _shape{other._shape}, vol{other.vol} {
            if (!this->vol) {
                if (this->_shape._r >= 1)
                    this->data = new T[this->vol];
                return;
            }
            this->data = new T[this->vol];
            gen::copy(this->data, other.data, this->vol);
        }
        explicit tensor(const char *tsr_path) /* : _shape{false} */ {
            process_tsr(tsr_path);
        }
        tensor<T> &assign(const T &value) {
            uint64_t counter = this->vol;
            T *ptr = this->data;
            while (counter --> 0)
                *ptr++ = value;
            return *this;
        }
        template <typename ...types> requires (std::is_integral<types>::value && ...)
        T &at(types... indices) { // no checking performed to not affect efficiency, so could produce SIGSEGV
            return *(this->data + this->_shape.offset_at(indices...));
        }
        template <typename ...types> requires (std::is_integral<types>::value && ...)
        const T &at(types... indices) const {
            return *(this->data + this->_shape.offset_at(indices...));
        }
        virtual T &at(const std::initializer_list<uint64_t> &indices) { // non-templated version, can be overriden
            return *(this->data + this->_shape.offset_at(indices));
        }
        virtual const T &at(const std::initializer_list<uint64_t> &indices) const {
            return *(this->data + this->_shape.offset_at(indices));
        }
        bool empty() const noexcept {
            /* The pending question is whether I should define an N-D tensor (where N > 0) which has one or more
             * dimensions equal to zero as an empty tensor or not. If so, it would be a different breed of empty tensor
             * from the ones in which `this->_shape._r == -1`, as they would have a positive rank, but hold no elements
             * because they have been completely squished down through one or more of their dimensions. */
            /* Decision taken: N-D tensors (where N > 0) with one or more zero-dimensions ARE classified as empty,
             * though they are treated slightly differently to empty tensors with `_r == -1`. Their `data` pointer is
             * allocated a size of zero, and their rank is preserved. Thus, their shape includes one or more zeros. */
            return !this->vol; // or `this->_shape._r == -1` ?
        }
        uint64_t volume() const noexcept {
            // if (this->data)
            //     return this->vol;//this->s.volume();
            // return 0;
            return this->vol;
        }
        int64_t rank() const noexcept {
            return this->_shape._r;
        }
        tensor<T> copy() const {
            return {*this};
        }
        virtual tensor &reshape(const tensor_shape &new_shape) {
            if (this->vol != new_shape.volume())
                throw exceptions::dimension_mismatch_error{"Error: cannot reshape to requested shape as this would "
                                                           "require modifying the number of elements.\n"};
            if (new_shape._r == this->_shape._r) { // case for new shape being same rank, more efficient not to realloc.
                gml::gen::memcopy(this->_shape._s, new_shape._s, sizeof(uint64_t), new_shape._r);
                this->_shape.calc_sizes();
                return *this;
            }
            this->_shape = new_shape;
            return *this;
        }
        const tensor_shape &shape() const noexcept {
            return this->_shape;
        }
        void to_zeros() noexcept { // set all elements to zero
            uint64_t counter = this->vol;
            T *ptr = this->data;
            while (counter --> 0)
                *ptr++ = 0;
        }
        void for_each(void (*func)(T&)) {
            uint64_t counter = this->vol;
            T *ptr = this->data;
            while (counter --> 0)
                func(*ptr++);
        }
        void for_each(T (*func)(T)) {
            uint64_t counter = this->vol;
            T *ptr = this->data;
            while (counter --> 0) {
                *ptr = func(*ptr);
                ++ptr;
            }
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
            out.write((char *) &this->_shape._r, sizeof(int64_t));
            if (this->_shape._r == -1) { // means the tensor is empty, so nothing more is written (just the rank of -1)
                out.close();
                return 4*sizeof(char) + sizeof(int64_t); // == 12
            }
            // if (!this->_shape._r) {
            //     // not needed, as this case is naturally addressed below
            // }
            out.write((char *) this->_shape._s, this->_shape._r*sizeof(uint64_t)); // no write if 0D
            uint64_t size = sizeof(T);
            out.write((char *) &size, sizeof(uint64_t));
            out.write((char *) this->data, sizeof(T)*this->vol); // writes the single element if 0D
            std::ofstream::pos_type bytes = out.tellp();
            out.close();
            return bytes;
        }
        virtual tensor &operator=(const tensor<T> &other) {
            if (&other == this)
                return *this;
            if (!other.vol) { // volume can never be zero for a 0D tensor (will always be 1)
                if (!this->vol) {
                    if (other._shape._r == -1) {
                        if (this->_shape._r == -1) {
                            return *this; // both tensors are truly empty tensors so must already be equivalent
                        } else {
                            delete [] this->data;
                            this->data = nullptr;
                        }
                    } else {
                        if (this->_shape._r == -1) {
                            this->data = new T[0];
                        } /* else {
                            // no need to free T[0] and allocate a new T[0] array
                        } */
                    }
                } else {
                    delete [] this->data;
                    if (other._shape._r >= 1)
                        this->data = new T[0];
                    else
                        this->data = nullptr; // case for a truly tempty tensor with `_r == -1`
                }
                // return *this;
            } else {
                if (this->vol != other.vol) {
                    delete [] this->data;
                    this->data = new T[other.vol];
                }
                gen::memcopy(this->data, other.data, sizeof(T), other.vol);
            }
            // this->data = new T[other.vol];
            // gen::memcopy(this->data, other.data, sizeof(T), other.vol);
            this->_shape = other._shape;
            this->vol = other.vol;
            return *this;
        }
        virtual tensor &operator=(tensor<T> &&other) { // not marked `noexcept` due to derived classes overriding
            if (&other == this)
                return *this;
            this->_shape = std::move(other._shape);
            this->vol = other.vol;
            delete [] this->data;
            this->data = other.data;
            // other._shape.make_empty();
            other.vol = 0;
            other.data = nullptr;
            return *this;
        }
        ~tensor() {
            delete [] data;
        }
        virtual tensor<T> &operator+=(const tensor<T> &other) {
            if (other._shape != this->_shape)
                throw exceptions::dimension_mismatch_error{"Error: cannot sum two tensors with different shapes.\n"};
            T *dptr = this->data;
            T *optr = other.data;
            uint64_t counter = this->vol;
            while (counter --> 0)
                *dptr++ += *optr++;
            return *this;
        }
        virtual tensor<T> &operator-=(const tensor<T> &other) {
            if (other._shape != this->_shape)
                throw exceptions::dimension_mismatch_error{"Error: cannot subtract tensor with different shape.\n"};
            T *dptr = this->data;
            T *optr = other.data;
            uint64_t counter = this->vol;
            while (counter --> 0)
                *dptr++ -= *optr++;
            return *this;
        }
        template <Numeric U, Numeric V>
        friend tensor<addComT<U, V>> operator+(const tensor<U> &t1, const tensor<V> &t2) {
            if (t1._shape != t2._shape)
                throw exceptions::dimension_mismatch_error{"Error: addition between tensors can only be performed "
                                                           "between tensors of equal shapes.\n"};
            if (t1._shape._r == -1) // case for addition between two empty tensors
                return {}; // empty tensor returned
            addComT<U, V> *_ndata = new addComT<U, V>[t1.vol];
            uint64_t counter = t1.vol;
            addComT<U, V> *nptr = _ndata;
            U *d1 = t1.data;
            V *d2 = t2.data;
            while (counter --> 0)
                *nptr++ = *d1++ + *d2++;
            return {_ndata, t1._shape._r, t1._shape._s, t1.vol, true}; // let ctor make a copy of the shape array
        }
        template <Numeric U, Numeric V>
        friend tensor<subComT<U, V>> operator-(const tensor<U> &t1, const tensor<V> &t2) {
            if (t1._shape != t2._shape)
                throw exceptions::dimension_mismatch_error{"Error: subtraction between tensors can only be performed "
                                                           "between tensors of equal shapes.\n"};
            if (t1._shape._r == -1)
                return {};
            subComT<U, V> *_ndata = new subComT<U, V>[t1.vol];
            uint64_t counter = t1.vol;
            subComT<U, V> *nptr = _ndata;
            U *d1 = t1.data;
            V *d2 = t2.data;
            while (counter --> 0)
                *nptr++ = *d1++ - *d2++;
            return {_ndata, t1._shape._r, t1._shape._s, t1.vol, true};
        }
        template <Numeric U, Numeric V>
        friend tensor<mulComT<U, V>> operator*(const U &scalar, const tensor<V> &tens) {
            if (tens._shape._r == -1)
                return {};
            mulComT<U, V> *_ndata = new mulComT<U, V>[tens.vol];
            uint64_t counter = tens.vol;
            mulComT<U, V> *nptr = _ndata;
            V *vptr = tens.data;
            while (counter --> 0)
                *nptr++ = scalar*(*vptr++);
            return {_ndata, tens._shape._r, tens._shape._s, tens.vol, true};
        }
        template <Numeric U, Numeric V>
        friend tensor<divComT<U, V>> operator/(const tensor<U> &tens, const V &scalar) {
            if (tens._shape._r == -1)
                return {};
            divComT<U, V> *_ndata = new divComT<U, V>[tens.vol];
            uint64_t counter = tens.vol;
            divComT<U, V> *nptr = _ndata;
            U *vptr = tens.data;
            while (counter --> 0)
                *nptr++ = (*vptr++)/scalar;
            return {_ndata, tens._shape._r, tens._shape._s, tens.vol, true};
        }
        template <Numeric U, Numeric V>
        friend tensor<mulComT<U, V>> tens_ops::hadamard(const tensor<U>&, const tensor<V>&);
        template <Numeric U>
        friend std::ostream &operator<<(std::ostream&, const tensor<U>&);
        template <Numeric U, Numeric V>
        friend bool operator==(const tensor<U>&, const tensor<V>&);
        template <Numeric U, Numeric V>
        friend bool operator!=(const tensor<U>&, const tensor<V>&);
        template <Numeric>
        friend class tensor; // so all `tensor<T>` classes are friends with each other
        template <Numeric>
        friend class matrix; // have to also declare these so different types can be used with each other in functions
        template <Numeric>
        friend class vector;
        template <Numeric>
        friend class ffnn;
    };
    template <Numeric T>
    std::ostream &operator<<(std::ostream &os, const tensor<T> &tens) {
        if (tens._shape._r == -1)
            return os << "[]";
        if (!tens._shape._r)
            return os << '[' << *tens.data << ']';
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
    namespace tens_ops {
        template <Numeric U, Numeric V>
        tensor<mulComT<U, V>> hadamard(const tensor<U> &t1, const tensor<V> &t2) {
            if (t1._shape != t2._shape)
                throw exceptions::dimension_mismatch_error{"Error: the Hadamard product can only be performed between "
                                                           "two tensors of equal shapes.\n"};
            if (t1._shape._r == -1) // case for hadamard product between two empty tensors
                return {}; // return empty tensor
            mulComT<U, V> *_ndata = new mulComT<U, V>[t1.vol];
            uint64_t counter = t1.vol;
            mulComT<U, V> *nptr = _ndata;
            U *d1 = t1.data;
            V *d2 = t2.data;
            while (counter --> 0)
                *nptr++ = (*d1++)*(*d2++);
            return {_ndata, t1._shape._r, t1._shape._s, t1.vol, true};
        }
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
