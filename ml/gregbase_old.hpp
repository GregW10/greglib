#include <iostream>
#include <vector>
#include <concepts>
#include <initializer_list>
#include <sstream>
#include <limits>

namespace gml {
    template <typename T>
    concept Numeric = requires (T value) {
        T{1};
        T{1.0l};
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
    typedef unsigned long long ull_t;
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
        char *o1 = static_cast<char*>(obj1);
        char *o2 = static_cast<char*>(obj2);
        size_t size = sizeof(T);
        char temp;
        while (size --> 0) {
            temp = *o1;
            *o1++ = *o2;
            *o2++ = temp;
        }
        return true;
    }
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
            out_of_bounds_error(const char *msg) : invalid_indexing_error{msg} {}

        };
    }
    template <Numeric T>
    class tensor {
        T *data;
    public:
        class shape {
        public:
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
                ull_t offset = *(ptr + _r - 1);
                ull_t *size = sizes;
                unsigned int counter = 1;
                while (counter++ < _r)
                    offset += (*ptr++)*(*size++);
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
        protected:
            // shape() {} // might need a dead ctor
        public:
            shape() : _r{0}, _s{new ull_t[0]} {}
            explicit shape(unsigned int rank) : _r{rank}, _s{new ull_t[rank]} {calc_sizes();}
            shape(unsigned int rank, ull_t *shape) : _r{rank}, _s{new ull_t[rank]} {
                memcopy(_s, shape, sizeof(ull_t), rank);
                calc_sizes();
            }
            shape(const std::initializer_list<ull_t> &shape) : _r{(unsigned int) shape.size()},
                                                                            _s{new ull_t[_r]} {
                memcopy(_s, shape.begin(), sizeof(ull_t), _r);
                calc_sizes();
            }
            shape(const shape &other) : _r{other._r}, _s{new ull_t[_r]} {
                memcopy(_s, other._s, sizeof(ull_t), _r);
                if (_r)
                    memcopy(sizes, other.sizes, sizeof(ull_t), _r - 1);
            }
            shape(shape &&other) noexcept : _r{other._r}, _s{other._s}, sizes{other.sizes} {
                other._s = nullptr;
                other.sizes = nullptr;
            }
            unsigned int rank() const noexcept {
                return _r;
            } /* I'll come back to this one, although I have realised it's impossible without recursion (so, not eff.)
            template <typename ... Args>
            ull_t index_at(Args ...args) const noexcept {
                if constexpr (!sizeof...(args))
                    return 0;
                if (sizeof...(args) != _r) {
                    constexpr size_t buff_size = 138 + std::numeric_limits<size_t>::digits10 +
                            std::numeric_limits<unsigned int>::digits10; // compile-time constant
                    char buff[buff_size]{};
                    snprintf(buff, buff_size, "Error: invalid use of function template \"index_at\" with %zu parameters"
                                              " (indices) on a tensor with rank %u (rank must match no. of indices).\n",
                                              sizeof...(args), _r);
                    throw std::invalid_argument{buff};
                }

            } */
            shape copy() const {
                return {*this};
            }
            shape &operator=(const shape &other) {
                if (&other == this)
                    return *this;
                if (other._r == this->_r) {
                    memcopy(_s, other._s, sizeof(ull_t), _r);
                    if (_r)
                        memcopy(sizes, other.sizes, sizeof(ull_t), _r - 1);
                    return *this;
                }
                this->_r = other._r;
                this->~shape();
                this->_s = new ull_t[this->_r]{};
                if (this->_r) {
                    this->sizes = new ull_t[this->_r - 1];
                    memcopy(this->_s, other._s, sizeof(ull_t), this->_r);
                    memcopy(this->sizes, other.sizes, sizeof(ull_t), this->_r - 1);
                }
                else
                    this->sizes = nullptr;
                return *this;
            }
            shape &operator=(shape &&other) noexcept {
                if (&other == this)
                    return *this;
                this->_r = other._r;
                this->~shape();
                this->_s = other._s;
                other._s = nullptr;
                this->sizes = other.sizes;
                other.sizes = nullptr;
                return *this;
            }
            ~shape() {
                delete [] _s;
                delete [] sizes;
            }
            friend std::ostream &operator<<(std::ostream &out, const shape &s) {
                if (!s._r)
                    return out << "()";
                unsigned int counter = 0;
                ull_t *ptr = s._s;
                out << '(';
                goto start;
                while (counter < s._r) {
                    out << ", ";
                    start:
                    out << *ptr++;
                    ++counter;
                }
                return out << ')';
            }
        };
    public:
        tensor() : data{1, 0} {
            // shape = new unsigned long long{1};

        }
    };
}
