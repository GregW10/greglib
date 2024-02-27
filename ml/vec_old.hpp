#ifndef GREGVCT_HPP
#define GREGVCT_HPP

#include "gregmtx.hpp"

namespace gml {
    template <Numeric T>
    class vector_base : public matrix<T> {
    public:
        vector_base &reshape(const tensor_shape &new_shape) override {
            /* For a vector, `reshape` can, essentially, only be used to transpose the vector to a row/column vector. */
            if (this->vol != new_shape.volume())
                throw exceptions::dimension_mismatch_error{"Error: cannot reshape to requested shape as this would "
                                                           "require modifying the number of elements.\n"};
            if (new_shape._r != 2)
                throw exceptions::dimension_mismatch_error{"Error: a matrix must have a rank of 2.\n"};
            if (*new_shape._s != 1 && *(new_shape._s + 1) != 1) {
                throw exceptions::dimension_mismatch_error{"Error: a vector must have one of its two dimensions "
                                                           "equalling 1.\n"};
            }
            this->_shape = new_shape;
            return *this;
        }
        T &at(const std::initializer_list<uint64_t> &indices) override {
#define AT_NC \
            typename std::initializer_list<uint64_t>::size_type _size = indices.size(); \
            if (_size > 2 || !_size) \
                throw exceptions::invalid_indexing_error{"Error: a vector can be indexed with only 1 or 2 indices.\n"};\
            uint64_t index{}; \
            if (_size == 2) { \
                uint64_t *ptr = indices.begin(); \
                if (*ptr == 1) \
                    index = *(ptr + 1); \
                else if (*(ptr + 1) == 1) \
                    index = *ptr; \
                else \
                    throw exceptions::invalid_indexing_error{"Error: cannot index a vector with two non-1 indices.\n"};\
            } else \
                index = *indices.begin(); \
            return *(matrix<T>::data + index);
            AT_NC
        }// check if 2, then mat, else, check == 1 in one
        const T &at(const std::initializer_list<uint64_t> &indices) const override {
            AT_NC
        }
        T &operator()(uint64_t _i, uint64_t _j) noexcept override { // doesn't throw, so must be careful when called
            return *(tensor<T>::data + _i*(*(tensor<T>::_shape._s + 1)) + _j);
        }
        const T &operator()(uint64_t _i, uint64_t _j) const noexcept override {
            return *(tensor<T>::data + _i*(*(tensor<T>::_shape._s + 1)) + _j);
        }
        T &operator[](uint64_t index) noexcept {
            return *(matrix<T>::data + index);
        }
        const T &operator[](uint64_t index) const noexcept {
            return *(matrix<T>::data + index);
        }
    };
}
#endif
