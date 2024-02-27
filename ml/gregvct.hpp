#ifndef GREGVCT_HPP
#define GREGVCT_HPP

#include "gregmtx.hpp"

namespace gml {
    template <Numeric T>
    class vector_base : public matrix<T> { // abstract base class for vectors
    public:
        mulComT<T, T> mag_sq() const noexcept { // squared magnitude of the vector
            mulComT<T, T> _tot{};
            T *ptr = matrix<T>::data;
            uint64_t counter = matrix<T>::vol;
            while (counter --> 0) {
                _tot += (*ptr)*(*ptr);
                ++ptr;
            }
            return _tot;
        }
        long double mag() const noexcept { // magnitude of the vector (in Euclidean space)
            return sqrtl(this->mag_sq());
        }
        vector_base &normalise() noexcept {
            long double _mag = this->mag();
            T *ptr = matrix<T>::data;
            uint64_t counter = matrix<T>::vol;
            while (counter --> 0)
                *ptr++ /= _mag;
            return *this;
        }
        divComT<T, long double> unitvec() const noexcept {
            return operator/(*this, this->mag());
        }
        vector_base &reshape(const tensor_shape &new_shape) override = 0;
        T &at(const std::initializer_list<uint64_t> &indices) override = 0;
        const T &at(const std::initializer_list<uint64_t> &indices) const override = 0;
        T &operator()(uint64_t _i, uint64_t _j) noexcept override = 0;
        const T &operator()(uint64_t _i, uint64_t _j) const noexcept override = 0;
        T &operator[](uint64_t index) noexcept {
            return *(matrix<T>::data + index);
        }
        const T &operator[](uint64_t index) const noexcept {
            return *(matrix<T>::data + index);
        }
    };
    template <Numeric T, bool isColVec = true>
    class vector : virtual public vector_base<T> {
    public:
        template <Numeric U> requires (std::is_convertible_v<mulComT<T, U>, T>)
        vector &apply(const matrix<U> &mat) { // might want to consider moving this into `vector<T>` class
            /* This method transforms the vector according to the transformation encoded by `mat`. This would be the
             * equivalent of multiplying `mat*(*this)` and then reassigning the result to `*this`. */
            if (*(mat._shape._s + 1) != matrix<T>::vol)
                throw exceptions::dimension_mismatch_error{"Error: the number of columns of the matrix doesn't match "
                                                           "the dimensionality of the vector.\n"};
            if (!*mat._shape._s) {
                *matrix<T>::_shape._s = 0;
                matrix<T>::vol = 0;
                delete [] matrix<T>::data;
                matrix<T>::data = new T[0];
                return *this;
            }
            T *ndata = new T[*mat._shape._s]; // need to allocate, even if vector doesn't change size
            uint64_t _i = *mat._shape._s; // new dimensionality of vector
            uint64_t _j;
            T *ptr;
            T *nptr = ndata;
            U *mptr = mat.data;
            while (_i --> 0) {
                ptr = matrix<T>::data;
                _j = matrix<T>::vol;
                while (_j --> 0)
                    *nptr += (*mptr++)*(*ptr++);
                ++nptr;
            }
            *matrix<T>::_shape._s = *mat._shape._s;
            matrix<T>::vol = *mat._shape._s;
            delete [] matrix<T>::data;
            matrix<T>::data = ndata;
            return *this;
        }
    };
}
#endif
