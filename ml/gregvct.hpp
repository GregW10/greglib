#ifndef GREGVCT_HPP
#define GREGVCT_HPP

#include "gregmtx.hpp"

namespace gml {
    template <Numeric T>
    class vector_base : public matrix<T> { // abstract base class for vectors
        vector_base(T *_data, int64_t rank, uint64_t *_s, bool copy_shape = true) :
        matrix<T>{_data, rank, _s, copy_shape} {}
        vector_base(T *_data, int64_t rank, uint64_t *_s, uint64_t _volume, bool copy_shape = true) :
        matrix<T>{_data, rank, _s, _volume, copy_shape} {}
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
        T &operator()(uint64_t _i, uint64_t _j) override = 0;
        const T &operator()(uint64_t _i, uint64_t _j) const override = 0;
        T &operator[](uint64_t index) noexcept {
            return *(matrix<T>::data + index);
        }
        const T &operator[](uint64_t index) const noexcept {
            return *(matrix<T>::data + index);
        }
        template <Numeric U>
        friend class vector_base;
    };
    template <Numeric T, bool isColVec = true>
    class vector : virtual public vector_base<T> {
        void init_empty() {
            vector_base<T>::_shape._r = 2;
            vector_base<T>::_shape._s = new uint64_t[2]{};
            if constexpr (isColVec) // even for empty vectors, one of their dimensions MUST be 1 (and the other 0)
                *(vector_base<T>::_shape._s + 1) = 1;
            else
                *vector_base<T>::_shape._s = 1;
            vector_base<T>::_shape.calc_sizes();
            vector_base<T>::vol = 0;
            vector_base<T>::data = new T[0];
        }
        vector(T *_data, int64_t rank, uint64_t *_s, bool copy_shape = true) :
        vector_base<T>{_data, rank, _s, copy_shape} {}
        vector(T *_data, int64_t rank, uint64_t *_s, uint64_t _volume, bool copy_shape = true) :
        vector_base<T>{_data, rank, _s, _volume, copy_shape} {}
    public:
        vector() {
            init_empty();
        }
        vector(const std::initializer_list<T> &values) {
            typename std::initializer_list<T>::size_type _size;
            if (!(_size = values.size())) {
                init_empty();
                return;
            }
            vector_base<T>::_shape._r = 2;
            if constexpr (isColVec)
                vector_base<T>::_shape._s = new uint64_t[2]{_size, 1};
            else
                vector_base<T>::_shape._s = new uint64_t[2]{1, _size};
            vector_base<T>::_shape.calc_sizes();
            vector_base<T>::vol = vector_base<T>::_shape.volume();
            vector_base<T>::data = new T[vector_base<T>::vol];
            gen::memcopy(vector_base<T>::data, values.begin(), sizeof(T), _size);
        }
        vector(const T *_data, uint64_t _size) {
            if (!_data) {
                if (_size)
                    throw std::invalid_argument{"Error: nullptr passed as data but with non-zero size.\n"};
                init_empty();
                return;
            }
            if (!_size) {
                init_empty();
                return;
            }
            vector_base<T>::_shape._r = 2;
            if constexpr (isColVec)
                vector_base<T>::_shape._s = new uint64_t[2]{_size, 1};
            else
                vector_base<T>::_shape._s = new uint64_t[2]{1, _size};
            vector_base<T>::_shape.calc_sizes();
            vector_base<T>::vol = vector_base<T>::_shape.volume();
            vector_base<T>::data = new T[vector_base<T>::vol];
            gen::memcopy(vector_base<T>::data, _data, sizeof(T), _size);
        }
        vector(const vector<T, isColVec> &other) : matrix<T>{other} {}
        vector(vector<T, isColVec> &&other) noexcept : matrix<T>{std::move(other)} {}
        /* vector &transpose() noexcept {
            uint64_t _temp = *vector_base<T>::_shape._s;
            *vector_base<T>::_shape._s = *(vector_base<T>::_shape._s + 1);
            *(vector_base<T>::_shape._s + 1) = _temp;
            return *this;
        } */ // gonna really have to think about this, and maybe even rethink using boolean template arguments
        vector &reshape(const tensor_shape &new_shape) override {
            if (new_shape._r != 2)
                throw std::invalid_argument{"Error: a vector must have a rank of 2.\n"};
            if constexpr (isColVec) {
                if (*(new_shape._s + 1) != 1)
                    throw std::invalid_argument{"Error: the shape of a column vector must be (N,1)"};
            } else {
                if (*new_shape._s != 1)
                    throw std::invalid_argument{"Error: the shape of a row vector must be (1,N)"};
            }
            return dynamic_cast<vector<T, isColVec>&>(tensor<T>::reshape(new_shape));
        }
        T &at(const std::initializer_list<uint64_t> &indices) override {
#define AT_NC \
            typename std::initializer_list<uint64_t>::size_type _size = indices.size(); \
            if (_size > 2 || !_size) \
                throw exceptions::invalid_indexing_error{"Error: a vector can be indexed by only 1 or 2 indices.\n"}; \
            if (_size == 1) \
                return *(vector_base<T>::data + *indices.begin()); \
            if constexpr (isColVec) { \
                if (*(indices.begin() + 1)) \
                    throw exceptions::dimension_mismatch_error{"Error: the second index for indexing a column vector " \
                                                               "must be 0.\n"}; \
                return *(vector_base<T>::data + *indices.begin()); \
            } \
            if (*indices.begin()) \
                throw exceptions::dimension_mismatch_error{"Error: the first index for indexing a row vector must be " \
                                                           "0.\n"}; \
            return *(vector_base<T>::data + *(indices.begin() + 1));
            AT_NC
        }
        const T &at(const std::initializer_list<uint64_t> &indices) const override {
            AT_NC
        }
#undef AT_NC
        T &operator()(uint64_t _i, uint64_t _j) override {
#define OP_NC \
            if constexpr (isColVec) { \
                if (_j) \
                    throw exceptions::dimension_mismatch_error{"Error: the second index for indexing a column vector " \
                                                               "must be 0.\n"}; \
                return *(vector_base<T>::data + _i); \
            } \
            if (_i) \
                throw exceptions::dimension_mismatch_error{"Error: the first index for indexing a row vector must be " \
                                                           "0.\n"}; \
            return *(vector_base<T>::data + _j);
            OP_NC
        }
        const T &operator()(uint64_t _i, uint64_t _j) const override {
            OP_NC
        }
#undef OP_NC
        template <Numeric U> requires (std::is_convertible_v<mulComT<T, U>, T> && isColVec)
        vector &apply(const matrix<U> &mat) {
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
        template <Numeric U> requires (std::is_convertible_v<mulComT<T, U>, T> && !isColVec)
        vector &apply(const matrix<U> &mat) {
            /* This overload of `apply` is the equivalent of the above, but for row vectors. It is the equivalent of
             * multiplying `(*this)*mat` and then reassigning the result to `*this`. */
            if (*mat._shape._s != matrix<T>::vol)
                throw exceptions::dimension_mismatch_error{"Error: the number of rows of the matrix doesn't match the "
                                                           "dimensionality of the vector.\n"};
            uint64_t _m; // new dimensionality of vector (equal to number of columns in `mat`)
            if (!(_m = *(mat._shape._s + 1))) {
                *(matrix<T>::_shape._s + 1) = 0;
                matrix<T>::vol = 0;
                delete [] matrix<T>::data;
                matrix<T>::data = new T[0];
                return *this;
            }
            T *ndata = new T[_m]{}; // need to allocate, even if vector doesn't change size
            uint64_t _i = _m; // new dimensionality of vector
            uint64_t _j;
            T *ptr;
            T *nptr = ndata;
            U *bmptr = mat.data;
            U *mptr;
            while (_i --> 0) {
                ptr = matrix<T>::data;
                _j = matrix<T>::vol;
                mptr = bmptr++;
                while (_j --> 0) {
                    *nptr += (*ptr++)*(*mptr);
                    mptr += _m;
                }
                ++nptr;
            }
            *(matrix<T>::_shape._s + 1) = _m;
            matrix<T>::vol = _m;
            delete [] matrix<T>::data;
            matrix<T>::data = ndata;
            return *this;
        }
        vector &operator=(const vector<T, isColVec> &other) {
            return dynamic_cast<vector<T, isColVec>&>(matrix<T>::operator=(other));
        }
        vector &operator=(vector<T, isColVec> &&other) {
            return dynamic_cast<vector<T, isColVec>&>(matrix<T>::operator=(std::move(other)));
        }
        vector &operator=(const matrix<T> &other) {
#define EX_MAT_EQ(text) \
            if constexpr (isColVec) { \
                if (*(other._shape._s + 1) != 1) \
                    throw exceptions::dimension_mismatch_error{"Error: a " text " with more than one column cannot be "\
                                                               "assigned to a column vector.\n"}; \
            } else { \
                if (*other._shape._s != 1) \
                    throw exceptions::dimension_mismatch_error{"Error: a " text " with more than one row cannot be " \
                                                               "assigned to a row vector.\n"}; \
            }
            EX_MAT_EQ("matrix")
            return dynamic_cast<vector<T, isColVec>&>(matrix<T>::operator=(other));
        }
        vector &operator=(matrix<T> &&other) {
            EX_MAT_EQ("matrix")
            return dynamic_cast<vector<T, isColVec>&>(matrix<T>::operator=(std::move(other)));
        }
        vector &operator=(const tensor<T> &other) {
            if (other._shape._r != 2)
                throw exceptions::dimension_mismatch_error{"Error: a tensor of rank greater than 2 cannot be assigned "
                                                           "to a vector.\n"};
            EX_MAT_EQ("2D tensor")
            return dynamic_cast<vector<T, isColVec>&>(tensor<T>::operator=(other));
        }
        vector &operator=(tensor<T> &&other) {
            if (other._shape._r != 2)
                throw exceptions::dimension_mismatch_error{"Error: a tensor of rank greater than 2 cannot be assigned "
                                                           "to a vector.\n"};
            EX_MAT_EQ("2D tensor")
            return dynamic_cast<vector<T, isColVec>&>(tensor<T>::operator=(std::move(other)));
        }
#undef EX_MAT_EQ
        template <Numeric U, Numeric V, bool I>
        friend vector<mulComT<U, V>, I> operator*(const matrix<U> &mat, const vector<V, I> &vec) {
            if (*(mat._shape._s + 1) != vec.vol)
                throw exceptions::dimension_mismatch_error{"Error: number of columns in matrix does not match "
                                                           "dimensionality of vector.\n"};
            uint64_t _n = *mat._shape._s; // new dimensionality of vector
            if (!_n)
                return {}; // return empty vector if multiplying by empty matrix
            uint64_t _i = _n;
            uint64_t _j;
            mulComT<U, V> *ndata = new mulComT<U, V>[_n]{};
            mulComT<U, V> *nptr = ndata;
            U *mptr = mat.data;
            V *vptr;
            while (_i --> 0) {
                _j = vec.vol;
                vptr = vec.data;
                while (_j --> 0)
                    *nptr += (*mptr++)*(*vptr++);
                ++nptr;
            }
            uint64_t *nshape = new uint64_t[2]{_n, 1};
            return {ndata, 2, nshape, _n, false};
        }
        template <Numeric U, Numeric V>
        friend mulComT<U, V> operator*(const vector<U, false> &rvec, const vector<V, true> &cvec) {
            if (*(rvec._shape._s + 1) != *cvec._shape._s)
                throw exceptions::dimension_mismatch_error{"Error: the dimensionality of the vectors do not match.\n"};
            uint64_t _i = *cvec._shape._s;
            U *rptr = rvec.data;
            V *cptr = cvec.data;
            mulComT<U, V> _dot{};
            while (_i --> 0)
                _dot += *(*rptr++)*(*cptr++);
            return _dot;
        };
        template <Numeric U, bool isCV>
        friend class vector;
    };
}
#endif
