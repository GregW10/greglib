#ifndef GREGVCT_HPP
#define GREGVCT_HPP

#include "gregmtx.hpp"

namespace gml {
    template <Numeric T>
    class vector_base : public matrix<T> { // abstract base class for vectors
    protected:
        vector_base(T *_data, int64_t rank, uint64_t *_s, bool copy_shape = true) :
        matrix<T>{_data, rank, _s, copy_shape} {}
        vector_base(T *_data, int64_t rank, uint64_t *_s, uint64_t _volume, bool copy_shape = true) :
        matrix<T>{_data, rank, _s, _volume, copy_shape} {}
        vector_base() = default; // protected as it only serves to delay construction in derived classes
    public:
        vector_base(const vector_base<T> &other) : matrix<T>{other} {}
        vector_base(vector_base<T> &&other) : matrix<T>{std::move(other)} {}
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
        template <Numeric U, Numeric V>
        friend mulComT<U, V> dot(const vector_base<U>&, const vector_base<V>&);
        template <Numeric>
        friend class vector_base;
    };
    template <Numeric U, Numeric V>
    mulComT<U, V> dot(const vector_base<U> &rvec, const vector_base<V> &cvec) {
        /* Dot-product between two vectors. */
        if (rvec.vol != cvec.vol)
            throw exceptions::dimension_mismatch_error{"Error: the dimensionality of the vectors do not match.\n"};
        uint64_t _i = rvec.vol;
        U *rptr = rvec.data;
        V *cptr = cvec.data;
        mulComT<U, V> _dot{};
        while (_i --> 0)
            _dot += (*rptr++)*(*cptr++);
        return _dot;
    }
    template <Numeric T>
    class vector : public vector_base<T> {
        bool _cv = true; // boolean indicating whether the vector is a column vector or not
        // having the above boolean is not 100% necessary, as the property can be checked, but for a given `vector<T>`
        // object, it can only be changed through the `transpose` method, so having a boolean avoids constant checking
        void init_empty() {
            vector_base<T>::_shape._r = 2;
            vector_base<T>::_shape._s = new uint64_t[2]{};
            if (_cv) { // even for empty vectors, one of their dimensions MUST be 1 (and the other 0)
                *(vector_base<T>::_shape._s + 1) = 1;
                vector_base<T>::_shape.sizes = new uint64_t[1]{1}; // one element per row, even if zero rows
            } else {
                *vector_base<T>::_shape._s = 1;
                vector_base<T>::_shape.sizes = new uint64_t[1]{}; // 0 elements per row, as there are zero columns
            }
            // vector_base<T>::_shape.calc_sizes();
            vector_base<T>::vol = 0;
            vector_base<T>::data = new T[0];
        }
        vector(T *_data, int64_t rank, uint64_t *_s, bool copy_shape = true, bool is_col_vec = true) :
        vector_base<T>{_data, rank, _s, copy_shape}, _cv{is_col_vec} {}
        vector(T *_data, int64_t rank, uint64_t *_s, uint64_t _volume, bool copy_shape = true, bool is_col_vec = true) :
        vector_base<T>{_data, rank, _s, _volume, copy_shape}, _cv{is_col_vec} {}
    public:
        vector() {
            init_empty();
        }
        explicit vector(bool is_col_vec) : _cv{is_col_vec} {
            init_empty();
        }
        explicit vector(uint64_t _dim, bool is_col_vec = true) : _cv{is_col_vec} { // creates a vector full of zeros
            if (!_dim) {
                init_empty();
                return;
            }
            vector_base<T>::_shape._r = 2;
            if (is_col_vec) {
                vector_base<T>::_shape._s = new uint64_t[2]{_dim, 1};
                vector_base<T>::_shape.sizes = new uint64_t[1]{1};
            } else {
                vector_base<T>::_shape._s = new uint64_t[2]{1, _dim};
                vector_base<T>::_shape.sizes = new uint64_t[1]{_dim};
            }
            this->vol = _dim;
            this->data = new T[_dim]{}; // full of zeros
        }
        vector(const std::initializer_list<std::initializer_list<T>> &values) {
            typename std::initializer_list<std::initializer_list<T>>::size_type _lsize;
            if (!(_lsize = values.size())) {
                init_empty();
                return;
            }
            const std::initializer_list<T> *ptr = values.begin();
            if (_lsize == 1) {
                typename std::initializer_list<T>::size_type _size;
                if (!(_size = ptr->size())) {
                    init_empty();
                    return;
                }
                if (_size == 1) { // column vectors are preferred, so it is interpreted as that when it is ambiguous
                    vector_base<T>::_shape._s = new uint64_t[2]{1, 1};
                    vector_base<T>::_shape.sizes = new uint64_t[1]{1};
                } else {
                    _cv = false; // if init. list contains only one init. list that is larger than 1, must be a row vec.
                    vector_base<T>::_shape._s = new uint64_t[2]{1, _size};
                    vector_base<T>::_shape.sizes = new uint64_t[1]{_size};
                }
                this->data = new T[_size];
                gen::memcopy(this->data, ptr->begin(), sizeof(T), _size);
                vector_base<T>::vol = _size;
            } else {
                typename std::initializer_list<std::initializer_list<T>>::size_type counter = 0;
                vector_base<T>::data = new T[_lsize];
                T *dptr = vector_base<T>::data;
                while (counter++ < _lsize) {
                    if (ptr->size() != 1) {
                        delete [] vector_base<T>::data;
                        throw exceptions::dimension_mismatch_error{"Error: all rows of a column vector must have 1 "
                                                                   "element.\n"};
                    }
                    *dptr++ = *ptr++->begin();
                }
                vector_base<T>::_shape._s = new uint64_t[2]{_lsize, 1};
                vector_base<T>::_shape.sizes = new uint64_t[1]{1};
                vector_base<T>::vol = _lsize;
            }
            vector_base<T>::_shape._r = 2;
        }
        vector(const T *_data, uint64_t _size, bool is_col_vec = true) : _cv{is_col_vec} {
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
            if (is_col_vec) {
                vector_base<T>::_shape._s = new uint64_t[2]{_size, 1};
                vector_base<T>::_shape.sizes = new uint64_t[1]{1};
            } else {
                vector_base<T>::_shape._s = new uint64_t[2]{1, _size};
                vector_base<T>::_shape.sizes = new uint64_t[1]{_size};
            }
            vector_base<T>::vol = _size;
            vector_base<T>::data = new T[_size];
            gen::memcopy(vector_base<T>::data, _data, sizeof(T), _size);
        }
        vector(const vector<T> &other) : vector_base<T>{other}, _cv{other._cv} {}
        vector(vector<T> &&other) noexcept : vector_base<T>{std::move(other)}, _cv{other._cv} {}
        vector &transpose() noexcept /* override */ { // `override` is commented until I write the method in `matrix<T>`
            if (_cv) {
                *vector_base<T>::_shape._s = 1;
                *(vector_base<T>::_shape._s + 1) = vector_base<T>::vol;
                *vector_base<T>::_shape.sizes = vector_base<T>::vol;
            } else {
                *vector_base<T>::_shape._s = vector_base<T>::vol;
                *(vector_base<T>::_shape._s + 1) = 1;
                *vector_base<T>::_shape.sizes = 1;
            }
            _cv = !_cv; // transforms a column to row vector or vice-versa
            return *this;
        }
        vector &reshape(const tensor_shape &new_shape) override {
            if (new_shape._r != 2)
                throw std::invalid_argument{"Error: a vector must have a rank of 2.\n"};
            if (*(new_shape._s + 1) == 1)
                _cv = true;
            else if (*new_shape._s == 1) // `else if` because column vectors are favoured
                _cv = false;
            else
                throw exceptions::dimension_mismatch_error{"Error: a vector must have a shape of either "
                                                           "(N,1) or (1,N).\n"};
            /* if constexpr (isColVec) {
                if (*(new_shape._s + 1) != 1)
                    throw std::invalid_argument{"Error: the shape of a column vector must be (N,1)"};
            } else {
                if (*new_shape._s != 1)
                    throw std::invalid_argument{"Error: the shape of a row vector must be (1,N)"};
            } */
            return dynamic_cast<vector<T>&>(tensor<T>::reshape(new_shape));
        }
        vector<T> copy() const { // again, can't do anything about the method hiding :(
            return {*this};
        }
        bool is_colvec() const noexcept {
            return _cv;
        }
        bool is_rowvec() const noexcept {
            return !_cv;
        }
        using matrix<T>::at; // again, as in `gregmtx.hpp`, have to bring `at` methods into scope due to name hiding
        T &at(const std::initializer_list<uint64_t> &indices) override {
#define AT_NC \
            typename std::initializer_list<uint64_t>::size_type _size = indices.size(); \
            if (_size > 2 || !_size) \
                throw exceptions::invalid_indexing_error{"Error: a vector can be indexed by only 1 or 2 indices.\n"}; \
            if (_size == 1) \
                return *(vector_base<T>::data + *indices.begin()); \
            if (_cv) { \
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
            if (_cv) { \
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
    private:
        template <Numeric U> requires (std::is_convertible_v<mulComT<T, U>, T>)
        vector<T> &apply_cv(const matrix<U> &mat) {
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
        template <Numeric U> requires (std::is_convertible_v<mulComT<T, U>, T>)
        vector<T> &apply_rv(const matrix<U> &mat) {
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
    public:
        template <Numeric U> requires (std::is_convertible_v<mulComT<T, U>, T>)
        vector &apply(const matrix<U> &mat) {
            if (_cv)
                return this->apply_cv(mat);
            return this->apply_rv(mat);
        }
        vector &operator=(const vector<T> &other) {
            _cv = other._cv;
            return static_cast<vector<T>&>(matrix<T>::operator=(other));
        }
        vector &operator=(vector<T> &&other) {
            _cv = other._cv;
            return static_cast<vector<T>&>(matrix<T>::operator=(std::move(other)));
        }
        vector &operator=(const matrix<T> &other) { /*
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
            // EX_MAT_EQ("matrix") */
#define EX_MAT_EQ(text) \
            if (*(other._shape._s + 1) == 1) \
                _cv = true; \
            else if (*other._shape._s == 1) \
                _cv = false; \
            else \
                throw exceptions::dimension_mismatch_error{"Error: for a " text " to be assigned to a vector, at least"\
                                                           " one of its dimensions must be 1.\n"};
            return static_cast<vector<T>&>(matrix<T>::operator=(other));
        }
        vector &operator=(matrix<T> &&other) {
            EX_MAT_EQ("matrix")
            return static_cast<vector<T>&>(matrix<T>::operator=(std::move(other)));
        }
        vector &operator=(const tensor<T> &other) {
            if (other._shape._r != 2)
                throw exceptions::dimension_mismatch_error{"Error: a tensor of rank greater than 2 cannot be assigned "
                                                           "to a vector.\n"};
            EX_MAT_EQ("2D tensor")
            return static_cast<vector<T>&>(tensor<T>::operator=(other));
        }
        vector &operator=(tensor<T> &&other) {
            if (other._shape._r != 2)
                throw exceptions::dimension_mismatch_error{"Error: a tensor of rank greater than 2 cannot be assigned "
                                                           "to a vector.\n"};
            EX_MAT_EQ("2D tensor")
            return static_cast<vector<T>&>(tensor<T>::operator=(std::move(other)));
        }
        vector<T> &operator+=(const tensor<T> &other) override {
            return static_cast<vector<T>&>(tensor<T>::operator+=(other));
        }
        vector<T> &operator-=(const tensor<T> &other) override {
            return static_cast<vector<T>&>(tensor<T>::operator-=(other));
        }
        template <Numeric U, Numeric V>
        friend vector<addComT<U, V>> operator+(const vector<U> &vec1, const vector<V> &vec2) {
            if (vec1._cv != vec2._cv || vec1.vol != vec2.vol)
                throw exceptions::dimension_mismatch_error{"Error: addition between vectors can only be performed "
                                                           "between vectors of equal shapes.\n"};
            if (vec1._shape._r == -1) // would only happen if vector has been moved from
                return {};
            addComT<U, V> *_ndata = new addComT<U, V>[vec1.vol];
            uint64_t counter = vec1.vol;
            addComT<U, V> *nptr = _ndata;
            U *d1 = vec1.data;
            V *d2 = vec2.data;
            while (counter --> 0)
                *nptr++ = *d1++ + *d2++;
            return {_ndata, vec1._shape._r, vec1._shape._s, vec1.vol, true}; // let ctor make a copy of the shape array
        }
        template <Numeric U, Numeric V>
        friend vector<subComT<U, V>> operator-(const vector<U> &vec1, const vector<V> &vec2) {
            if (vec1._cv != vec2._cv || vec1.vol != vec2.vol)
                throw exceptions::dimension_mismatch_error{"Error: subtraction between vectors can only be performed "
                                                           "between vectors of equal shapes.\n"};
            if (vec1._shape._r == -1)
                return {};
            subComT<U, V> *_ndata = new subComT<U, V>[vec1.vol];
            uint64_t counter = vec1.vol;
            subComT<U, V> *nptr = _ndata;
            U *d1 = vec1.data;
            V *d2 = vec2.data;
            while (counter --> 0)
                *nptr++ = *d1++ - *d2++;
            return {_ndata, vec1._shape._r, vec1._shape._s, vec1.vol, true};
        }
        template <Numeric U, Numeric V>
        friend vector<mulComT<U, V>> operator*(const U &scalar, const vector<V> &vec) {
            if (vec._shape._r == -1)
                return {};
            mulComT<U, V> *_ndata = new mulComT<U, V>[vec.vol];
            uint64_t counter = vec.vol;
            mulComT<U, V> *nptr = _ndata;
            V *vptr = vec.data;
            while (counter --> 0)
                *nptr++ = scalar*(*vptr++);
            return {_ndata, vec._shape._r, vec._shape._s, vec.vol, true};
        }
        template <Numeric U, Numeric V>
        friend vector<divComT<U, V>> operator/(const vector<U> &vec, const V &scalar) {
            if (vec._shape._r == -1)
                return {};
            divComT<U, V> *_ndata = new divComT<U, V>[vec.vol];
            uint64_t counter = vec.vol;
            divComT<U, V> *nptr = _ndata;
            U *vptr = vec.data;
            while (counter --> 0)
                *nptr++ = (*vptr++)/scalar;
            return {_ndata, vec._shape._r, vec._shape._s, vec.vol, true};
        }
#undef EX_MAT_EQ
        template <Numeric U, Numeric V>
        friend vector<mulComT<U, V>> matcvecmul(const matrix<U> &mat, const vector<V> &vec) {
            /* This function provides dedicated multiplication between an `N x M` matrix and an `M x 1` vector. I could
             * not define this as another `operator*` overload because of the exception I have to throw if `vec` is a
             * row vector, which would make any multiplication between an `N x 1` matrix and `1 x M` vector impossible,
             * since that would require returning a `matrix<T>` with shape `N x M`, and C++ is statically typed. */
            if (!vec._cv)
                throw exceptions::matmul_error{"Error: column vector expected.\n"};
            if (*(mat._shape._s + 1) != *vec._shape._s)
                throw exceptions::dimension_mismatch_error{"Error: number of columns in matrix does not match "
                                                           "dimensionality of vector.\n"};
            uint64_t _n = *mat._shape._s; // dimensionality of new vector
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
        friend vector<mulComT<U, V>> tens_ops::hadamard(const vector<U>&, const vector<V>&);
        template <Numeric>
        friend class vector;
        template <Numeric>
        friend class ffnn;
    };
    namespace tens_ops {
        template <Numeric U, Numeric V>
        vector<mulComT<U, V>> hadamard(const vector<U> &vec1, const vector<V> &vec2) {
            if (vec1._cv != vec2._cv || vec1.vol != vec2.vol)
                throw exceptions::dimension_mismatch_error{"Error: the Hadamard product can only be performed between "
                                                           "two vectors of equal shapes.\n"};
            if (vec1._shape._r == -1) // case for hadamard product between two empty vectors
                return {}; // return empty vector
            mulComT<U, V> *_ndata = new mulComT<U, V>[vec1.vol];
            uint64_t counter = vec1.vol;
            mulComT<U, V> *nptr = _ndata;
            U *d1 = vec1.data;
            V *d2 = vec2.data;
            while (counter --> 0)
                *nptr++ = (*d1++)*(*d2++);
            return {_ndata, vec1._shape._r, vec1._shape._s, vec1.vol, true, vec1._cv};
        }
    }
}
#endif
