#ifndef GREGMTX_HPP
#define GREGMTX_HPP

#include "gregtens.hpp"

namespace gml {
    namespace exceptions {
        class matmul_error : public tensor_error {
            char *_msg{};
        public:
            matmul_error() : tensor_error{"Error: matrices cannot be multiplied together.\n"} {}
            explicit matmul_error(const char *msg) : tensor_error{msg} {}
            matmul_error(const tensor_shape &s1, const tensor_shape &s2) {
                std::ostringstream oss;
                oss << "Error: matrix with shape " << s1 << " cannot be multiplied by matrix with shape " << s2
                    << ". The number of columns of the first matrix must equal the number of rows of the second.\n";
#ifndef __APPLE__
                std::string_view view = oss.view();
                typename std::string_view::size_type _size = view.size();
                _msg = new char[_size + 1];
                gen::memcopy(_msg, view.data(), sizeof(char), _size); // string view not guaranteed to have null term.
                *(_msg + _size) = 0;
#else
                _msg = new char[static_cast<std::streamsize>(oss.tellp()) + 1];
                gen::strcpy_c(_msg, oss.str().c_str()); // temporary solution until Apple fixes clang
#endif
            }
            const char *what() const noexcept override {
                if (_msg)
                    return _msg;
                return tensor_error::what();
            }
            ~matmul_error() override { // `std::exception` has a virtual destructor, thank god
                delete [] _msg;
            }
        };
    }
    template <Numeric T>
    class matrix : public tensor<T> {
        void load_mtsr(std::ifstream &in, size_t st_size) {
            uint64_t shape[2];
            in.read((char *) shape, 2*sizeof(uint64_t));
            if (!shape[0] || !shape[1]) { // case for matrix with zero size
                if (st_size > 4 + 3*sizeof(uint64_t))
                    throw exceptions::invalid_mtsr_format{"Error: invalid .mtsr format - file too large for empty "
                                                          "matrix.\n"};
                tensor<T>::_shape._r = 2;
                tensor<T>::_shape._s = new ull_t[2]{shape[0], shape[1]};
                tensor<T>::_shape.sizes = new uint64_t[1]{shape[1]};
                tensor<T>::vol = 0;
                tensor<T>::data = new T[0];
                in.close();
                return;
            }
            if (st_size <= 4 + 3*sizeof(uint64_t))
                throw exceptions::invalid_mtsr_format{"Error: invalid .mtsr format - file too small for non-empty "
                                                      "matrix.\n"};
            uint64_t T_size;
            in.read((char *) &T_size, sizeof(uint64_t));
            if (T_size != sizeof(T))
                throw exceptions::invalid_mtsr_format{"Error: reported size of data type in .mtsr file does not match "
                                                      "\"sizeof(T)\".\n"};
            tensor<T>::vol = shape[0]*shape[1];
            uint64_t fsize = 4 + 3*sizeof(uint64_t) + tensor<T>::vol*sizeof(T);
            if (st_size != fsize)
                throw exceptions::invalid_mtsr_format{"Error: invalid file size.\n"};
            tensor<T>::data = new T[tensor<T>::vol];
            in.read((char *) tensor<T>::data, tensor<T>::vol*sizeof(T));
            tensor<T>::_shape._r = 2;
            tensor<T>::_shape._s = new uint64_t[2]{shape[0], shape[1]};
            tensor<T>::_shape.sizes = new uint64_t[1]{shape[1]};
            in.close();
        }
        void delegate_mtsr(const char *path) {
            if (!path || !*path)
                throw std::invalid_argument{"Error: null or empty string passed as path.\n"};
            struct stat buff{};
            if (stat(path, &buff) == -1)
                throw std::ios_base::failure{"Error: could not access file information.\n"};
            if (buff.st_size < 8)
                throw exceptions::invalid_tsr_format{"Error: invalid file size (too small).\n"};
            std::ifstream in{path, std::ios_base::in | std::ios_base::binary};
            if (!in.good())
                throw std::ios_base::failure{"Error: could not read from .mtsr file.\n"};
            char header[5]{};
            in.read(header, 4*sizeof(char));
            if (gen::equals(header, "mtsr")) {
                if (buff.st_size < 4 + 2*sizeof(uint64_t))
                    throw exceptions::invalid_mtsr_format{"Error: invalid .mtsr format (file too small).\n"};
                load_mtsr(in, buff.st_size);
            } else if (!header[0] && gen::equals(((char *) header) + 1, "tsr")) {
                tensor<T>::load_tsr(in, buff.st_size);
                if (tensor<T>::_shape._r != 2 && tensor<T>::data)
                    throw exceptions::invalid_tsr_format{"Error: the .tsr file used to load the matrix was for an "
                                                         "N-dimensional tensor. Load the file using a "
                                                         "\"gml::tensor<T>\" object instead.\n"};
            }
            else
                throw exceptions::invalid_tsr_format{"Error: invalid .tsr format (invalid header).\n"};
        }
    protected:
        matrix(T *_data, int64_t rank, uint64_t *_s, bool copy_shape = true) :
        tensor<T>{_data, rank, _s, copy_shape} {}
        matrix(T *_data, int64_t rank, uint64_t *_s, uint64_t _volume, bool copy_shape = true) :
        tensor<T>{_data, rank, _s, _volume, copy_shape} {}
    public:
        matrix() : tensor<T>{tensor_shape{0, 0}} {
            // CHECK THIS
        }
        explicit matrix(const char *tsr_path) /* : tensor<T>{false} */ {
            delegate_mtsr(tsr_path);
        }
        matrix(uint64_t rows, uint64_t columns) : tensor<T>{tensor_shape{rows, columns}} {} // matrix full of zeros
        matrix(const std::vector<T> &_data, uint64_t rows, uint64_t columns) :
        tensor<T>{_data, tensor_shape{rows, columns}} {}
        matrix(const T *_data, uint64_t rows, uint64_t columns) : tensor<T>{_data, tensor_shape{rows, columns}} {}
        matrix(const std::initializer_list<std::initializer_list<T>> &li) /* : tensor<T>{false} */ {
            typename std::initializer_list<T>::size_type li_size = li.size();
            tensor<T>::_shape._r = 2; // rank will always be 2 for a matrix
            tensor<T>::_shape._s = new ull_t[2]{}; // default initialised to zeros
            tensor<T>::_shape.sizes = new ull_t[1]{}; // zero
            if (!li_size) { // emtpy initializer list, so no action taken
                this->data = new T[tensor<T>::vol]; // zero elements
                return; // matrix is left as a 2D empty tensor
            }
            const std::initializer_list<T> *ptr = li.begin();
            typename std::initializer_list<T>::size_type sub_elems = ptr->size();
            if (!sub_elems) { // list of empty initializer lists
                this->data = new T[tensor<T>::vol]; // zero elements
                return; // again, left as a 2D empty tensor
            }
            tensor<T>::_shape._s[0] = li_size;
            tensor<T>::_shape._s[1] = sub_elems;
            tensor<T>::vol = li_size*sub_elems;
            tensor<T>::data = new T[tensor<T>::vol];
            T *dptr = tensor<T>::data;
            while (li_size --> 0) {
                if (ptr->size() != sub_elems) {
                    delete [] tensor<T>::data; // it is necessary to free the memory before throwing
                    throw std::invalid_argument{"Error: sizes of nested initializer lists do not match.\n"};
                }
                gen::memcopy(dptr, ptr++->begin(), sizeof(T), sub_elems);
                dptr += sub_elems;
            }
            tensor<T>::_shape.sizes[0] = sub_elems;
        }
        matrix(const matrix<T> &other) : tensor<T>{other} {}
        // only when a matrix is moved from will it be allowed to be a properly empty tensor with `_r == -1`:
        matrix(matrix<T> &&other) noexcept : tensor<T>{std::move(other)} { /*
            tensor<T>::_shape._r = 2;
            tensor<T>::_shape._s = other._shape._s;
            tensor<T>::_shape.sizes = other._shape.sizes;
            tensor<T>::vol = other.vol;
            tensor<T>::data = other.data; */
        }
        template <Numeric U> requires (std::is_convertible<U, T>::value)
        matrix(const matrix<U> &other) : tensor<T>{other} {}
        matrix &reshape(const tensor_shape &new_shape) override {
            if (new_shape._r != 2)
                throw exceptions::dimension_mismatch_error{"Error: a matrix must have a rank of 2.\n"};
            // be aware that this results in an extra redundant boolean check:
            return dynamic_cast<matrix<T>&>(tensor<T>::reshape(new_shape));
        }
        using tensor<T>::at; // have to bring the `at` methods into scope due to name hiding
        T &at(const std::initializer_list<uint64_t> &indices) override {
            if (indices.size() != 2)
                throw exceptions::indexing_dimension_error{"Error: exactly two indices must be given to index a "
                                                           "matrix.\n"};
            return *(tensor<T>::data + (*indices.begin())*(*(tensor<T>::_shape._s + 1)) + *(indices.end() - 1));
        }
        const T &at(const std::initializer_list<uint64_t> &indices) const override {
            if (indices.size() != 2)
                throw exceptions::indexing_dimension_error{"Error: exactly two indices must be given to index a "
                                                           "matrix.\n"};
            return *(tensor<T>::data + (*indices.begin())*(*(tensor<T>::_shape._s + 1)) + *(indices.end() - 1));
        }
        matrix<T> copy() const { // nothing doing with the method hiding :(
            return {*this};
        }
        std::ofstream::pos_type to_mtsr(const char *path) {
            if (!path || !*path)
                return 0;
            std::ofstream out{path, std::ios_base::out | std::ios_base::trunc};
            if (!out.good())
                return 0;
            const char *header = "mtsr";
            out.write(header, 4*sizeof(char));
            out.write((char *) tensor<T>::_shape._s, tensor<T>::_shape._r*sizeof(uint64_t));
            uint64_t T_size = sizeof(T);
            out.write((char *) &T_size, sizeof(uint64_t));
            out.write((char *) tensor<T>::data, tensor<T>::vol*sizeof(T));
            std::ofstream::pos_type pos = out.tellp();
            out.close();
            return pos;
        }
        virtual T &operator()(uint64_t _i, uint64_t _j) { // `noexcept` not specified as derived classes throw here
            return *(tensor<T>::data + _i*(*(tensor<T>::_shape._s + 1)) + _j);
        }
        virtual const T &operator()(uint64_t _i, uint64_t _j) const {
            return *(tensor<T>::data + _i*(*(tensor<T>::_shape._s + 1)) + _j);
        }
        virtual matrix &operator=(const matrix<T> &other) {
            return static_cast<matrix<T>&>(tensor<T>::operator=(other));
        }
        virtual matrix &operator=(matrix<T> &&other) {
            return static_cast<matrix<T>&>(tensor<T>::operator=(std::move(other)));
        }
        matrix &operator=(const tensor<T> &other) {
            if (other._shape._r != 2)
                throw exceptions::dimension_mismatch_error{"Error: only 2D tensors can be assigned to matrices.\n"};
            return static_cast<matrix<T>&>(tensor<T>::operator=(other));
        }
        matrix &operator=(tensor<T> &&other) {
            if (other._shape._r != 2)
                throw exceptions::dimension_mismatch_error{"Error: only 2D tensors can be assigned to matrices.\n"};
            return static_cast<matrix<T>&>(tensor<T>::operator=(std::move(other)));
        }
        matrix<T> &operator+=(const tensor<T> &other) override { // can be made more efficient
            return static_cast<matrix<T>&>(tensor<T>::operator+=(other));
        }
        matrix<T> &operator-=(const tensor<T> &other) override {
            return static_cast<matrix<T>&>(tensor<T>::operator-=(other));
        }
        template <Numeric U, Numeric V>
        friend matrix<addComT<U, V>> operator+(const matrix<U>&, const matrix<V>&);
        template <Numeric U, Numeric V>
        friend matrix<subComT<U, V>> operator-(const matrix<U>&, const matrix<V>&);
        template <Numeric U, Numeric V>
        friend matrix<mulComT<U, V>> operator*(const U&, const matrix<V>&);
        template <Numeric U, Numeric V>
        friend matrix<divComT<U, V>> operator/(const matrix<U>&, const V&);
        template <Numeric U, Numeric V>
        friend matrix<mulComT<U, V>> operator*(const matrix<U>&, const matrix<V>&);
        template <Numeric U, Numeric V>
        friend vector<mulComT<U, V>> matcvecmul(const matrix<U>&, const vector<V>&);
        template <Numeric U, Numeric V>
        friend matrix<mulComT<U, V>> tens_ops::hadamard(const matrix<U>&, const matrix<V>&);
        template <Numeric>
        friend class matrix;
        template <Numeric>
        friend class vector;
        template <Numeric>
        friend class ffnn;
    };
    template <Numeric U, Numeric V>
    matrix<addComT<U, V>> operator+(const matrix<U> &mat1, const matrix<V> &mat2) {
        if (mat1._shape._r == -1) { // would only happen if the matrix has been moved from
            if (mat2._shape._r != -1)
                throw exceptions::empty_tensor_error{"Error: addition between an empty and a non-empty tensor "
                                                     "cannot be performed.\n"};
            return {};
        }
        else if (mat2._shape._r == -1)
            throw exceptions::empty_tensor_error{"Error: addition between an empty and a non-empty tensor "
                                                 "cannot be performed.\n"};
        if (*mat1._shape._s != *mat2._shape._s || *(mat1._shape._s + 1) != *(mat2._shape._s + 1))
            throw exceptions::dimension_mismatch_error{"Error: addition between matrices can only be performed "
                                                       "between matrices of equal shapes.\n"};
        addComT<U, V> *_ndata = new addComT<U, V>[mat1.vol];
        uint64_t counter = mat1.vol;
        addComT<U, V> *nptr = _ndata;
        U *d1 = mat1.data;
        V *d2 = mat2.data;
        while (counter --> 0)
            *nptr++ = *d1++ + *d2++;
        return {_ndata, mat1._shape._r, mat1._shape._s, mat1.vol, true}; // let ctor make a copy of the shape array
    }
    template <Numeric U, Numeric V>
    matrix<subComT<U, V>> operator-(const matrix<U> &mat1, const matrix<V> &mat2) {
        if (mat1._shape._r == -1) {
            if (mat2._shape._r != -1)
                throw exceptions::empty_tensor_error{"Error: addition between an empty and a non-empty tensor "
                                                     "cannot be performed.\n"};
            return {};
        }
        else if (mat2._shape._r == -1)
            throw exceptions::empty_tensor_error{"Error: addition between an empty and a non-empty tensor "
                                                 "cannot be performed.\n"};
        if (*mat1._shape._s != *mat2._shape._s || *(mat1._shape._s + 1) != *(mat2._shape._s + 1))
            throw exceptions::dimension_mismatch_error{"Error: subtraction between matrices can only be performed "
                                                       "between matrices of equal shapes.\n"};
        subComT<U, V> *_ndata = new subComT<U, V>[mat1.vol];
        uint64_t counter = mat1.vol;
        subComT<U, V> *nptr = _ndata;
        U *d1 = mat1.data;
        V *d2 = mat2.data;
        while (counter --> 0)
            *nptr++ = *d1++ - *d2++;
        return {_ndata, mat1._shape._r, mat1._shape._s, mat1.vol, true};
    }
    template <Numeric U, Numeric V>
    matrix<mulComT<U, V>> operator*(const U &scalar, const matrix<V> &mat) {
        if (mat._shape._r == -1)
            return {};
        mulComT<U, V> *_ndata = new mulComT<U, V>[mat.vol];
        uint64_t counter = mat.vol;
        mulComT<U, V> *nptr = _ndata;
        V *vptr = mat.data;
        while (counter --> 0)
            *nptr++ = scalar*(*vptr++);
        return {_ndata, mat._shape._r, mat._shape._s, mat.vol, true};
    }
    template <Numeric U, Numeric V>
    matrix<divComT<U, V>> operator/(const matrix<U> &mat, const V &scalar) {
        if (mat._shape._r == -1)
            return {};
        divComT<U, V> *_ndata = new divComT<U, V>[mat.vol];
        uint64_t counter = mat.vol;
        divComT<U, V> *nptr = _ndata;
        U *vptr = mat.data;
        while (counter --> 0)
            *nptr++ = (*vptr++)/scalar;
        return {_ndata, mat._shape._r, mat._shape._s, mat.vol, true};
    }
    template <Numeric U, Numeric V>
    matrix<mulComT<U, V>> operator*(const matrix<U> &m1, const matrix<V> &m2) {
        uint64_t _n = *m1._shape._s; // num. rows of resulting matrix
        uint64_t _m = *(m2._shape._s + 1); // num. columns of resulting matrix
        uint64_t _c = *m2._shape._s; // common dimension
        if (*(m1._shape._s + 1) != _c) // case for undefined matrix multiplication
            throw exceptions::matmul_error{m1._shape, m2._shape};
        uint64_t _volume = _n*_m;
        mulComT<U, V> *pdata = new mulComT<U, V>[_volume]{};
        mulComT<U, V> *dptr = pdata;
        U *dpm1 = m1.data;
        U *m1_inner{};
        V *dpm2;
        V *m2_inner{};
        uint64_t _i = 0;
        uint64_t _j;
        uint64_t _k;
        for (; _i < _n; ++_i) {
            // dpm1 = m1.data + _i*_c;
            dpm2 = m2.data;
            for (_j = 0; _j < _m; ++_j) {
                m1_inner = dpm1;
                m2_inner = dpm2;
                for (_k = 0; _k < _c; ++_k) {
                    *dptr += (*m1_inner++)*(*m2_inner);
                    m2_inner += _m;
                }
                ++dptr;
                ++dpm2; // this gets `dpm2` pointing to the first element of the correct column
            }
            dpm1 += _c;
        }
        uint64_t *nshape = new uint64_t[2];
        *nshape = _n;
        *(nshape + 1) = _m;
        return {pdata, m1._shape._r, nshape, _volume, false}; // ctor WON'T make a copy of shape array
    }
    namespace tens_ops {
        template <Numeric U, Numeric V>
        matrix<mulComT<U, V>> hadamard(const matrix<U> &mat1, const matrix<V> &mat2) {
            if (mat1._shape._r == -1) {
                if (mat2._shape._r != -1)
                    throw exceptions::empty_tensor_error{"Error: addition between an empty and a non-empty tensor "
                                                         "cannot be performed.\n"};
                return {};
            }
            if (mat2._shape._r == -1)
                throw exceptions::empty_tensor_error{"Error: addition between an empty and a non-empty tensor "
                                                     "cannot be performed.\n"};
            if (*mat1._shape._s != *mat2._shape._s || *(mat1._shape._s + 1) != *(mat2._shape._s + 1))
                throw exceptions::dimension_mismatch_error{"Error: the Hadamard product can only be performed between "
                                                           "two matrices of equal shapes.\n"};
            mulComT<U, V> *_ndata = new mulComT<U, V>[mat1.vol];
            uint64_t counter = mat1.vol;
            mulComT<U, V> *nptr = _ndata;
            U *d1 = mat1.data;
            V *d2 = mat2.data;
            while (counter --> 0)
                *nptr++ = (*d1++)*(*d2++);
            return {_ndata, mat1._shape._r, mat1._shape._s, mat1.vol, true};
        }
    }
}
#endif
