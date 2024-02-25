#ifndef GREGMTX_HPP
#define GREGMTX_HPP

#include "gregtens.hpp"

namespace gml {
    template <Numeric T>
    class matrix : public tensor<T> {
        void load_mtsr(std::ifstream &in, size_t st_size) {
            uint64_t shape[2];
            in.read((char *) shape, 2*sizeof(uint64_t));
            if (!shape[0] || !shape[1]) { // case for matrix with zero size
                if (st_size > 4 + 3*sizeof(uint64_t))
                    throw exceptions::invalid_mtsr_format{"Error: invalid .mtsr format - file too large for empty "
                                                          "matrix.\n"};
                tensor<T>::_shape._r = 0;
                tensor<T>::_shape._s = new ull_t[0];
                return;
            } else {
                if (st_size <= 4 + 3*sizeof(uint64_t))
                    throw exceptions::invalid_mtsr_format{"Error: invalid .mtsr format - file too small for non-empty "
                                                          "matrix.\n"};
            }
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
            tensor<T>::_shape._s = new uint64_t[2]{shape[0], shape[1]};
            tensor<T>::_shape._r = 2;
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
    public:
        matrix() = default;
        explicit matrix(const char *tsr_path) : tensor<T>{false} {
            delegate_mtsr(tsr_path);
        }
        matrix(uint64_t rows, uint64_t columns) : tensor<T>{tensor_shape{rows, columns}} {} // matrix full of zeros
        matrix(const std::vector<T> &_data, uint64_t rows, uint64_t columns) :
        tensor<T>{_data, tensor_shape{rows, columns}} {}
        matrix(const T *_data, uint64_t rows, uint64_t columns) : tensor<T>{_data, tensor_shape{rows, columns}} {}
        matrix(const std::initializer_list<std::initializer_list<T>> &li) : tensor<T>{false} {
            typename std::initializer_list<T>::size_type li_size = li.size();
            tensor<T>::_shape._r = 2; // rank will always be 2 for a matrix
            tensor<T>::_shape._s = new ull_t[2]{}; // default initialised to zeros
            tensor<T>::_shape.sizes = new ull_t[1]{}; // zero
            if (!li_size) { // emtpy initializer list, so no action taken
                return; // does any other action have to be taken? check this
            }
            const std::initializer_list<T> *ptr = li.begin();
            typename std::initializer_list<T>::size_type sub_elems = ptr->size();
            if (!sub_elems) // list of empty initializer lists
                return;
            tensor<T>::_shape._s[0] = li_size;
            tensor<T>::_shape._s[1] = sub_elems;
            tensor<T>::vol = li_size*sub_elems;
            tensor<T>::data = new T[tensor<T>::vol];
            T *dptr = tensor<T>::data;
            std::cout << "li_size: " << li_size << "\nsub_elems: " << sub_elems << "\nvol: " << this->vol << '\n';
            while (li_size --> 0) {
                if (ptr->size() != sub_elems) {
                    delete [] tensor<T>::data; // is this safe? - because of stack unwinding in exception throwing
                    throw std::invalid_argument{"Error: sizes of nested initializer lists do not match.\n"};
                }
                // memcopy(dptr++, ptr++->begin(), sizeof(T), sub_elems);
                gen::memcopy(dptr, ptr++->begin(), sizeof(T), sub_elems);
                dptr += sub_elems;
            }
            tensor<T>::_shape.sizes[0] = sub_elems;
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
        T &operator()(uint64_t _i, uint64_t _j) noexcept { // doesn't throw, so must be careful when called
            return *(tensor<T>::data + _i*(*(tensor<T>::_shape._s + 1)) + _j);
        }
        const T &operator()(uint64_t _i, uint64_t _j) const noexcept {
            return *(tensor<T>::data + _i*(*(tensor<T>::_shape._s + 1)) + _j);
        }
    };
}
#endif
