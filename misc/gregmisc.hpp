#ifndef GREGMISC_HPP
#define GREGMISC_HPP

#include <vector>
#include <string>
#include <dirent.h>
#include <stdexcept>
#include <ios>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <cmath>

/* C++ header file containing miscellaneous function definitions. */

#include <iostream>

#ifdef __CUDACC__
#warning "Compiling with nvcc.\n"
#include "cuda_runtime.h"
#define HOST_DEVICE __host__ __device__
#define DEVICE __device__
#define GLOBAL __global__
#include "cuda_runtime.h"
#define PI 3.14159265358979323846264338327950288419716939937510582097494459
#define ROOT2 1.41421356237309504880
#define LIGHT_SPEED 299'792'458.0 // m/s
#define PERMITTIVITY 0.0000000000088541878188 // F/m
#define PERMEABILITY 0.00000125663706127 // F/m
#else
#define HOST_DEVICE
#define DEVICE
#define GLOBAL
#define PI 3.14159265358979323846264338327950288419716939937510582097494459l
#define ROOT2 1.41421356237309504880l
#define LIGHT_SPEED 299'792'458.0l // m/s
#define PERMITTIVITY 0.0000000000088541878188l // F/m
#define PERMEABILITY 0.00000125663706127 // F/m
#endif

namespace gtd {
    template <typename T>
    concept numeric = requires (T a, T b) {
        {a + b} -> std::convertible_to<T>;
        {a + b} -> std::convertible_to<T>;
        {a * b} -> std::convertible_to<T>;
        {a / b} -> std::convertible_to<T>;
    };
    template <typename C>
    concept callable = requires (C a) {
        {a()};
    };
    template <typename C, typename T, typename R = T>
    concept callret = requires (C a, T t) {
        {a(t)} -> std::same_as<R>;
    };
    template <typename C, typename T, typename R = T>
    concept calldblret = requires (C a, T t1, T t2) {
        {a(t1, t2)} -> std::same_as<R>;
    };
    HOST_DEVICE bool str_eq(const char *s1, const char *s2) {
        if (!s1 || !s2)
            return false;
        while (*s1 || *s2)
            if (*s1++ != *s2++)
                return false;
        return true;
    }
#ifndef GREGALG_HPP
    HOST_DEVICE bool memcopy(void *dst, const void *src, size_t bytes) {
        if (!bytes || !dst || !src)
            return false;
        char *cdst = (char*) dst;
        const char *csrc = (char*) src;
        while (bytes --> 0)
            *cdst++ = *csrc++;
        return true;
    }
#endif
#ifndef GREGSTR_HPP
    HOST_DEVICE size_t strlen_c(const char *str) {
        if (!str || !*str)
            return -1;
        size_t len = 0;
        while (*str++) ++len;
        return len;
    }
    HOST_DEVICE char *strcpy_c(char *dst, const char *src) {
        if (!dst || !src)
            return nullptr;
        char *org = dst;
        while (*src)
            *dst++ = *src++;
        *dst = 0;
        return org;
    }
    HOST_DEVICE char *strcat_c(char *dst, const char *src) {
        if (!dst || !src)
            return nullptr;
        char *org = dst;
        while (*dst) ++dst;
        while (*src)
            *dst++ = *src++;
        *dst = 0;
        return org;
    }
    HOST_DEVICE bool endswith(const char *str, const char *with) {
        if (!str || !with)
            return false;
        uint64_t _slen = strlen_c(str);
        uint64_t _wlen = strlen_c(with);
        if (_wlen > _slen)
            return false;
        str += _slen - _wlen;
        while (*str)
            if (*str++ != *with++)
                return false;
        return true;
    }
#endif
    template <std::integral T>
    void to_string(T val, char *buffer) {
        if (!val) {
            *buffer++ = 48;
            *buffer = 0;
            return;
        }
        if (val < 0) {
            *buffer++ = '-';
            val *= -1;
        }
        uint64_t i = log10l(val);
        while (i --> 0) ++buffer;
        *(buffer + 1) = 0;
        while (val > 0) {
            *buffer-- = (val % 10) + 48;
            val /= 10;
        }
    }
    template <std::unsigned_integral T>
    bool to_uint(const char *str, T *val) { // leaves val unchanged in case of first error (obviously), -1 in second er.
        if (!str || !*str || !val)
            return false;
        T result = 0;
        goto start;
        do {
            result *= 10;
            start:
            if (*str < 48 || *str > 57) {
                *val = (T) -1;
                return false;
            }
            result += *str++ - 48;
        } while (*str);
        *val = result;
        return true;
    }
    template <std::floating_point T>
    bool to_float(const char *str, T *out) {
        if (!str || !*str || !out)
            return false;
        char *endptr{};
        T val;
        if constexpr (std::same_as<T, float>)
            val = strtof(str, &endptr);
        if constexpr (std::same_as<T, double>)
            val = strtod(str, &endptr);
        if constexpr (std::same_as<T, long double>)
            val = strtold(str, &endptr);
        if (endptr == str)
            return false;
        *out = val;
        return true;
    }
    uint64_t write_all(int fd, const void *buff, size_t count) {
        size_t bwritten;
        size_t rem = count;
        const char *rembuff = (const char *) buff;
        uint64_t tot_written = 0;
        while (rem) {
            if ((bwritten = write(fd, rembuff, rem)) == -1)
                return tot_written;
            tot_written += bwritten;
            rem -= bwritten;
            rembuff += bwritten;
        }
        return tot_written;
    }
    uint64_t read_all(int fd, void *buff, size_t count = -1) {
        size_t bread;
        size_t rem = count;
        char *rembuff = (char *) buff;
        uint64_t tot_read = 0;
        if (count == ((size_t) -1)) {
            while ((bread = read(fd, rembuff, rem)) > 0) {
                tot_read += bread;
                rem -= bread;
                rembuff += bread;
            }
        } else {
            while (tot_read != count) {
                if ((bread = read(fd, rembuff, rem)) <= 0)
                    return tot_read;
                tot_read += bread;
                rem -= bread;
                rembuff += bread;
            }
        }
        return tot_read;
    }
    [[nodiscard]] std::vector<std::string> *find_files(const char *dirpath,
                                                       const char *extension,
                                                       bool prep_dpath) {
        if (!dirpath || !extension)
            throw std::invalid_argument{"Error: path to directory and extension cannot be nullptr.\n"};
        DIR *dir;
        if (!(dir = opendir(dirpath))) {
            std::string error = "Error: could not open directory \"";
            error += dirpath;
            error += "\".\n";
            throw std::ios_base::failure{error};
        }
        struct dirent *entry;
        std::vector<std::string> *files = new std::vector<std::string>{};
        if (prep_dpath) {
            if (*(dirpath + strlen_c(dirpath) - 1) == '/') {
                while ((entry = readdir(dir)))
                    if (gtd::endswith(entry->d_name, extension)) {
                        files->emplace_back(dirpath);
                        files->back() += entry->d_name;
                    }
            } else
                while ((entry = readdir(dir)))
                    if (gtd::endswith(entry->d_name, extension)) {
                        files->emplace_back(dirpath);
                        files->back().push_back('/');
                        files->back() += entry->d_name;
                    }
        } else
            while ((entry = readdir(dir)))
                if (gtd::endswith(entry->d_name, extension))
                    files->emplace_back(entry->d_name);
        closedir(dir);
        if (files->empty()) {
            delete files;
            files = nullptr;
        }
        return files;
    }
    std::vector<std::string> *find_files(const char *dirpath,
                                         const char *extension,
                                         std::vector<std::string> *ptr = nullptr) {
        /* Recursive version of the above function. Finds all files within a given directory and all its subdirectories
         * with the given extension. Always prepends the dir. path, given the recursive traversal of the dir. tree. */
        /* Not thread-safe due to the use of `strerror`, and could also cause problems due to changing the CWD. */
        if (!dirpath)
            throw std::invalid_argument{"Error: path to directory cannot be nullptr.\n"};
        if (!*dirpath)
            throw std::invalid_argument{"Error: path to directory cannot be an empty string.\n"};
        char cwd[PATH_MAX];
        if (!getcwd(cwd, PATH_MAX)) {
            std::string error = "Error: could not obtain current working directory.\nReason: ";
            error += strerror(errno);
            error.push_back('\n');
            throw std::ios_base::failure{error};
        }
        std::unique_ptr<char[]> _new_dpath{};
        if (*dirpath == '.' && !*(dirpath + 1))
            dirpath = cwd;
        else {
            _new_dpath.reset(new char[PATH_MAX]);
            if (!realpath(dirpath, _new_dpath.get())) {
                _new_dpath.reset();
                std::string error = "Error: could not obtain the fully-expanded path for \"";
                error += dirpath;
                error += "\".\nReason: ";
                error += strerror(errno);
                error.push_back('\n');
                throw std::ios_base::failure{error};
            }
            dirpath = _new_dpath.get();
        }
        if (chdir(dirpath) == -1) {
            std::string error = "Error: the current working directory could not be changed to \"";
            error += dirpath;
            error += "\".\nCurrent working directory: \"";
            error += cwd;
            error += "\"\nReason: ";
            error += strerror(errno);
            error.push_back('\n');
            throw std::ios_base::failure{error};
        }
        DIR *dir;
        if (!(dir = opendir("." /* dirpath */))) {
            std::string error = "Error: could not open directory \"";
            error += dirpath;
            error += "\".\nCurrent working directory: \"";
            error += cwd;
            error += "\"\nReason: ";
            error += strerror(errno);
            error.push_back('\n');
            throw std::ios_base::failure{error};
        }
        struct dirent *entry{};
        struct stat buff{};
        if (!ptr)
            ptr = new std::vector<std::string>{};
        uint64_t _dlen;
        std::unique_ptr<char[]> _sub_dpath{};
        if (*(dirpath + (_dlen = strlen_c(dirpath)) - 1) == '/') {
            while ((entry = readdir(dir))) {
                if (str_eq(entry->d_name, ".") || str_eq(entry->d_name, ".."))
                    continue;
                if (stat(entry->d_name, &buff) == -1) // ignore
                    continue;
                if (S_ISDIR(buff.st_mode)) {
                    _sub_dpath.reset(new char[_dlen + strlen_c(entry->d_name) + 1]);
                    strcpy_c(_sub_dpath.get(), dirpath);
                    strcpy_c(_sub_dpath.get() + _dlen, entry->d_name);
                    find_files(_sub_dpath.get(), extension, ptr);
                    _sub_dpath.reset();
                    continue;
                }
                if (endswith(entry->d_name, extension)) {
                    ptr->emplace_back(dirpath);
                    ptr->back() += entry->d_name;
                }
            }
        } else {
            while ((entry = readdir(dir))) {
                if (str_eq(entry->d_name, ".") || str_eq(entry->d_name, ".."))
                    continue;
                if (stat(entry->d_name, &buff) == -1)
                    continue;
                if (S_ISDIR(buff.st_mode)) {
                    _sub_dpath.reset(new char[_dlen + strlen_c(entry->d_name) + 2]);
                    strcpy_c(_sub_dpath.get(), dirpath);
                    *(_sub_dpath.get() + _dlen) = '/';
                    strcpy_c(_sub_dpath.get() + _dlen + 1, entry->d_name);
                    find_files(_sub_dpath.get(), extension, ptr);
                    _sub_dpath.reset();
                    continue;
                }
                if (endswith(entry->d_name, extension)) {
                    ptr->emplace_back(dirpath);
                    ptr->back().push_back('/');
                    ptr->back() += entry->d_name;
                }
            }
        }
        closedir(dir);
        _new_dpath.reset();
        if (chdir(cwd) == -1) {
            std::string error = "Error: the current working directory could not be changed back to \"";
            error += cwd;
            error += "\".\nReason: ";
            error += strerror(errno);
            error.push_back('\n');
            throw std::ios_base::failure{error};
        }
        return ptr;
    }
    template <numeric T>
    HOST_DEVICE T rad_to_deg(const T &rad) {
        return rad*180.0l/PI;
    }
    template <numeric T>
    HOST_DEVICE T deg_to_rad(const T &deg) {
        return deg*PI/180.0l;
    }
    template <numeric T>
    HOST_DEVICE T abs(const T &x) {
        return x >= 0 ? x : -x;
    }
    template <typename T>
    void print_array(const T *arr, uint64_t size) {
        if (!arr) {
            printf("[(null)]\n");
            return;
        }
        if (!size) {
            printf("[]\n");
            return;
        }
        putchar('[');
        while (size --> 0)
            std::cout << *arr++ << ", ";
        printf("\b\b]\n");
    }
    template <numeric T, uint64_t N>
    class point final { // point in N-dimensional space
        T data[N];
    public:
        template <typename ...Args> requires (std::same_as<Args, T> && ...)
        HOST_DEVICE explicit point(Args&& ...args) : data{args...} {
            static_assert(sizeof...(args) == N, "Number of constructor arguments does not match N.\n");
        }
        HOST_DEVICE explicit point(const T *_data) {
            if (!_data)
                throw std::invalid_argument{"Error: nullptr passed as pointer to data.\n"};
            memcopy(this->data, _data, sizeof(T)*N);
        }
        HOST_DEVICE point(const std::initializer_list<T> &_li) {
            if (_li.size() != N)
                throw std::invalid_argument{"Error: initialiser list size does not match N.\n"};
            memcopy(this->data, _li.begin(), sizeof(T)*N);
        }
        HOST_DEVICE T mag_sq() const noexcept {
            if constexpr (N == 0)
                return 0;
            if constexpr (N == 1)
                return data[0]*data[0];
            if constexpr (N == 2)
                return data[0]*data[0] + data[1]*data[1];
            if constexpr (N == 3)
                return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
            if constexpr (N == 4)
                return data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3];
            T *ptr = data;
            T _sum = 0;
            uint64_t counter = N;
            while (counter --> 0) {
                _sum += (*ptr)*(*ptr);
                ++ptr;
            }
            return _sum;
        }
        HOST_DEVICE T mag() const noexcept {
            return std::sqrt(this->mag_sq());
        }
        template <numeric U, uint64_t M>
        HOST_DEVICE friend U distance(const point<U, M> &p1, const point<U, M> &p2);
        template <numeric U, uint64_t M>
        HOST_DEVICE friend U max_dist(const point<U, M> *points, uint64_t _num);
    };
    template <numeric U, uint64_t M>
    HOST_DEVICE U distance(const point<U, M> &p1, const point<U, M> &p2) {
        const U *ptr1 = p1.data;
        const U *ptr2 = p2.data;
        U _sum = 0;
        U _diff;
        uint64_t counter = M;
        while (counter --> 0) {
            _diff = *ptr1++ - *ptr2++;
            _sum += _diff*_diff;
        }
        return std::sqrt(_sum);
    }
    template <numeric U, uint64_t M>
    HOST_DEVICE U max_dist(const point<U, M> *points, uint64_t _num) {
        if constexpr (M == 0 || M == 1)
            return 0;
        if (!points || !_num)
            return 0;
        static constexpr uint64_t Mm1 = M - 1;
        U _maxd = 0;
        U _d;
        const point<U, M> *outer = points;
        const point<U, M> *inner;
        uint64_t _i = 0;
        uint64_t _j;
        while (_i++ < Mm1) {
            inner = outer + 1;
            _j = _i;
            while (_j++ < M)
                if ((_d = distance(*outer, *inner++)) > _maxd)
                    _maxd = _d;
            ++outer;
        }
        return _maxd;
    }
#ifdef __CUDACC__
    class cuda_error : public std::exception {
        char *msg{};
    public:
        explicit cuda_error(int line, const char *file, const char *error_s) {
            msg = new char[40 + gtd::strlen_c(file) + gtd::strlen_c(error_s)];
            gtd::strcpy_c(msg, "Error \"");
            gtd::strcat_c(msg, error_s);
            gtd::strcat_c(msg, "\" occurred in file \"");
            gtd::strcat_c(msg, file);
            gtd::strcat_c(msg, "\" on line ");
            gtd::to_string(line, msg);
            gtd::strcat_c(msg, ".\n");
        }
        const char *what() const noexcept override {
            return this->msg;
        }
        ~cuda_error() {
            delete [] msg;
        }
    };
#define CUDA_ERROR(func_call) \
    { cudaError_t err; \
    if ((err = func_call) != cudaSuccess) { \
        fprintf(stderr, "Error: %s\n", cudaGetErrorString(err)); \
        throw gtd::cuda_error{__LINE__, __FILE__, cudaGetErrorString(err)}; \
    } }
#endif
    const char *get_time() {
        static char _str[26]{};
        time_t _now = time(nullptr);
        if (_now == ((time_t) -1))
            return nullptr;
        char *_s = ctime(&_now);
        if (!_s)
            return nullptr;
        strcpy_c(_str, _s);
        _str[3] = '_';
        _str[7] = '_';
        if (_str[8] == ' ') // in case the day of the month is < 10
            _str[8] = '0'; // just found this error as it is the 1st of August today lol
        _str[10] = '_';
        _str[13] = 'h';
        _str[16] = 'm';
        _str[19] = 's';
        _str[24] = _str[23];
        _str[23] = _str[22];
        _str[22] = _str[21];
        _str[21] = _str[20];
        _str[20] = '_';
        return _str;
    }
}
#endif
