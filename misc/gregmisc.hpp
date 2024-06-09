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
    HOST_DEVICE void to_string(T val, char *buffer) {
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
}
#endif
