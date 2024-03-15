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

/* C++ header file containing miscellaneous function definitions. */

#include <iostream>

namespace gtd {
    bool str_eq(const char *s1, const char *s2) {
        if (!s1 || !s2)
            return false;
        while (*s1 || *s2)
            if (*s1++ != *s2++)
                return false;
        return true;
    }
#ifndef GREGALG_HPP
    bool memcopy(void *dst, const void *src, size_t bytes) {
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
    size_t strlen_c(const char *str) {
        if (!str || !*str)
            return -1;
        size_t len = 0;
        while (*str++) ++len;
        return len;
    }
    char *strcpy_c(char *dst, const char *src) {
        if (!dst || !src)
            return nullptr;
        char *org = dst;
        while (*src)
            *dst++ = *src++;
        *dst = 0;
        return org;
    }
    bool endswith(const char *str, const char *with) {
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
        if (!dirpath)
            throw std::invalid_argument{"Error: path to directory cannot be nullptr.\n"};
        DIR *dir;
        if (!(dir = opendir(dirpath))) {
            std::string error = "Error: could not open directory \"";
            error += dirpath;
            error += "\".\n";
            throw std::ios_base::failure{error};
        }
        char cwd[PATH_MAX];
        getcwd(cwd, PATH_MAX);
        chdir(dirpath);
        struct dirent *entry{};
        struct stat buff{};
        if (!ptr)
            ptr = new std::vector<std::string>{};
        if (*(dirpath + strlen_c(dirpath) - 1) == '/') {
            while ((entry = readdir(dir))) {
                if (str_eq(entry->d_name, ".") || str_eq(entry->d_name, ".."))
                    continue;
                if (stat(entry->d_name, &buff) == -1) // ignore
                    continue;
                if (S_ISDIR(buff.st_mode)) {
                    find_files(entry->d_name, extension, ptr);
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
                    find_files(entry->d_name, extension, ptr);
                    continue;
                }
                if (endswith(entry->d_name, extension)) {
                    ptr->emplace_back(dirpath);
                    ptr->back().push_back('/');
                    ptr->back() += entry->d_name;
                }
            }
        }
        chdir(cwd);
        return ptr;
    }
}
#endif
