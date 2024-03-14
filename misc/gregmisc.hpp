#ifndef GREGMISC_HPP
#define GREGMISC_HPP

#include <vector>
#include <string>
#include <dirent.h>
#include "../nbod/gregstr.hpp"

/* C++ header file containing miscellaneous function definitions. */

namespace gtd {
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
        return ptr;
    }
}
#endif
