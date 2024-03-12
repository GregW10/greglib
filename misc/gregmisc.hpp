#ifndef GREGMISC_HPP
#define GREGMISC_HPP

#include <vector>
#include <string>
#include <dirent.h>
#include "../nbod/gregstr.hpp"

/* C++ header file containing miscellaneous function definitions. */

namespace gtd {
    template <bool PrependDirpath = true>
    std::vector<std::string> find_files(const char *dirpath, const char *extension) {
        if (!dirpath || !extension)
            throw std::invalid_argument{"Error: path to directory and extension cannot be nullptr.\n"};
        DIR *dir;
        if (!(dir = opendir(dirpath)))
            throw std::ios_base::failure{"Error: could not open directory.\n"};
        struct dirent *entry;
        std::vector<std::string> *files = new std::vector<std::string>{};
        while ((entry = readdir(dir))) {
            if (gtd::endswith(entry->d_name, extension)) {
                if constexpr (PrependDirpath) {
                    files->emplace_back(dirpath);
                    files.back().push_back('/');
                    files.back() += entry->d_name;
                } else {
                    files->emplace_back(entry->d_name);
                }
            }
        }
        closedir(dir);
        if (files->empty()) {
            delete files;
            files = nullptr;
        }
        return files;
    }
}

#endif