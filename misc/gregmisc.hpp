#ifndef GREGMISC_HPP
#define GREGMISC_HPP

#include <vector>
#include <string>
#include <dirent.h>
#include "../nbod/gregstr.hpp"

/* C++ header file containing miscellaneous function definitions. */

namespace gtd {
    std::vector<std::string> *find_files(const char *dirpath, const char *extension, bool prep_dpath = true) {
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
                        files.back() += entry->d_name;
                    }
            } else
                while ((entry = readdir(dir)))
                    if (gtd::endswith(entry->d_name, extension)) {
                        files->emplace_back(dirpath);
                        files.back().push_back('/');
                        files.back() += entry->d_name;
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
}
#endif
