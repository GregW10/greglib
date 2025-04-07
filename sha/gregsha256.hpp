#ifndef GREGSHA256_HPP
#define GREGSHA256_HPP

#include "gregshabase.hpp"

namespace gtd {
    class sha256 { // : public shasum {
        uint64_t _size = 0;
        bool _done = false;
        uint32_t _state[8] = {
            
        };
    public:
        sha256() = default;
    };
}

#endif

