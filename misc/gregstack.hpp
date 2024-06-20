#ifndef GREGSTACK_HPP
#define GREGSTACK_HPP

#include "gregmisc.hpp"
#include <iostream>
#include <cstdint>
#include <cinttypes>
#include <stdexcept>
#include <string>
#include <cassert>

#define DEF_STACK_SIZE 32

#ifndef GREGCOMPLEX_HPP
namespace diff {
    template <gtd::numeric T, gtd::numeric R, gtd::callret<T, R> F, bool = false>
    HOST_DEVICE R simpquad(const F&,
                           T,
                           T,
                           T = 0.0625l/1024.0l,
                           T = 0.0625l/1024.0l,
                           T = 1/(1024.0l*1024.0l),
                           uint64_t* = nullptr);
    template <gtd::numeric T, gtd::numeric R, gtd::calldblret<T, R> F, gtd::callret<T> GH, bool = false>
    HOST_DEVICE R simpdblquad(const F&,
                              T,
                              T,
                              const GH&,
                              const GH&,
                              T = 0.0625l/1024.0l,
                              T = 0.0625l/1024.0l,
                              T = 1/(1024.0l*1024.0l),
                              uint64_t* = nullptr,
                              T = 0.0625l/1024.0l,
                              T = 0.0625l/1024.0l,
                              T = 1/(1024.0l*1024.0l),
                              uint64_t* = nullptr);
}
#endif

namespace gtd {
#ifdef __CUDACC__
    class bad_gpu_alloc : std::runtime_error {
    public:
        bad_gpu_alloc() : std::runtime_error{"Error: could not allocate memory block on the GPU.\n"} {}
        explicit bad_gpu_alloc(const std::string &_s) : std::runtime_error{_s} {}
    };
    class bad_gpu_dealloc : std::runtime_error {
    public:
        bad_gpu_dealloc() : std::runtime_error{"Error: could not deallocate memory block on the GPU.\n"} {}
        explicit bad_gpu_dealloc(const std::string &_s) : std::runtime_error{_s} {}
    };
#endif
    class empty_stack_error : public std::logic_error {
    public:
        empty_stack_error() : std::logic_error{"Error: operation cannot be performed on an empty stack.\n"} {}
        explicit empty_stack_error(const std::string &_s) : std::logic_error{_s} {}
    };
    template <typename T>
    class stack {
        T *_data{}; // pointer to entire data block
        uint64_t _asize{}; // total size allocated for elements (max. possible num. of elements)
        uint64_t _size{}; // actual number of elements stored
        T *_top{}; // pointer to top of stack, one past the last element
    public:
#ifndef __CUDACC__
        stack() : _data{new T[DEF_STACK_SIZE]}, _asize{DEF_STACK_SIZE}, _top{_data} {}
        explicit stack(uint64_t size) : _data{size ? new T[size] : nullptr}, _asize{size}, _top{_data} {}
#else
        __host__ __device__ explicit stack(uint64_t size = DEF_STACK_SIZE) : _asize{size} {
            if (!size)
                return;
            // cudaError_t err;
            // if ((err = cudaMalloc(&_data, _asize*sizeof(T))) != cudaSuccess) {
            //     std::string error = "Error: GPU memory could not be allocated because \"";
            //     error += cudaGetErrorString(err);
            //     error += "\".\n";
            //     throw bad_gpu_alloc{error};
            // }
            assert(cudaMalloc(&_data, _asize*sizeof(T)) == cudaSuccess);
            _top = _data;
        }
#endif
        HOST_DEVICE T &push(const T &val = T{}) {
            if (_size == _asize) {
                if (!_asize) {
                    _data = new T[DEF_STACK_SIZE];
                    _asize = DEF_STACK_SIZE;
                    _size = 1;
                    _top = _data;
                    return *_top++ = val;
                }
                _asize *= 2;
#ifndef __CUDACC__
                T *ndata = new T[_asize];
#else
                T *ndata;
                // if (cudaMalloc(&ndata, _asize*sizeof(T)) != cudaSuccess)
                //     throw bad_gpu_alloc{};
                assert(cudaMalloc(&ndata, _asize*sizeof(T)) == cudaSuccess);
#endif
                T *nptr = ndata;
                T *dptr = _data;
                uint64_t _i = _size;
                while (_i --> 0)
                    *nptr++ = *dptr++;
#ifndef __CUDACC__
                delete [] _data;
#else
                // if (cudaFree(_data) != cudaSuccess)
                //     throw bad_gpu_dealloc{};
                assert(cudaFree(_data) == cudaSuccess);
#endif
                _data = ndata;
                _top = _data + _size++;
                return *_top++ = val;
            }
            ++_size;
            return *_top++ = val;
        }
        HOST_DEVICE void pop() noexcept { // performs no checking whatsoever, must be done with caution
            (--_top)->~T();
            --_size;
            // if (_top < _data) {
            //     fprintf(stderr, "You done fucked up.\n");
            //     abort();
            // }
            // if (!_size)
            //     _top = nullptr;
        }
        HOST_DEVICE T pop_r() {
#ifndef __CUDACC__
            if (!_size)
                throw empty_stack_error{"Error: pop_r cannot be called on an empty stack.\n"};
#else
            assert(_size);
#endif
            T val = *--_top;
            _top->~T();
            --_size;
            return val;
        }
        //HOST_DEVICE bool reserve(uint64_t _nsize) {
        //    if (_nsize < _asize)
        //        return false;
        //
        //    return true;
        // }
        HOST_DEVICE T &top() noexcept { // performs no checking
            return *(_top - 1);
        }
        HOST_DEVICE const T &top() const {
            return *(_top - 1);
        }
        HOST_DEVICE bool empty() const noexcept {
            return !_size;
        }
        HOST_DEVICE uint64_t size() const noexcept {
            return _size;
        }
        HOST_DEVICE uint64_t capacity() const noexcept {
            return _asize;
        }
        HOST_DEVICE explicit operator bool() const noexcept {
            return _size;
        }
        HOST_DEVICE ~stack() {
#ifndef __CUDACC__
            delete [] _data;
#else
            assert(cudaFree(_data) == cudaSuccess);
#endif
            _asize = 0;
            _size = 0;
            _data = nullptr;
            _top = nullptr;
        }
        template <gtd::numeric U, gtd::numeric R, gtd::callret<U, R> F, bool prog>
        friend HOST_DEVICE R diff::simpquad(const F&, U, U, U, U, U, uint64_t*);
        template <gtd::numeric U, gtd::numeric R, gtd::calldblret<U, R> F, gtd::callret<U> GH, bool prog>
        friend HOST_DEVICE R diff::simpdblquad(const F&, U, U, const GH&, const GH&,U,U,U,uint64_t*,U,U,U,uint64_t*);
    };
}

#endif
