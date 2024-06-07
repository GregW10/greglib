#ifndef GREGCOMPLEX_HPP
#define GREGCOMPLEX_HPP

// Nothing special, just a GPU-friendly complex number class

#include "gregmisc.hpp"

#ifndef GREGSTACK_HPP
#ifdef __CUDACC__
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif
#endif

namespace gtd {
    class invalid_complex_operation : public std::logic_error {
    public:
        invalid_complex_operation() : std::logic_error{"Error: this operation is invalid for a complex number.\n"} {}
        invalid_complex_operation(const char *msg) : std::logic_error{msg} {}
    };
    template <numeric T>
    class complex {
        T _real{};
	T _imag{};
    public:
	HOST_DEVICE complex() = default;
	HOST_DEVICE complex(const T &_real_, const T &_imag_) : _real{_real_}, _imag{_imag_} {}
        HOST_DEVICE complex(const complex<T> &other) : _real{other._real}, _imag{other._imag} {}
        HOST_DEVICE complex(complex<T> &&other) : _real{std::move(other._real)}, _imag{std::move(other._imag)} {} // moving is almost certainly overkill
	HOST_DEVICE T real() const noexcept {
            return this->_real;
        }
        HOST_DEVICE T imag() const noexcept {
            return this->_imag;
        }
        HOST_DEVICE void real(const T &nreal) noexcept {
            this->_real = nreal;
        }
        HOST_DEVICE void imag(const T &nimag) noexcept {
            this->_imag = nimag;
        }
	HOST_DEVICE T mag() const noexcept {
            return std::sqrt(this->_real*this->_real + this->_imag*this->_imag);
        }
        HOST_DEVICE T mag_sq() const noexcept {
            return this->_real*this->_real + this->_imag*this->_imag;
        }
        HOST_DEVICE T arg() const {
            if (this->_real == 0 && this->_imag == 0)
                throw invalid_complex_operation{"Error: the argument of the zero complex number is undefined.\n"};
            return std::acos(this->_real/this->mag_sq());
        }
        HOST_DEVICE friend complex<T> operator*(const T &s, const complex<T> &c) {
            return {s*c._real, s*c._imag};
        }
        HOST_DEVICE friend complex<T> operator*(const complex<T> &c, const T &s) {
            return {s*c._real, s*c._imag};
        }
        HOST_DEVICE friend complex<T> operator+(const complex<T> &c1, const complex<T> &c2) {
            return {c1._real + c2._real, c1._imag + c2._imag};
        }
        HOST_DEVICE friend complex<T> operator-(const complex<T> &c1, const complex<T> &c2) {
            return {c1._real - c2._real, c1._imag - c2._imag};
        }
        HOST_DEVICE friend complex<T> operator*(const complex<T> &c1, const complex<T> &c2) {
            return {c1._real*c2._real - c1._imag*c2._imag, c1._real*c2._imag + c1._imag*c2._real};
        }
    };
}
#endif

