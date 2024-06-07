#ifndef GREGCOMPLEX_HPP
#define GREGCOMPLEX_HPP

// Nothing special, just a GPU-friendly complex number class

#include "gregmisc.hpp"

namespace gtd {
    class invalid_complex_operation : public std::logic_error {
    public:
        invalid_complex_operation() : std::logic_error{"Error: this operation is invalid for a complex number.\n"} {}
        explicit invalid_complex_operation(const char *msg) : std::logic_error{msg} {}
    };
#ifndef __CUDACC__
    template <numeric T = long double>
#else
    template <numeric T = double>
#endif
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
	        if (this->_real == 0) {
	            if (this->_imag == 0)
	                throw invalid_complex_operation{"Error: the argument of the zero complex number is undefined.\n"};
	            return this->_imag > 0 ? this->_imag : -this->_imag;
	        }
	        if (this->_real > 0) {
	            if (this->_imag >= 0)
	                return std::atan(this->_imag/this->_real);
	            return 2*PI + std::atan(this->_imag/this->_real);
	        }
	        return PI + std::atan(this->_imag/this->_real);
            /* if (this->_real == 0 && this->_imag == 0)
                throw invalid_complex_operation{"Error: the argument of the zero complex number is undefined.\n"};
            return std::acos(this->_real/this->mag()); */
        }
    	HOST_DEVICE complex<T> &operator+=(const complex<T> &other) {
		    this->_real += other._real;
	    	this->_imag += other._imag;
	    	return *this;
	    }
    	HOST_DEVICE complex<T> &operator-=(const complex<T> &other) {
	    	this->_real -= other._real;
	    	this->_imag -= other._imag;
	    	return *this;
	    }
    	HOST_DEVICE complex<T> &operator*=(const complex<T> &other) {
	    	this->_real = this->_real*other._real - this->_imag*other._imag;
	    	this->_imag = this->_real*other._imag + this->_imag*other._real;
	    	return *this;
	    }
    	HOST_DEVICE complex<T> &operator*=(const T &scalar) {
	    	this->_real *= scalar;
	    	this->_imag *= scalar;
	    	return *this;
	    }
    	HOST_DEVICE complex<T> &operator+=(const T &scalar) { // scalar is interpreted as a complex number `scalar + 0i`
	    	this->_real += scalar;
	    	return *this;
	    }
        HOST_DEVICE friend complex<T> operator*(const T &s, const complex<T> &c) {
            return {s*c._real, s*c._imag};
        }
        HOST_DEVICE friend complex<T> operator*(const complex<T> &c, const T &s) {
            return {s*c._real, s*c._imag};
        }
    	HOST_DEVICE friend complex<T> operator+(const complex<T> &c, const T &s) {
	    	return {c._real + s, c._imag};
	    }
    	HOST_DEVICE friend complex<T> operator+(const T &s, const complex<T> &c) {
	    	return {c._real + s, c._imag};
	    }
        HOST_DEVICE friend complex<T> operator/(const complex<T> &c, const T &s) {
	        return {c._real/s, c._imag/s};
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
        HOST_DEVICE complex<T> &operator=(const complex<T> &other) {
	        if (&other == this)
	            return *this;
	        this->_real = other._real;
	        this->_imag = other._imag;
	        return *this;
	    }
        HOST_DEVICE complex<T> &operator=(complex<T> &&other) {
	        if (&other == this)
	            return *this;
	        this->_real = std::move(other._real);
	        this->_imag = std::move(other._imag);
	        return *this;
	    }
        friend std::ostream &operator<<(std::ostream &os, const complex<T> &c) {
	        os << +c._real;
	        if (c._imag >= 0)
	            return os << " + " << +c._imag << "i";
	        return os << " - " << -c._imag << "i";
	    }
    };
}
#endif
