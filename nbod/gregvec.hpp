//
// Created by mario on 01/09/2022.
//

#ifndef GREGVEC_H
#define GREGVEC_H

#ifndef __cplusplus
#error "The gregvec.hpp header file is a C++ header file only.\n"
#endif

#include "gregstr.hpp"
#include "gregmat.hpp"

#include <iostream>
#include <type_traits>
#include <utility>

namespace gtd { // forward declarations, to be able to use the functions inside the classes
    template <typename U, typename V>
    concept bitwiseOperands = requires (U val1, V val2) {
        {val1 << 1}; // must have all of these operators overloaded
        {val1 >> 1};
        {val1 <<= val2};
        {val1 >>= val2};
        {val1 | val2};
        {val1 & val2};
        {val1 ^ val2};
        {~val1};
        {val1 |= val2};
        {val1 &= val2};
        {val1 ^= val2};
        {val2 <<= val1};
        {val2 >>= val1};
        {val2 | val1};
        {val2 & val1};
        {val2 ^ val1};
        {~val2};
        {val2 |= val1};
        {val2 &= val1};
        {val2 ^= val1};
    };
    template <typename T>
    concept isIntegralNumWrapper = isNumWrapper<T> && requires (T val, T val2, size_t l, long double f) {
        {val % val2} -> isConvertible<T>; // modulo operator
        {val % l} -> isConvertible<T>;
        {l % val} -> isConvertible<T>;
        {val %= val2} -> isConvertible<T>;
        {val %= l} -> isConvertible<T>;
        {val % f} -> isConvertible<T>;
        {f % val} -> isConvertible<T>;
        {val %= f} -> isConvertible<T>;
    };
    namespace vec_ops {
        template <isNumWrapper U, isNumWrapper V>
        auto distance(const vector2D<U>& vec1, const vector2D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        auto distance(const vector2D<U>& vec1, const vector3D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        auto distance(const vector3D<U>& vec1, const vector2D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        auto distance(const vector3D<U>& vec1, const vector3D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        auto cross(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        auto cross(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        auto cross(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        auto cross(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V,long double>
        long double angle_between(const vector2D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V,long double>
        long double angle_between(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V,long double>
        long double angle_between(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V,long double>
        long double angle_between(const vector3D<U> &vec1, const vector3D<V> &vec2);
    }
    class division_by_zero : public std::invalid_argument {
    public:
        division_by_zero() : std::invalid_argument("Cannot divide by zero.") {}
        division_by_zero(const char *message) : std::invalid_argument(message) {}
    };
    template <isNumWrapper T>
    class vector {
    public:
        constexpr vector() noexcept = default;
        constexpr vector(const vector<T> &other) = default;
        constexpr vector(vector<T> &&other) noexcept = default;
        // template <isNumWrapper U> requires isConvertible<U, T>
        // constexpr vector(const vector<U> &other) noexcept {}
        // template <isNumWrapper U> requires isConvertible<U, T>
        // constexpr vector(const vector<U> &&other) noexcept {}
        virtual String str(unsigned char f_p_dec_places = 5) const = 0;
        virtual long double magnitude() const noexcept = 0;
        virtual T &operator[](unsigned char index) = 0;
        virtual const T &operator[](unsigned char index) const = 0;
        virtual vector<T> &operator++() noexcept = 0; //can only declare reference-returning func. for an abstract class
        virtual vector<T> &operator--() noexcept = 0;
        virtual vector<T> &operator=(const vector<T> &other) noexcept = 0; // copies other
        virtual vector<T> &operator=(vector<T> &&other) noexcept = 0; // copies other
        virtual ~vector() = default;
        template <isNumWrapper U>
        friend class vector;
    };
    template <isNumWrapper T>
    class vector3D;
    template <isNumWrapper T>
    class vector2D : public vector<T> {
    public:
        static inline const vector2D<T> zero{T{}, T{}};
        static inline const vector2D<T> left{T{-1}, T{}};
        static inline const vector2D<T> right{T{1}, T{}};
        static inline const vector2D<T> forwards{T{}, T{1}};
        static inline const vector2D<T> backwards{T{}, T{-1}};
    protected:
        T x{0};
        T y{0};
    public:
        constexpr vector2D() noexcept = default;
        constexpr vector2D(const vector2D<T> &other) : x{other.x}, y{other.y} {}
        constexpr vector2D(vector2D<T> &&other) noexcept : x{std::move(other.x)}, y{std::move(other.y)} {}
        template <isConvertible<T> U>
        constexpr vector2D(const vector2D<U> &other) noexcept : x{(T) other.x}, y{(T) other.y} {}
        template <isConvertible<T> U>
        constexpr vector2D(const vector2D<U> &&other) noexcept : x{(T) other.x}, y{(T) other.y} {}
        constexpr vector2D(const T &x_component, const T &y_component) noexcept : x{x_component}, y{y_component} {}
        constexpr vector2D(T &&x_component, T &&y_component) noexcept :
        x{std::move(x_component)}, y{std::move(y_component)} {}
        String str(unsigned char f_p_dec_places = 5) const override {
            String s;
            s.append_back(x, f_p_dec_places);
            if (y < T{0}) {
                s.append_back("i - ").append_back(-y, f_p_dec_places).append_back("j");
            }
            else {
                s.append_back("i + ").append_back(y, f_p_dec_places).append_back("j");
            }
            return s;
        }
        virtual vector2D<T> &set_x(T value) noexcept {
            x = value;
            return *this;
        }
        virtual vector2D<T> &set_y(T value) noexcept {
            y = value;
            return *this;
        }
        virtual vector2D<T> &assign(const T& xcomp, const T& ycomp) {
            this->x = xcomp;
            this->y = ycomp;
            return *this;
        }
        virtual vector2D<T> &assign(T&& xcomp, T&& ycomp) noexcept {
            this->x = std::move(xcomp);
            this->y = std::move(ycomp);
            return *this;
        }
        const T &get_x() const noexcept {
            return x;
        }
        const T &get_y() const noexcept {
            return y;
        }
        virtual vector2D<T> &make_zero() {
            this->x = T{0};
            this->y = T{0};
            return *this;
        }
        virtual vector2D<T> &set_length(const T &new_length) noexcept {
            long double frac = new_length/this->magnitude();
            this->x *= frac;
            this->y *= frac;
            return *this;
        }
        virtual vector2D<T> &set_length(const T &&new_length) noexcept {
            long double frac = new_length/this->magnitude();
            this->x *= frac;
            this->y *= frac;
            return *this;
        }
        std::pair<T, T> to_pair() const {
            return {this->x, this->y};
        }
        long double magnitude() const noexcept override {
            return sqrtl(static_cast<long double>(x)*static_cast<long double>(x) + // best to avoid call to
                         static_cast<long double>(y)*static_cast<long double>(y)); // std::pow() where possible
        }
        vector2D<long double> unit_vector() const noexcept {
            if (this->is_zero()) // no place for a division_by_zero exception if x and y are zero themselves
                return {};
            return *this / this->magnitude();
        }
        virtual vector2D<T> &normalise() noexcept {
            if (this->is_zero())
                return *this;
            long double mag = this->magnitude();
            this->x /= mag;
            this->y /= mag;
            return *this;
        }
        vector2D<T> x_projection() const {
            return {this->x, T{0}};
        }
        vector2D<T> y_projection() const {
            return {T{0}, this->y};
        }
        template <isNumWrapper U, isNumWrapper V>
        bool between(const vector2D<U> &v1, const vector2D<V> &v2) const noexcept {
            /* Tests if the vector lies within the 2D square spanned by v1 & v2 */
            return v1.x <= this->x && this->x < v2.x &&
                   v1.y <= this->y && this->y < v2.y;
        }
        virtual vector2D<T> &rotate(const long double &&angle_in_rad = _PI_) {
            return this->apply(matrix<long double>::get_2D_rotation_matrix(angle_in_rad));
        }
        virtual vector2D<T> &rotate(const long double &angle_in_rad = PI) {
            return this->apply(matrix<long double>::get_2D_rotation_matrix(angle_in_rad));
        }
        // template <isNumWrapper U = T>
        virtual vector2D<T> &rotate_to(const vector2D<T> &new_direction) noexcept {
            if (this->is_zero())
                return *this;
            return this->rotate(vec_ops::angle_between(*this, new_direction));
        }
        // template <isNumWrapper U = T>
        virtual vector2D<T> &rotate_to(const vector2D<T> &&new_direction) noexcept {
            return this->rotate_to(new_direction);
        }
        virtual bool is_zero() const noexcept {
            return this->x == T{0} && this->y == T{0};
        }
        virtual vector2D<T> &apply(const matrix<T> &transform) {
            if (transform.mat.size() != 2 || transform.mat[0].size() != 2)
                throw invalid_matrix_format("Only 2x2 matrices can be applied to a vector2D object.");
            T org_x = this->x;
            T org_y = this->y;
            this->x = transform.mat[0][0]*org_x + transform.mat[0][1]*org_y;
            this->y = transform.mat[1][0]*org_x + transform.mat[1][1]*org_y;
            return *this;
        }
        virtual vector2D<T> &apply(const matrix<T> &&transform) {
            return this->apply(transform);
        }
        vector2D<T> copy() const {
            return {this->x, this->y};
        }
        T &operator[](unsigned char index) override {
            switch (index) {
                case 0:
                    return this->x;
                case 1:
                    return this->y;
                default:
                    throw std::invalid_argument("Only the indices '0' and '1' are possible.\n");
            }
        }
        const T &operator[](unsigned char index) const override {
            switch (index) { // must repeat above code, as a non-const method cannot be called from a const method
                case 0:
                    return this->x;
                case 1:
                    return this->y;
                default:
                    throw std::invalid_argument("Only the indices '0' and '1' are possible.\n");
            }
        }
        vector2D<T> operator-() const {
            return {-this->x, -this->y};
        }
        virtual vector2D<T> &operator+=(const vector2D<T> &other)noexcept {//cannot requires-constrain virtual functions
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        virtual vector2D<T> &operator-=(const vector2D<T> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        virtual vector2D<T> &operator*=(const vector2D<T> &other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        virtual vector2D<T> &operator/=(const vector2D<T> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        vector2D<T> &operator%=(const vector2D<T> &other) requires isIntegralNumWrapper<T> { // cannot be virtual
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        virtual vector2D<T> &operator+=(const vector2D<T> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        virtual vector2D<T> &operator-=(const vector2D<T> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        virtual vector2D<T> &operator*=(const vector2D<T> &&other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        virtual vector2D<T> &operator/=(const vector2D<T> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        vector2D<T> &operator%=(const vector2D<T> &&other) requires isIntegralNumWrapper<T> {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(const vector2D<U> &other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(const vector2D<U> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(const vector2D<U> &other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator/=(const vector2D<U> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector2D<T> &operator%=(const vector2D<U> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator<<=(const vector2D<U> &other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator>>=(const vector2D<U> &other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator|=(const vector2D<U> &other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator&=(const vector2D<U> &other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator^=(const vector2D<U> &other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(const vector2D<U> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(const vector2D<U> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(const vector2D<U> &&other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator/=(const vector2D<U> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector2D<T> &operator%=(const vector2D<U> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator<<=(const vector2D<U> &&other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator>>=(const vector2D<U> &&other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator|=(const vector2D<U> &&other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator&=(const vector2D<U> &&other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator^=(const vector2D<U> &&other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        template <isConvertible<T> U>
        vector3D<T> &operator*=(const U &scalar) {
            this->x *= scalar;
            this->y *= scalar;
            return *this;
        }
        template <isConvertible<T> U>
        vector3D<T> &operator/=(const U &scalar) {
            this->x /= scalar;
            this->y /= scalar;
            return *this;
        }
        vector2D<T> operator++(int) noexcept { // cannot make virtual due to differing return types
            vector2D<T> retvec{this->x, this->y};
            ++this->x;
            ++this->y;
            return retvec;
        }
        vector2D<T> &operator++() noexcept override {
            ++this->x;
            ++this->y;
            return *this;
        }
        vector2D<T> operator--(int) noexcept { // cannot make virtual due to differing return types
            vector2D<T> retvec{this->x, this->y};
            --this->x;
            --this->y;
            return retvec;
        }
        vector2D<T> &operator--() noexcept override {
            --this->x;
            --this->y;
            return *this;
        }
        vector2D<T> operator~() requires requires (T val) {~val;} {
            return vector2D<T>(~this->x, ~this->y);
        }
        vector2D<T> &operator=(const vector<T> &other) noexcept override {
            if (&other == this)
                return *this;
            try {
                const vector2D<T> &oth = dynamic_cast<const vector2D<T>&>(other);
                this->x = oth.x;
                this->y = oth.y;
            } catch (const std::bad_cast &bc) {} // no action taken in case of std::bad_cast - vector obj. is unchanged
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator=(const vector2D<U> &other) noexcept {
            this->x = other.x;
            this->y = other.y;
            return *this;
        }
        vector2D<T> &operator=(const vector2D<T> &other) noexcept {
            this->x = other.x;
            this->y = other.y;
            return *this;
        }
        vector2D<T> &operator=(vector<T> &&other) noexcept override {
            if (&other == this)
                return *this;
            try {
                const vector2D<T> &oth = dynamic_cast<const vector2D<T>&>(other);
                this->x = std::move(oth.x);
                this->y = std::move(oth.y);
            } catch (const std::bad_cast &bc) {} // no action taken in case of std::bad_cast - vector obj. is unchanged
            return *this;
        }
        vector2D<T> &operator=(vector2D<T> &&other) noexcept {
            this->x = std::move(other.x); // in case T is an object and not a primitive type
            this->y = std::move(other.y);
            return *this;
        }
        virtual vector2D<T> &operator=(const T &value) {
            this->x = value;
            this->y = value;
            return *this;
        }
        virtual vector2D<T> &operator=(T &&value) noexcept {
            this->x = value;
            this->y = std::move(value);
            return *this;
        }
        template <isConvertible<T> U>
        vector2D<T> &operator=(const U &value) {
            this->x = value;
            this->y = value;
            return *this;
        }
        virtual operator bool() const noexcept { // to allow a vector object to be used as a boolean
            return this->x != T{0} || this->y != T{0};
        }
        template <isNumWrapper U>
        friend std::ostream &operator<<(std::ostream&, const vector2D<U>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator==(const vector2D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator!=(const vector2D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<(const vector2D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>(const vector2D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<=(const vector2D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>=(const vector2D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator==(const vector2D<U>&, const vector3D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator!=(const vector2D<U>&, const vector3D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<(const vector2D<U>&, const vector3D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>(const vector2D<U>&, const vector3D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<=(const vector2D<U>&, const vector3D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>=(const vector2D<U>&, const vector3D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator==(const vector3D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator!=(const vector3D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<(const vector3D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>(const vector3D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<=(const vector3D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>=(const vector3D<U>&, const vector2D<V>&);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x & value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x ^ value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value % vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m, const vector2D<V> &v)
        -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::cross(const vector2D<U> &vec1, const vector2D<V> &vec2) ->
        vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::cross(const vector3D<U> &vec1, const vector2D<V> &vec2) ->
        vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::cross(const vector2D<U> &vec1, const vector3D<V> &vec2) ->
        vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::cross(const vector3D<U> &vec1, const vector3D<V> &vec2) ->
        vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::distance(const vector2D<U>& vec1, const vector2D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::distance(const vector2D<U>& vec1, const vector3D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::distance(const vector3D<U>& vec1, const vector2D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::distance(const vector3D<U>& vec1, const vector3D<V>& vec2);
        template <isNumWrapper>
        friend class vector2D;
        template <isNumWrapper>
        friend class vector3D;
        template <isNumWrapper, isNumWrapper, isNumWrapper, uint64_t>
        friend class body;
        template <isNumWrapper, isNumWrapper, isNumWrapper>
        friend class ray;
        template <isNumWrapper, isNumWrapper, isNumWrapper>
        friend class camera;
    };
    template <isNumWrapper U>
    std::ostream &operator<<(std::ostream &out, const vector2D<U> &vec) {
        if (vec.y >= 0) {
            return out << +vec.x << "i + " << +vec.y << "j"; // '+' for always printing out numerical value
        }
        return out << +vec.x << "i - " << +(-vec.y) << "j";
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator==(const vector2D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.x == vec2.x && vec1.y == vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator!=(const vector2D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.x != vec2.x || vec1.y != vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator<(const vector2D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.magnitude() < vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator>(const vector2D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.magnitude() > vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator<=(const vector2D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.magnitude() <= vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator>=(const vector2D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.magnitude() >= vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x + vec2.x)> {
        return vector2D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x - vec2.x)> {
        return vector2D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x << vec2.x)> {
        return vector2D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x >> vec2.x)> {
        return vector2D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x | vec2.x)> {
        return vector2D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x & vec2.x)> {
        return vector2D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x ^ vec2.x)> {
        return vector2D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x + value)> {
        return vector2D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x - value)> {
        return vector2D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x * value)> {
        return vector2D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x << value)> {
        return vector2D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x >> value)> {
        return vector2D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x | value)> {
        return vector2D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x & value)> {
        return vector2D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x ^ value)> {
        return vector2D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value + vec1.x)> { // could just return vec op value,
        return vector2D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y); // but would be an extra func. call
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value - vec1.x)> {
        return vector2D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value * vec1.x)> {
        return vector2D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value ^ vec1.x)> {
        return vector2D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &m, const vector2D<V> &v)
    -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
            requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)> {
        return vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>(v).apply(m);
    }
    template <isNumWrapper T>
    class vector3D : public vector2D<T> {
    public:
        static inline const vector3D<T> zero{T{}, T{}, T{}};
        static inline const vector3D<T> left{T{-1}, T{}};
        static inline const vector3D<T> right{T{1}, T{}};
        static inline const vector3D<T> forwards{T{}, T{1}};
        static inline const vector3D<T> backwards{T{}, T{-1}};
        static inline const vector3D<T> up{T{}, T{}, T{1}};
        static inline const vector3D<T> down{T{}, T{}, T{-1}};
    protected:
        T z{};
    public:
        constexpr vector3D() noexcept = default;
        constexpr vector3D(const vector2D<T> &other) noexcept : vector2D<T>(other) {}
        constexpr vector3D(vector2D<T> &&other) noexcept : vector2D<T>(std::move(other)) {}
        constexpr vector3D(const vector3D<T> &other) noexcept : vector2D<T>(other), z{other.z} {}
        constexpr vector3D(vector3D<T> &&other) noexcept : vector2D<T>(std::move(other)), z{std::move(other.z)} {}
        template <isConvertible<T> U>
        constexpr vector3D(const vector2D<U> &other) noexcept : vector2D<T>(other) {}
        template <isConvertible<T> U>
        constexpr vector3D(const vector2D<U> &&other) noexcept : vector2D<T>(other) {}
        template <isConvertible<T> U>
        constexpr vector3D(const vector3D<U> &other) noexcept : vector2D<T>(other), z{(T) other.z} {}
        template <isConvertible<T> U>
        constexpr vector3D(const vector3D<U> &&other) noexcept : vector2D<T>(other), z{(T) other.z} {}
        constexpr vector3D(const T &x_component, const T &y_component, const T &z_component) noexcept :
        vector2D<T>{x_component, y_component}, z{z_component} {}
        constexpr vector3D(T &&x_component, T &&y_component, T &&z_component) noexcept :
        vector2D<T>{std::move(x_component), std::move(y_component)}, z{std::move(z_component)} {}
        constexpr vector3D(const T &x_component, const T &y_component) noexcept : vector2D<T>{x_component,y_component}{}
        constexpr vector3D(T &&x_component, T &&y_component) noexcept :
        vector2D<T>{std::move(x_component), std::move(y_component)} {}
        String str(unsigned char f_p_dec_places = 5) const override {
            String s;
            s.append_back(this->x, f_p_dec_places);
            if (this->y < 0) {
                s.append_back("i - ").append_back(-this->y, f_p_dec_places).append_back("j");
            }
            else {
                s.append_back("i + ").append_back(this->y, f_p_dec_places).append_back("j");
            }
            if (this->z < 0) {
                s.append_back(" - ").append_back(-this->z, f_p_dec_places).append_back("k");
            }
            else {
                s.append_back(" + ").append_back(this->z, f_p_dec_places).append_back("k");
            }
            return s;
        }
        vector3D<T> &set_x(T value) noexcept override {
            this->x = value;
            return *this;
        }
        vector3D<T> &set_y(T value) noexcept override {
            this->y = value;
            return *this;
        }
        virtual vector3D<T> &set_z(T value) noexcept {
            this->z = value;
            return *this;
        }
        vector3D<T> &assign(const T& xcomp, const T& ycomp) override {
            this->x = xcomp;
            this->y = ycomp;
            this->z = T{};
            return *this;
        }
        vector3D<T> &assign(T&& xcomp, T&& ycomp) noexcept override {
            this->x = std::move(xcomp);
            this->y = std::move(ycomp);
            this->z = T{};
            return *this;
        }
        vector3D<T> &assign(const T& xcomp, const T& ycomp, const T& zcomp) {
            this->x = xcomp;
            this->y = ycomp;
            this->z = zcomp;
            return *this;
        }
        vector3D<T> &assign(T&& xcomp, T&& ycomp, T&& zcomp) noexcept {
            this->x = std::move(xcomp);
            this->y = std::move(ycomp);
            this->z = std::move(zcomp);
            return *this;
        }
        const T &get_z() const {
            return this->z;
        }
        vector3D<T> &make_zero() override {
            this->x = T{0};
            this->y = T{0};
            this->z = T{0};
            return *this;
        }
        // template <isNumWrapper U>
        vector3D<T> &set_length(const T &new_length) noexcept override {
            long double frac = new_length/this->magnitude();
            this->x *= frac;
            this->y *= frac;
            this->z *= frac;
            return *this;
        }
        // template <isNumWrapper U>
        vector3D<T> &set_length(const T &&new_length) noexcept override {
            long double frac = new_length/this->magnitude();
            this->x *= frac;
            this->y *= frac;
            this->z *= frac;
            return *this;
        }
        inline std::tuple<T, T, T> to_tuple() const {
            return {this->x, this->y, this->z};
        }
        long double magnitude() const noexcept override {
            return sqrtl(static_cast<long double>(this->x)*static_cast<long double>(this->x) +
                         static_cast<long double>(this->y)*static_cast<long double>(this->y) +
                         static_cast<long double>(this->z)*static_cast<long double>(this->z));
        }
        vector3D<long double> unit_vector() const noexcept { // no way around the method hiding apart from changing the
            if (this->is_zero())                             // function's name
                return {};
            return *this / this->magnitude();
        }
        vector3D<T> &normalise() noexcept override { // makes the vector a unit vector
            if (this->is_zero())
                return *this;
            long double mag = this->magnitude();
            this->x /= mag;
            this->y /= mag;
            this->z /= mag;
            return *this;
        }
        vector3D<T> z_projection() const {
            return {T{0}, T{0}, this->z};
        }
        vector2D<T> xy_projection() const {
            return {this->x, this->y};
        }
        vector3D<T> xz_projection() const {
            return {this->x, T{0}, this->z};
        }
        vector3D<T> yz_projection() const {
            return {T{0}, this->y, this->z};
        }
        template <isNumWrapper U, isNumWrapper V>
        bool between(const vector2D<U> &v1, const vector3D<V> &v2) const noexcept {
            /* Tests if the vector lies within the 3D cube spanned by v1 & v2 (v1 == bottom corner, v2 == top corner) */
            return v1.x <= this->x && this->x < v2.x &&
                   v1.y <= this->y && this->y < v2.y &&
                   this->z >= T{} && this->z < v2.z;
        }
        template <isNumWrapper U, isNumWrapper V>
        bool between(const vector3D<U> &v1, const vector2D<V> &v2) const noexcept {
            /* Tests if the vector lies within the 3D cube spanned by v1 & v2 (v1 == bottom corner, v2 == top corner) */
            return v1.x <= this->x && this->x < v2.x &&
                   v1.y <= this->y && this->y < v2.y &&
                   this->z >= v1.z && this->z < T{};
        }
        template <isNumWrapper U, isNumWrapper V>
        bool between(const vector3D<U> &v1, const vector3D<V> &v2) const noexcept {
            /* Tests if the vector lies within the 3D cube spanned by v1 & v2 (v1 == bottom corner, v2 == top corner) */
            return v1.x <= this->x && this->x < v2.x &&
                   v1.y <= this->y && this->y < v2.y &&
                   v1.z <= this->z && this->z < v2.z;
        }
        virtual vector3D<T> &rotate(const long double &angle_in_rad = PI, char about = 'z') {
            return this->apply(matrix<long double>::get_3D_rotation_matrix(angle_in_rad, about));
        }
        virtual vector3D<T> &rodrigues_rotate(const vector3D<T> &about, long double angle = PI) noexcept {
            if (about.is_zero())
                return *this;
            vector3D<long double> about_unit = about.unit_vector();
            vector3D<long double> cpy(*this); // this assignment is done in case T is not a long double
            return *this = cosl(angle)*cpy + sinl(angle)*vec_ops::cross(about_unit, cpy) +
                    (about_unit*cpy)*(1 - cosl(angle))*about_unit;
        }
        virtual vector3D<T> &rotate_to(const vector3D<T> &new_direction) noexcept {
            if (this->is_zero())
                return *this;
            return this->rodrigues_rotate(vec_ops::cross(*this, new_direction),
                                          vec_ops::angle_between(*this, new_direction));
        }
        virtual vector3D<T> &rodrigues_rotate(const vector2D<T> &about, long double angle = PI) noexcept {
            if (about.is_zero())
                return *this;
            vector2D<long double> about_unit = about.unit_vector();
            vector3D<long double> cpy(*this);
            return *this = cosl(angle)*cpy + sinl(angle)*vec_ops::cross(about_unit, cpy) +
                    (about_unit*cpy)*((1 - cosl(angle))*about_unit);
        }
        vector3D<T> &rotate_to(const vector2D<T> &new_direction) noexcept override {
            if (this->is_zero())
                return *this;
            return this->rodrigues_rotate(vec_ops::cross(*this, new_direction),
                                          vec_ops::angle_between(*this, new_direction));
        }
        bool is_zero() const noexcept override {
            return this->x == T{0} && this->y == T{0} && this->z == T{0};
        }
        vector3D<T> &apply(const matrix<T> &transform) override {
            if (transform.mat.size() != 3 || transform.mat[0].size() != 3)
                throw invalid_matrix_format("Only 3x3 matrices can be applied to a vector3D object.");
            T org_x = this->x;
            T org_y = this->y;
            T org_z = this->z;
            this->x = transform.mat[0][0]*org_x + transform.mat[0][1]*org_y + transform.mat[0][2]*org_z;
            this->y = transform.mat[1][0]*org_x + transform.mat[1][1]*org_y + transform.mat[1][2]*org_z;
            this->z = transform.mat[2][0]*org_x + transform.mat[2][1]*org_y + transform.mat[2][2]*org_z;
            return *this;
        }
        vector3D<T> copy() const {
            return {this->x, this->y, this->z};
        }
        T &operator[](unsigned char index) override {
            switch (index) {
                case 0:
                    return this->x;
                case 1:
                    return this->y;
                case 2:
                    return this->z;
                default:
                    throw std::invalid_argument("Only the indices '0', '1' and '2' are possible.\n");
            }
        }
        const T &operator[](unsigned char index) const override {
            switch (index) { // must repeat above code, as a non-const method cannot be called from a const method
                case 0:
                    return this->x;
                case 1:
                    return this->y;
                case 2:
                    return this->z;
                default:
                    throw std::invalid_argument("Only the indices '0', '1' and '2' are possible.\n");
            }
        }
        vector3D<T> operator-() const {
            return {-this->x, -this->y, -this->z};
        }
        vector3D<T> &operator+=(const vector2D<T> &other) noexcept override {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        vector3D<T> &operator-=(const vector2D<T> &other) noexcept override {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        vector3D<T> &operator*=(const vector2D<T> &other) noexcept override {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0}; // given this is the scalar product, mul. of a 3D vec. by a 2D vec. yields a 2D vec (or a 3D
            return *this; // vec. with a zero z-component)
        }
        vector3D<T> &operator/=(const vector2D<T> &other) override {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x; // division between vectors makes no physical sense anyway, so there is no reason to
            this->y /= other.y; // 'divide' by the zero z-component of the 2D vector
            return *this; // these 'non-physical' methods are only here for programming convenience... it should be
        }                 // remembered that they do not represent valid mathematical operations
        vector3D<T> &operator%=(const vector2D<T> &other) requires isIntegralNumWrapper<T> { // hides the parent method,
            if (this->x == T{0} || this->y == T{0}) { // but nothing can be done, since virtual functions cannot have
                throw division_by_zero();             // requires clauses
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        vector3D<T> &operator+=(const vector2D<T> &&other) noexcept override {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        vector3D<T> &operator-=(const vector2D<T> &&other) noexcept override {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        vector3D<T> &operator*=(const vector2D<T> &&other) noexcept override {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0};
            return *this;
        }
        vector3D<T> &operator/=(const vector2D<T> &&other) override {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        vector3D<T> &operator%=(const vector2D<T> &&other) requires isIntegralNumWrapper<T> {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator+=(const vector2D<U> &other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator-=(const vector2D<U> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator*=(const vector2D<U> &other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0};
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator/=(const vector2D<U> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector3D<T> &operator%=(const vector2D<U> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator<<=(const vector2D<U> &other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y; // z << 0 would just equal z, so not included
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator>>=(const vector2D<U> &other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y; // z >> 0 would just equal z, so not included
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator|=(const vector2D<U> &other) noexcept {
            this->x |= other.x;
            this->y |= other.y; // z | 0 would equal z, so not included
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator&=(const vector2D<U> &other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            this->z = T{0}; // bitwise AND with zero will always equal zero
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator^=(const vector2D<U> &other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y; // z ^ 0 would equal z, so left unchanged
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator+=(const vector2D<U> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator-=(const vector2D<U> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator*=(const vector2D<U> &&other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0};
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator/=(const vector2D<U> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector3D<T> &operator%=(const vector2D<U> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator<<=(const vector2D<U> &&other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator>>=(const vector2D<U> &&other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator|=(const vector2D<U> &&other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator&=(const vector2D<U> &&other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            this->z = T{0};
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator^=(const vector2D<U> &&other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        vector3D<T> &operator+=(const vector3D<T> &other) noexcept {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }
        vector3D<T> &operator-=(const vector3D<T> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }
        vector3D<T> &operator*=(const vector3D<T> &other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }
        vector3D<T> &operator/=(const vector3D<T> &other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            this->z /= other.z;
            return *this;
        }
        vector3D<T> &operator%=(const vector3D<T> &other) requires isIntegralNumWrapper<T> {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            this->z %= other.z;
            return *this;
        }
        vector3D<T> &operator+=(const vector3D<T> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }
        vector3D<T> &operator-=(const vector3D<T> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }
        vector3D<T> &operator*=(const vector3D<T> &&other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }
        vector3D<T> &operator/=(const vector3D<T> &&other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            this->z /= other.z;
            return *this;
        }
        vector3D<T> &operator%=(const vector3D<T> &&other) requires isIntegralNumWrapper<T> {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            this->z %= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator+=(const vector3D<U> &other) noexcept {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator-=(const vector3D<U> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator*=(const vector3D<U> &other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator/=(const vector3D<U> &other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            this->z /= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector3D<T> &operator%=(const vector3D<U> &other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            this->z %= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator<<=(const vector3D<U> &other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            this->z <<= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator>>=(const vector3D<U> &other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            this->z >>= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator|=(const vector3D<U> &other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            this->z |= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator&=(const vector3D<U> &other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            this->z &= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator^=(const vector3D<U> &other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            this->z ^= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator+=(const vector3D<U> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator-=(const vector3D<U> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator*=(const vector3D<U> &&other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator/=(const vector3D<U> &&other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            this->z /= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector3D<T> &operator%=(const vector3D<U> &&other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            this->z %= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator<<=(const vector3D<U> &&other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            this->z <<= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator>>=(const vector3D<U> &&other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            this->z >>= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator|=(const vector3D<U> &&other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            this->z |= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator&=(const vector3D<U> &&other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            this->z &= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator^=(const vector3D<U> &&other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            this->z ^= other.z;
            return *this;
        }
        template <isConvertible<T> U>
        vector3D<T> &operator*=(const U &scalar) {
            this->x *= scalar;
            this->y *= scalar;
            this->z *= scalar;
            return *this;
        }
        template <isConvertible<T> U>
        vector3D<T> &operator/=(const U &scalar) {
            this->x /= scalar;
            this->y /= scalar;
            this->z /= scalar;
            return *this;
        }
        vector3D<T> operator++(int) noexcept {
            vector3D<T> retvec{this->x, this->y, this->z};
            ++this->x;
            ++this->y;
            ++this->z;
            return retvec;
        }
        vector3D<T> &operator++() noexcept override {
            ++this->x;
            ++this->y;
            ++this->z;
            return *this;
        }
        vector3D<T> operator--(int) noexcept {
            vector3D<T> retvec{this->x, this->y, this->z};
            --this->x;
            --this->y;
            --this->z;
            return retvec;
        }
        vector3D<T> &operator--() noexcept override {
            --this->x;
            --this->y;
            --this->z;
            return *this;
        }
        vector3D<T> &operator=(const vector<T> &other) noexcept override {
            if (&other == this)
                return *this;
            try {
                const vector3D<T> &oth = dynamic_cast<const vector3D<T>&>(other);
                this->x = oth.x;
                this->y = oth.y;
                this->z = oth.z;
            } catch (const std::bad_cast &bc) { // if cast to vector3D fails, try to vector2D
                try {
                    const vector2D<T> &oth_v2 = dynamic_cast<const vector2D<T>&>(other);
                    this->x = oth_v2.x;
                    this->y = oth_v2.y;
                } catch (std::bad_cast &bc2) {} // this would only occur if a user has created a subclass of vector<T>
            }
            return *this;
        }
        virtual vector3D<T> &operator=(const vector3D<T> &other) noexcept {
            if (&other == this)
                return *this;
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            return *this;
        }
        template <isConvertible<T> U>
        vector3D<T> &operator=(const vector3D<U> &other) noexcept {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            return *this;
        }
        template <isConvertible<T> U>
        vector3D<T> &operator=(const vector3D<U> &&other) noexcept {
            return *this = other;
        }
        vector3D<T> &operator=(vector<T> &&other) noexcept override {
            if (&other == this)
                return *this;
            try {
                const vector3D<T> &oth = dynamic_cast<const vector3D<T>&>(other);
                this->x = std::move(oth.x);
                this->y = std::move(oth.y);
                this->z = std::move(oth.z);
            } catch (const std::bad_cast &bc) {
                try {
                    const vector2D<T> &oth_v2 = dynamic_cast<const vector2D<T>&>(other);
                    this->x = std::move(oth_v2.x);
                    this->y = std::move(oth_v2.y);
                } catch (std::bad_cast &bc2) {}
            }
            return *this;
        }
        virtual vector3D<T> &operator=(vector3D<T> &&other) noexcept {
            this->x = std::move(other.x);
            this->y = std::move(other.y);
            this->z = std::move(other.z);
            return *this;
        }
        vector3D<T> &operator=(const T &value) override {
            this->x = value;
            this->y = value;
            this->z = value;
            return *this;
        }
        vector3D<T> &operator=(T &&value) noexcept override {
            this->x = value;
            this->y = value;
            this->z = std::move(value);
            return *this;
        }
        template <isConvertible<T> U>
        vector3D<T> &operator=(const U &value) {
            this->x = value;
            this->y = value;
            this->z = value;
            return *this;
        }
        operator bool() const noexcept override {
            return this->x != T{0} || this->y != T{0} || this->z != T{0};
        }
        template <isNumWrapper U>
        friend std::ostream &operator<<(std::ostream &out, const vector3D<U> &vec);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator==(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator!=(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<=(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>=(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator==(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator!=(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<=(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>=(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator==(const vector3D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator!=(const vector3D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<(const vector3D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>(const vector3D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator<=(const vector3D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator>=(const vector3D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x & value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x ^ value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value % vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m, const vector3D<V> &v)
        -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::cross(const vector2D<U> &vec1, const vector2D<V> &vec2) ->
        vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::cross(const vector3D<U> &vec1, const vector2D<V> &vec2) ->
        vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::cross(const vector2D<U> &vec1, const vector3D<V> &vec2) ->
        vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::cross(const vector3D<U> &vec1, const vector3D<V> &vec2) ->
        vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::distance(const vector2D<U>& vec1, const vector2D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::distance(const vector2D<U>& vec1, const vector3D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::distance(const vector3D<U>& vec1, const vector2D<V>& vec2);
        template <isNumWrapper U, isNumWrapper V>
        friend auto vec_ops::distance(const vector3D<U>& vec1, const vector3D<V>& vec2);
        template <isNumWrapper>
        friend class vector2D;
        template <isNumWrapper>
        friend class vector3D;
        template <isNumWrapper, isNumWrapper, isNumWrapper, uint64_t>
        friend class body;
        template <isNumWrapper, isNumWrapper, isNumWrapper>
        friend class ray;
        template <isNumWrapper, isNumWrapper, isNumWrapper>
        friend class camera;
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, uint64_t, uint64_t, bool>
        friend class system;
        template <isNumWrapper, isNumWrapper, isNumWrapper, uint64_t>
        friend class bh_cube;
        template <isNumWrapper, isNumWrapper, isNumWrapper, uint64_t>
        friend class bh_tree;
    };
    template <isNumWrapper U>
    std::ostream &operator<<(std::ostream &out, const vector3D<U> &vec) {
        out << vec.x << "i ";
        if (vec.y < 0) {
            out << "- " << +(-vec.y) << "j ";
        }
        else {
            out << "+ " << +vec.y << "j ";
        }
        if (vec.z < 0) {
            out << "- " << +(-vec.z) << "k";
        }
        else {
            out << "+ " << +vec.z << "k";
        }
        return out;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator==(const vector2D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.x == vec2.x && vec1.y == vec2.y && vec2.z == 0;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator!=(const vector2D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.x != vec2.x || vec1.y != vec2.y || vec2.z != 0;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator<(const vector2D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.magnitude() < vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator>(const vector2D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.magnitude() > vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator<=(const vector2D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.magnitude() <= vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator>=(const vector2D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.magnitude() >= vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator==(const vector3D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.x == vec2.x && vec1.y == vec2.y && vec1.z == 0;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator!=(const vector3D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.x != vec2.x || vec1.y != vec2.y || vec1.z != 0;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator<(const vector3D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.magnitude() < vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator>(const vector3D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.magnitude() > vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator<=(const vector3D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.magnitude() <= vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator>=(const vector3D<U> &vec1, const vector2D<V> &vec2) {
        return vec1.magnitude() >= vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator==(const vector3D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.x == vec2.x && vec1.y == vec2.y && vec1.z == vec2.z;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator!=(const vector3D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.x != vec2.x || vec1.y != vec2.y || vec1.z != vec2.z;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator<(const vector3D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.magnitude() < vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator>(const vector3D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.magnitude() > vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator<=(const vector3D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.magnitude() <= vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator>=(const vector3D<U> &vec1, const vector3D<V> &vec2) {
        return vec1.magnitude() >= vec2.magnitude();
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z / vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z % vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z << vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z >> vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z | vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, vec1.z & vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z ^ vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x + value)> {
        return vector3D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value, vec1.z + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x - value)> {
        return vector3D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value, vec1.z - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x * value)> {
        return vector3D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value, vec1.z * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value, vec1.z / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value, vec1.z % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x << value)> {
        return vector3D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value, vec1.z << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x >> value)> {
        return vector3D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value, vec1.z >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x | value)> {
        return vector3D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value, vec1.z | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x & value)> {
        return vector3D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value, vec1.z & value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x ^ value)> {
        return vector3D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value, vec1.z ^ value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value + vec1.x)> { // could just return vec op value,
        return vector3D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y, value + vec1.z); // but would be an extra func. call
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value - vec1.x)> {
        return vector3D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y, value - vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value * vec1.x)> {
        return vector3D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y, value * vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y, value / vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y, value % vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value << vec1.x)> {
        return vector3D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y, value << vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value >> vec1.x)> {
        return vector3D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y, value >> vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value | vec1.x)> {
        return vector3D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y, value | vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value & vec1.x)> {
        return vector3D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y, value & vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value ^ vec1.x)> {
        return vector3D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y, value ^ vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &m, const vector3D<V> &v)
    -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
    requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)> {
        return vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>(v).apply(m);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec2.z);
    }
    namespace vec_ops {
        template <isNumWrapper U, isNumWrapper V>
        auto distance(const vector2D<U>& vec1, const vector2D<V>& vec2) {
            auto dist_x = vec2.x - vec1.x;
            auto dist_y = vec2.y - vec1.y;
            return sqrtl(dist_x*dist_x + dist_y*dist_y);
        }
        template <isNumWrapper U, isNumWrapper V>
        auto distance(const vector2D<U>& vec1, const vector3D<V>& vec2) {
            auto dist_x = vec2.x - vec1.x;
            auto dist_y = vec2.y - vec1.y;
            return sqrtl(dist_x*dist_x + dist_y*dist_y + vec2.z*vec2.z);
        }
        template <isNumWrapper U, isNumWrapper V>
        auto distance(const vector3D<U>& vec1, const vector2D<V>& vec2) {
            auto dist_x = vec2.x - vec1.x;
            auto dist_y = vec2.y - vec1.y;
            return sqrtl(dist_x*dist_x + dist_y*dist_y + vec1.z*vec1.z);
        }
        template <isNumWrapper U, isNumWrapper V>
        auto distance(const vector3D<U>& vec1, const vector3D<V>& vec2) {
            auto dist_x = vec2.x - vec1.x;
            auto dist_y = vec2.y - vec1.y;
            auto dist_z = vec2.z - vec1.z;
            return sqrtl(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
        }
        template <isNumWrapper U, isNumWrapper V>
        auto cross(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
            return vector3D<decltype(vec1.x * vec2.y)>(0, 0, vec1.x*vec2.y - vec1.y*vec2.x);
        }
        template <isNumWrapper U, isNumWrapper V>
        auto cross(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
            return vector3D<decltype(vec1.x * vec2.y)>(-vec1.z*vec2.y, // vec2.z is zero
                                                       vec1.z*vec2.x,
                                                       vec1.x*vec2.y - vec1.y*vec2.x);
        }
        template <isNumWrapper U, isNumWrapper V>
        auto cross(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
            return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z, // vec1.z is zero
                                                       -vec1.x*vec2.z,
                                                       vec1.x*vec2.y - vec1.y*vec2.x);
        }
        template <isNumWrapper U, isNumWrapper V>
        auto cross(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
            return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z - vec1.z*vec2.y,
                                                       vec1.z*vec2.x - vec1.x*vec2.z,
                                                       vec1.x*vec2.y - vec1.y*vec2.x);
        }
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V,long double>
        long double angle_between(const vector2D<U> &vec1, const vector2D<V> &vec2) {
            long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
            return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : acosl(cosine_of_angle));
        } // have to check for less than -1 or greater than 1 because of f. p. rounding errors
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V,long double>
        long double angle_between(const vector3D<U> &vec1, const vector2D<V> &vec2) {
            long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
            return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : acosl(cosine_of_angle));
        }
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V,long double>
        long double angle_between(const vector2D<U> &vec1, const vector3D<V> &vec2) {
            long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
            return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : acosl(cosine_of_angle));
        }
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V,long double>
        long double angle_between(const vector3D<U> &vec1, const vector3D<V> &vec2) {
            long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
            return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : acosl(cosine_of_angle));
        }
    }
    typedef vector2D<long double> vec2;
    typedef vector3D<long double> vec3;
}
#endif
