//
// Created by Gregor Hartl Watters on 04/11/2022.
//

#ifndef GREGASTRO_H
#define GREGASTRO_H

#ifndef __cplusplus
#error "The gregastro.hpp header file is a C++ header file only."
#endif

#include <functional>
#include <random>
#include <thread>
#include <latch>
#include <condition_variable>

#include "gregsys.hpp"
#include "gregbmp.hpp"

namespace gtd {
    template <typename T>
    concept container = requires (T val, T val2, const T val_c) {
        typename T::value_type;
        typename T::reference;
        typename T::const_reference;
        typename T::iterator;
        typename T::const_iterator;
        typename T::difference_type;
        typename T::size_type;
        requires std::copy_constructible<typename T::value_type>;
        requires std::default_initializable<typename T::value_type>;
        requires std::signed_integral<typename T::difference_type>;
        requires std::unsigned_integral<typename T::size_type>;
        requires std::destructible<typename T::value_type>;
        requires std::convertible_to<typename T::iterator, typename T::const_iterator>;
        requires std::same_as<typename T::difference_type,
                              typename std::iterator_traits<typename T::iterator>::difference_type>;
        requires std::same_as<typename T::difference_type,
                              typename std::iterator_traits<typename T::const_iterator>::difference_type>;
        {val.begin()} -> std::same_as<typename T::iterator>;
        {val.end()} -> std::same_as<typename T::iterator>;
        {val.cbegin()} -> std::same_as<typename T::const_iterator>;
        {val.cend()} -> std::same_as<typename T::const_iterator>;
        {val_c.begin()} -> std::same_as<typename T::const_iterator>;
        {val_c.end()} -> std::same_as<typename T::const_iterator>;
        {val_c.cbegin()} -> std::same_as<typename T::const_iterator>;
        {val_c.cend()} -> std::same_as<typename T::const_iterator>;
        {val == val_c} -> std::convertible_to<bool>;
        {val != val_c} -> std::convertible_to<bool>;
        val.swap(val2);
        std::swap(val, val2);
        {val.size()} -> std::same_as<typename T::size_type>;
        {val.max_size()} -> std::same_as<typename T::size_type>;
        {val.empty()} -> std::convertible_to<bool>;
    };
    template <typename T>
    concept colorGenerator = requires (T func) {
        {func()} -> std::same_as<color>;
    };
    class astro_error : public std::logic_error {
    public:
        astro_error() : std::logic_error{"An error occurred due to invalid evaluation of conditions.\n"} {}
        explicit astro_error(const char *message) : std::logic_error{message} {}
    };
    class camera_error : public astro_error {
    public:
        camera_error() : astro_error{"Error occured with camera.\n"} {}
        explicit camera_error(const char *message) : astro_error{message} {}
    };
    class no_direction : public astro_error {
    public:
        no_direction() : astro_error{"Zero vector error.\n"} {}
        explicit no_direction(const char *message) : astro_error{message} {}
    };
    class zero_image_distance : public camera_error {
    public:
        zero_image_distance() : camera_error{"The image distance of a camera cannot be zero.\n"} {}
        explicit zero_image_distance(const char *message) : camera_error{message} {}
    };
    class no_intersection : public astro_error {
    public:
        no_intersection() : astro_error{"The ray has not intersected with any body.\n"} {}
        explicit no_intersection(const char *message) : astro_error{message} {}
    };
    class camera_dimensions_error : public camera_error {
    public:
        camera_dimensions_error() : camera_error{"The dimensions of the camera are invalid.\n"} {}
        explicit camera_dimensions_error(const char *message) : camera_error{message} {}
    };
    struct image_dimensions {
        unsigned int x;
        unsigned int y;
    };
    template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper PosT, isNumWrapper, isNumWrapper, isNumWrapper,
              uint64_t> requires isConvertible<long double, PosT>
    class star;
    template <isNumWrapper PosT = long double, isNumWrapper DirT = long double, isNumWrapper LenT = long double>
    class ray : public vector3D<DirT> {
        /* A ray object represents a ray in simulated 3D space - whether this be an actual ray of light (from a light
         * source), or a "line-of-sight" ray, such as those shot out from a camera object (see below). A ray has an
         * origin (its position) a direction, and a parameter named 'l', which represents the length of the ray. A ray
         * with length zero is taken as one which continues on to infinity. Since the ray class is a child class of
         * vector3D<DirT>, the internal parent vector3D<DirT> object is made to be the direction, whilst a ray's
         * position is represented by a separate vector3D<PosT> object. */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        using ray_t = ray<PosT, DirT, LenT>;
        static inline const vector3D<DirT> zero{};
        vector3D<PosT> pos; // origin of the ray
        LenT l{0}; // length of the ray
        uint64_t ibod_id = -1; // intersected body id - a ray stores the id of the body it intersects with
        bool intersected = false;
        mutable vector3D<DirT> normal; // set during the reaches() method - used internally by astro_scene
        ray<PosT, DirT, LenT> &calc() {
            if (*this == zero)
                throw no_direction("A ray's direction cannot be a zero vector, "
                                   "or else the ray would have nowhere to point to.\n");
            this->normalise();
            return *this;
        }
        template <isNumWrapper T>
        static inline const T &&minimum(const T &&first, const T &&second) {
            return first <= second ? std::move(first) : std::move(second);
        }
        template <typename ...Args>
        bool contains_id(const uint64_t &id, const uint64_t &first, Args... args){
            if constexpr (sizeof...(args) == 0)
                return id == first;
            return id == first || contains_id(id, args...);
        };
        /* the (seemingly redundant) function below is used in a fold expression further down (the syntax requires it)*/
        static inline auto equal =
                [](const uint64_t &arg1, const uint64_t &arg2){return arg1 == arg2;};
    public:
        ray() : vector3D<DirT>{0, 0, -1} {calc();}
        ray(const vector2D<DirT> &direction) noexcept : vector3D<DirT>(direction) {calc();}
        ray(vector2D<DirT> &&direction) noexcept : vector3D<DirT>(std::move(direction)) {calc();}
        ray(const vector3D<DirT> &direction) noexcept : vector3D<DirT>(direction) {calc();}
        ray(vector3D<DirT> &&direction) noexcept : vector3D<DirT>(std::move(direction)) {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector2D<U> &direction) noexcept : vector3D<DirT>(direction) {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector2D<U> &&direction) noexcept : vector3D<DirT>(direction) {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector3D<U> &direction) noexcept : vector3D<DirT>(direction) {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector3D<U> &&direction) noexcept : vector3D<DirT>(direction) {calc();}
        ray(const DirT &x_component, const DirT &y_component, const DirT &z_component) noexcept :
                vector3D<DirT>{x_component, y_component, z_component} {calc();}
        ray(DirT &&x_component, DirT &&y_component, DirT &&z_component) noexcept :
        vector3D<DirT>{std::move(x_component), std::move(y_component), std::move(z_component)} {calc();}
        ray(const DirT &x_component, const DirT &y_component) noexcept :
        vector3D<DirT>{x_component, y_component} {calc();}
        ray(DirT &&x_component, DirT &&y_component) noexcept :
        vector3D<DirT>{std::move(x_component), std::move(y_component)} {calc();}
        ray(const vector2D<PosT> &origin, const vector2D<DirT> &direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        ray(vector2D<PosT> &&origin, vector2D<DirT> &&direction) noexcept :
        vector3D<DirT>(std::move(direction)), pos{std::move(origin)} {calc();}
        ray(const vector2D<PosT> &origin, vector2D<DirT> &&direction) noexcept :
                vector3D<DirT>(std::move(direction)), pos{origin} {calc();}
        ray(vector2D<PosT> &&origin, const vector2D<DirT> &direction) noexcept :
                vector3D<DirT>(direction), pos{std::move(origin)} {calc();}
        ray(const vector3D<PosT> &origin, const vector3D<DirT> &direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        ray(vector3D<PosT> &&origin, vector3D<DirT> &&direction) noexcept :
        vector3D<DirT>(std::move(direction)), pos{std::move(origin)} {calc();}
        ray(const vector3D<PosT> &origin, vector3D<DirT> &&direction) noexcept :
                vector3D<DirT>(std::move(direction)), pos{origin} {calc();}
        ray(vector3D<PosT> &&origin, const vector3D<DirT> &direction) noexcept :
                vector3D<DirT>(direction), pos{std::move(origin)} {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector2D<U> &origin, const vector2D<U> &direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector2D<U> &&origin, const vector2D<U> &&direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector3D<U> &origin, const vector3D<U> &direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector3D<U> &&origin, const vector3D<U> &&direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        ray(const DirT &x_pos_comp, const DirT &y_pos_comp, const DirT &z_pos_comp,
            const DirT &x_dir_comp, const DirT &y_dir_comp, const DirT &z_dir_comp) noexcept :
                vector3D<DirT>{x_dir_comp, y_dir_comp, z_dir_comp}, pos{x_pos_comp, y_pos_comp, z_pos_comp} {calc();}
        ray(DirT &&x_pos_comp, DirT &&y_pos_comp, DirT &&z_pos_comp,
            DirT &&x_dir_comp, DirT &&y_dir_comp, DirT &&z_dir_comp) noexcept :
                vector3D<DirT>{std::move(x_dir_comp), std::move(y_dir_comp), std::move(z_dir_comp)},
                pos{std::move(x_pos_comp), std::move(y_pos_comp), std::move(z_pos_comp)} {calc();}
        ray(const DirT &x_pos_comp, const DirT &y_pos_comp, const DirT &x_dir_comp, const DirT &y_dir_comp) noexcept :
                vector3D<DirT>{x_dir_comp, y_dir_comp}, pos{x_pos_comp, y_pos_comp} {calc();}
        ray(DirT &&x_pos_comp, DirT &&y_pos_comp, DirT &&x_dir_comp, DirT &&y_dir_comp) noexcept :
                vector3D<DirT>{std::move(x_dir_comp), std::move(y_dir_comp)},
                pos{std::move(x_pos_comp), std::move(y_pos_comp)} {calc();}
        ray(const vector3D<PosT> &origin, DirT &&x_dir_comp, DirT &&y_dir_comp, DirT &&z_dir_comp) noexcept :
                vector3D<DirT>{std::move(x_dir_comp), std::move(y_dir_comp), std::move(z_dir_comp)},
                pos{origin} {calc();}
        void set_pos(const vector3D<PosT> &new_pos) {
            pos = new_pos;
        }
        void set_pos(vector3D<PosT> &&new_pos) {
            pos = std::move(new_pos);
        }
        void set_dir(const vector3D<DirT> &new_dir) { // set_dir() methods need not be called (reassign ray directly)
            *this = new_dir;
            calc();
        }
        void set_dir(vector3D<DirT> &&new_dir) {
            *this = std::move(new_dir);
            calc();
        }
        vector3D<DirT> &get_dir() {
            return *this;
        }
        const vector3D<PosT> &get_pos() const noexcept { // useless method - no need to use
            return pos;
        }
        const LenT &length() const noexcept {
            return l;
        }
        uint64_t intersected_body_id() const {
            if (!intersected)
                throw no_intersection();
            return ibod_id;
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t recFreq>
        bool intersects(const body<M, R, T, recFreq> &bod) {
            /* This method calculates whether the ray intersects with a given body, and if it does, stores the
             * intersection distance in l. The 'b' and 'c' variables below are quadratic formula coefficients - the 'a'
             * coefficient is always equal to 1 (since a = (*this)*(*this)), so has been inlined. */
            this->calc();
            auto b = 2*(pos - bod.curr_pos)*(*this); // operator* between two vectors is treated as scalar product
            auto c = pos*pos - 2*bod.curr_pos*pos + bod.curr_pos*bod.curr_pos - bod.radius*bod.radius;
            auto discriminant = b*b - 4*c; // recall that a = 1
            if (discriminant < 0) {
                ibod_id = -1;
                return (intersected = false);
            }
            // l = minimum((-b + sqrtl(discriminant))/2, (-b - sqrtl(discriminant))/2);
            l = (-b - sqrtl(discriminant))/2; // this will always be the shorter distance
            if (l < 0) { // if l is negative, that means the ray projected backwards intersected, which is not valid
                ibod_id = -1;
                l = 0;
                return (intersected = false);
            }
            ibod_id = bod.id;
            return (intersected = true);
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t recFreq>
        bool intersects(const body<M, R, T, recFreq> &&bod) {
            return intersects(bod);
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t recFreq, container C, typename ... except>
        requires ((std::same_as<typename C::mapped_type, const body<M, R, T, recFreq>*> ||
                  isConvertible<typename C::mapped_type, const body<M, R, T, recFreq>*>) &&
                  std::same_as<typename C::key_type, const uint64_t>)
        bool intersects(const C &bodies, except ...ids) {
            /* This method is similar to the other intersects() overloads, but it will set l to the distance between the
             * ray's origin and the closest body, and return true (if the ray intersects with at least one body). This
             * is intuitive, since a ray in real space will "terminate" once it hits an object, and will (in general)
             * not continue through the object and then intersect with other objects behind it. */
            intersected = false;
            bool one_intersection = false;
            uint64_t closest_body_id;
            LenT min_l;
            if constexpr (sizeof...(ids) == 0) {
                for (const auto &[id, ptr] : bodies) {
                    if (this->intersects(*ptr)) {
                        if (!one_intersection) {
                            min_l = l;
                            one_intersection = true;
                            closest_body_id = id;// ptr->id;
                            continue;
                        }
                        if (l < min_l) {
                            min_l = l;
                            closest_body_id = id;// ptr->id;
                        }
                    }
                }
            }
            else {
                for (const auto &[id, ptr] : bodies) {
                    // if (contains_id(id, ids...))
                    if (!(equal(id, ids) || ...)) { //required, as the == and || ops cannot both be used in a fold expr.
                        if (this->intersects(*ptr)) {
                            if (!one_intersection) {
                                min_l = l;
                                one_intersection = true;
                                closest_body_id = id;// ptr->id;
                                continue;
                            }
                            if (l < min_l) {
                                min_l = l;
                                closest_body_id = id;// ptr->id;
                            }
                        }
                    }
                }
            }
            if (!one_intersection) {
                l = LenT{0};
                ibod_id = -1;
                return false;
            }
            l = min_l;
            ibod_id = closest_body_id;
            return (intersected = true);
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t rF>
        bool reaches(const body<M, R, T, rF>* b, const vector3D<PosT> &point_on_body) const {
            return ibod_id == b->id && (this->normal = (point_on_body - b->curr_pos))*(*this) <= 0;
        }
        auto end_point() const {
            return pos + (*this)*l;
        }
        ray_t &set_x(DirT value) noexcept override {
            vector3D<DirT>::x = value;
            return *this;
        }
        ray_t &set_y(DirT value) noexcept override {
            vector3D<DirT>::y = value;
            return *this;
        }
        ray_t &set_z(DirT value) noexcept override {
            vector3D<DirT>::z = value;
            return *this;
        }
        ray_t &make_zero() override {
            vector3D<DirT>::x = DirT{0};
            vector3D<DirT>::y = DirT{0};
            vector3D<DirT>::z = DirT{0};
            return *this;
        }
        ray_t &set_length(const DirT &new_length) noexcept override {
            long double frac = new_length/this->magnitude();
            vector3D<DirT>::x *= frac;
            vector3D<DirT>::y *= frac;
            vector3D<DirT>::z *= frac;
            return *this;
        }
        ray_t &set_length(const DirT &&new_length) noexcept override {
            long double frac = new_length/this->magnitude();
            vector3D<DirT>::x *= frac;
            vector3D<DirT>::y *= frac;
            vector3D<DirT>::z *= frac;
            return *this;
        }
        ray_t &normalise() noexcept override { // makes the vector a unit vector
            if (this->is_zero()) {
                return *this;
            }
            long double mag = this->magnitude();
            vector3D<DirT>::x /= mag;
            vector3D<DirT>::y /= mag;
            vector3D<DirT>::z /= mag;
            return *this;
        }
        ray_t &rotate(const long double &angle_in_rad = PI, char about = 'z') override {
            this->apply(matrix<long double>::get_3D_rotation_matrix(angle_in_rad, about));
            return *this;
        }
        // template <isNumWrapper U = T>
        ray_t &rodrigues_rotate(const vector3D<DirT> &about, long double angle = PI) noexcept override {
            vector3D<DirT>::rodrigues_rotate(about, angle);
            return *this;
        }
        // template <isNumWrapper U = T>
        ray_t &rotate_to(const vector3D<DirT> &new_direction) noexcept override {
            if (this->is_zero())
                return *this;
            return this->rodrigues_rotate(vec_ops::cross(*this, new_direction),
                                          vec_ops::angle_between(*this, new_direction));
        }
        // template <isNumWrapper U = T>
        ray_t &rodrigues_rotate(const vector2D<DirT> &about, long double angle = PI) noexcept override {
            vector3D<DirT>::rodrigues_rotate(about, angle);
            return *this;
        }
        // template <isNumWrapper U = T>
        ray_t &rotate_to(const vector2D<DirT> &new_direction) noexcept override {
            if (this->is_zero())
                return *this;
            return this->rodrigues_rotate(vec_ops::cross(*this, new_direction),
                                          vec_ops::angle_between(*this, new_direction));
        }
        ray_t &apply(const matrix<DirT> &transform) override {
            vector3D<DirT>::apply(transform);
            return *this;
        }
        ray_t &operator=(const vector<DirT> &other) noexcept override {
            vector3D<DirT>::operator=(other);
            this->calc();
            return *this;
        }
        ray_t &operator=(const vector3D<DirT> &other) noexcept override {
            vector3D<DirT>::operator=(other);
            this->calc();
            return *this;
        }
        template <isConvertible<DirT> U>
        ray_t &operator=(const vector3D<U> &other) noexcept {
            vector3D<DirT>::operator=(other);
            this->calc();
            return *this;
        } /*
        template <isConvertible<DirT> U>
        ray_t &operator=(const vector3D<U> &&other) noexcept {
            vector3D<DirT>::operator=(other); // no need to use std::move() here
            this->calc();
            return *this;
        } */
        ray_t &operator=(vector<DirT> &&other) noexcept override {
            vector3D<DirT>::operator=(std::move(other));
            this->calc();
            return *this;
        }
        ray_t &operator=(vector3D<DirT> &&other) noexcept override {
            vector3D<DirT>::operator=(std::move(other));
            this->calc();
            return *this;
        }
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU>
        friend std::ostream &operator<<(std::ostream&, const ray<PosU, DirU, LenU>&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU>
        friend std::ostream &operator<<(std::ostream&, const ray<PosU, DirU, LenU>&&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                  isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator==(const ray<PosU, DirU, LenU>&, const ray<PosV, DirV, LenV>&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator==(const ray<PosU, DirU, LenU>&, const ray<PosV, DirV, LenV>&&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator==(const ray<PosU, DirU, LenU>&&, const ray<PosV, DirV, LenV>&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator==(const ray<PosU, DirU, LenU>&&, const ray<PosV, DirV, LenV>&&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator!=(const ray<PosU, DirU, LenU>&, const ray<PosV, DirV, LenV>&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator!=(const ray<PosU, DirU, LenU>&, const ray<PosV, DirV, LenV>&&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator!=(const ray<PosU, DirU, LenU>&&, const ray<PosV, DirV, LenV>&);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator!=(const ray<PosU, DirU, LenU>&&, const ray<PosV, DirV, LenV>&&);
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper>
        friend class light_src;
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                  isNumWrapper, bool, uint64_t>
        friend class astro_scene;
    };
    template <isNumWrapper PosT = long double, isNumWrapper DirT = long double,
              isNumWrapper LenT = long double, isNumWrapper LumT = long double>
    class light_src {
        /* A light_src object represents a point source of light in simulated 3D space, i.e. it is a zero-dimensional
         * source of rays that travel in all directions. */
        vector3D<PosT> pos; // position of the light source
        LumT lum{1}; // luminosity of the light source (W)
    public:
        light_src() = default;
        explicit light_src(const vector3D<PosT> &position) : pos{position} {}
        explicit light_src(vector3D<PosT> &&position) : pos{std::move(position)} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        explicit light_src(const vector3D<U> &position) : pos{position} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        explicit light_src(vector3D<U> &&position) : pos{position} {}
        light_src(const PosT &x, const PosT &y, const PosT &z) : pos{x, y, z} {}
        light_src(PosT &&x, PosT &&y, PosT &&z) : pos{std::move(x), std::move(y), std::move(z)} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(const U &x, const U &y, const U &z) : pos{x, y, z} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(const U &&x, const U &&y, const U &&z) : pos{x, y, z} {}
        light_src(const vector3D<PosT> &position, const LumT &luminosity) : pos{position}, lum{luminosity} {}
        light_src(vector3D<PosT> &&position, const LumT &luminosity) : pos{std::move(position)}, lum{luminosity} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(const vector3D<U> &position, const LumT &luminosity) : pos{position}, lum{luminosity} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(vector3D<U> &&position, const LumT &luminosity) : pos{position}, lum{luminosity} {}
        light_src(const PosT &x, const PosT &y, const PosT &z, const LumT &luminosity) :
        pos{x, y, z}, lum{luminosity} {}
        light_src(PosT &&x, PosT &&y, PosT &&z, const LumT &luminosity) :
        pos{std::move(x), std::move(y), std::move(z)}, lum{luminosity} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(const U &x, const U &y, const U &z, const LumT &luminosity) : pos{x, y, z}, lum{luminosity} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(const U &&x, const U &&y, const U &&z, const LumT &luminosity) : pos{x, y, z}, lum{luminosity} {}
        light_src(const vector3D<PosT> &position, LumT &&luminosity) : pos{position}, lum{std::move(luminosity)} {}
        light_src(vector3D<PosT> &&position, LumT &&luminosity) :
        pos{std::move(position)}, lum{std::move(luminosity)} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(const vector3D<U> &position, LumT &&luminosity) : pos{position}, lum{std::move(luminosity)} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(vector3D<U> &&position, LumT &&luminosity) : pos{position}, lum{std::move(luminosity)} {}
        light_src(const PosT &x, const PosT &y, const PosT &z, LumT &&luminosity) :
        pos{x, y, z}, lum{std::move(luminosity)} {}
        light_src(PosT &&x, PosT &&y, PosT &&z, LumT &&luminosity) :
        pos{std::move(x), std::move(y), std::move(z)}, lum{std::move(luminosity)} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(const U &x, const U &y, const U &z, LumT &&luminosity) : pos{x, y, z}, lum{std::move(luminosity)} {}
        template <isNumWrapper U> requires isConvertible<U, PosT>
        light_src(const U &&x, const U &&y, const U &&z, LumT &&luminosity) :
        pos{x, y, z}, lum{std::move(luminosity)} {}
        vector3D<PosT> position() const {
            return pos;
        }
        LumT luminosity() const {
            return lum;
        }
        ray<PosT, DirT, LenT> cast_ray(const vector3D<PosT> &end_point) const requires isConvertible<long double, LenT>{
            return {pos, end_point - pos}; // length zero as its intersection with bodies has not yet been calculated
        }
        ray<PosT, DirT, LenT> cast_ray(const vector3D<PosT> &&end_point) const requires isConvertible<long double,LenT>{
            return {pos, end_point - pos};
        }
        LumT flux_at_point(const vector3D<PosT> &pnt) const requires isConvertible<long double,LenT> {
            auto &&d = (pnt - pos).magnitude();
            return lum/(/* 4*__PI__* */d*d); // 4PI removed since all brightnesses are relative anyway
        }
        LumT flux_at_point(const vector3D<PosT> &&pnt) const requires isConvertible<long double,LenT> {
            auto &&d = (pnt - pos).magnitude();
            return lum/(/* 4*__PI__* */d*d);
        }
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU, isNumWrapper LumU>
        friend std::ostream &operator<<(std::ostream &os, const light_src<PosU, DirU, LenU, LumU> &src);
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                  isNumWrapper, bool, uint64_t>
        friend class astro_scene;
    };
    template <isNumWrapper PosT = long double, isNumWrapper DirT = long double, isNumWrapper DistT = long double>
    class camera {
        /* The "camera" class is one which represents a pinhole camera in 3D space with a certain position, orientation
         * (direction it's pointing in) and rotation. A camera is essentially a box with a certain width and height
         * (which are directly proportional to the real image's dimensions) that are "normalised" in order for the
         * shorter of the two to always equal 2. The length of the camera "box" is the image distance. The pinhole lies
         * on one of the sides that are perpendicular to the lengthwise dimension, and the location of the pinhole in
         * simulated 3D space is taken as the location of the camera. On the face opposite the pinhole (aperture) lies
         * the camera's "receptor", which is where the image in simulated 3D space is projected onto. This is done by
         * sending a ray out from every single "pixel" on the receptor through the pinhole, and calculating what it
         * intersects with - although the intersection calculation is not done in the camera class (see astro_scene). */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        static inline const vector3D<DirT> zero{}; // as this is repeatedly used
        static inline const vector3D<DirT> up{DirT{0}, DirT{0}, DirT{1}}; // same with all these
        static inline const vector3D<DirT> down{DirT{0}, DirT{0}, DirT{-1}};
        static inline const vector3D<DirT> right{DirT{1}, DirT{0}, DirT{0}};
        static inline const vector3D<DirT> forwards{DirT{0}, DirT{1}, DirT{0}};
        image_dimensions dims{1000, 1000}; // the real image dimensions in pixels
        vector3D<PosT> pos; // position of the camera
        vector3D<DirT> dir{0, 0, 1}; // direction the camera is pointing in (doesn't have to be a unit vector)
        long double rot = 0; /* this specifies the roll of the camera (see pitch, yaw & roll) in terms of an angle in
        * radians - in the case of the camera's direction being perpendicular to the x-y plane, this angle is with
        * respect to the positive y-axis */
        DistT dist = DistT{2}; /* distance from the camera's "receptor" - i.e., the image projected onto its
        * back - to the camera's aperture (the pinhole) -> greater image distance = greater zoom */
        vector3D<DirT> perp; // vector around which rays are rotated to be in the correct orientation
        long double angle_with_down = 0; // the angle the camera direction makes with the negative z-axis
        long double rotation_correction = 0; // a correction factor needed to due to rays being rotated
        bool points_up = true; // true if the camera is pointing directly up (direction parallel to positive z-axis)
        bool points_down = false;// true if the camera is pointing directly down (direction parallel to negative z-axis)
        long double w_factor = 1; // these factors dictate how the x- and y-components of the ray vectors should be
        long double h_factor = 1; // ... scaled based on the true dimensions of the image
        bool invert_ray = true; /* this parameter dictates whether a ray emanating from the pinhole has its y-component
        * flipped or not - this, essentially, avoids the image having to be flipped later - although for images that
        * have their origin in the top-left corner (unlike BMPs), this should be set to false */
        void calc_all() {
            check_non_zero();
            if (dir.x == DirT{0} && dir.y == DirT{0}) {
                perp.set_x(DirT{1});
                perp.set_y(DirT{0});
                perp.set_z(DirT{0});
                rotation_correction = 0;
                if ((points_down = !(points_up = dir.z > DirT{0})))
                    angle_with_down = 0;
                else
                    angle_with_down = PI;
                goto end;
            }
            perp = vec_ops::cross(down, dir);
            angle_with_down = vec_ops::angle_between(down, dir);
            if (dir.x < 0)
                rotation_correction = vec_ops::angle_between(forwards, dir.xy_projection());
            else if (dir.x > 0)
                rotation_correction = -vec_ops::angle_between(forwards, dir.xy_projection());
            else if (dir.y < 0)
                rotation_correction = PI;
            end:
            if (dims.x >= dims.y) {
                w_factor = ((long double) dims.x)/dims.y;
                h_factor = 1;
            }
            else {
                w_factor = 1;
                h_factor = ((long double) dims.y)/dims.x;
            }
            if (!invert_ray)
                h_factor = -h_factor;
        }
        void check_non_zero() { // the camera cannot, however, not be pointing in any direction
            if (dir == zero)
                throw no_direction("The camera's direction is a zero vector. The camera must point in a direction.\n");
            if (dist == DistT{0})
                throw zero_image_distance("A camera object's image distance cannot be zero, or else simulated bodies"
                                          "directly\nin front of the camera would appear to be infinitely far away.\n");
        }
    public:
        camera() = default; // for a camera at the origin pointing up (taking the z-axis as up)
        explicit camera(const image_dimensions &image_dim) : dims{image_dim} {calc_all();}
        explicit camera(const vector3D<PosT> &position) : pos{position} {calc_all();} // for a camera pointing up
        explicit camera(const vector3D<PosT> &&position) : pos{std::move(position)} {calc_all();}
        camera(const vector3D<PosT> &position, const vector3D<DirT> &direction) :
        pos{position}, dir{direction} {calc_all();}
        camera(const vector3D<PosT> &&position, const vector3D<DirT> &&direction) :
        pos{std::move(position)}, dir{std::move(direction)} {calc_all();}
        camera(const vector3D<PosT> &position, const vector3D<DirT> &direction, const long double &rotation,
               const DistT &image_distance) : pos{position}, dir{direction}, rot{rotation}, dist{image_distance}
               {calc_all();}
        camera(const vector3D<PosT> &&position, const vector3D<DirT> &&direction, const long double &&rotation,
               const DistT &&image_distance) :
                pos{std::move(position)}, dir{std::move(direction)}, rot{rotation}, dist{std::move(image_distance)}
        {calc_all();}
        camera(const vector3D<PosT> &position, const image_dimensions &image_dim) : pos{position}, dims{image_dim}
        {calc_all();}
        camera(const vector3D<PosT> &&position, const image_dimensions &&image_dim) :
                pos{std::move(position)}, dims{image_dim} {calc_all();}
        camera(const vector3D<PosT> &position, const vector3D<DirT> &direction, const image_dimensions &image_dim) :
                pos{position}, dir{direction}, dims{image_dim} {calc_all();}
        camera(const vector3D<PosT> &&position, const vector3D<DirT> &&direction, const image_dimensions &&image_dim) :
                pos{std::move(position)}, dir{std::move(direction)}, dims{image_dim} {calc_all();}
        camera(const vector3D<PosT> &position, const vector3D<DirT> &direction, const long double &rotation,
               const DistT &image_distance, const image_dimensions &image_dim) : pos{position}, dir{direction},
               rot{rotation}, dist{image_distance}, dims{image_dim} {calc_all();}
        camera(const vector3D<PosT> &&position, const vector3D<DirT> &&direction, const long double &&rotation,
               const DistT &&image_distance, const image_dimensions &&image_dim) : pos{std::move(position)},
               dir{std::move(direction)}, rot{rotation}, dist{std::move(image_distance)}, dims{image_dim}
        {calc_all();}
        const vector3D<PosT> &position() const noexcept {
            return pos;
        }
        const vector3D<DirT> &direction() const noexcept {
            return dir;
        }
        const long double &rotation() const noexcept {
            return rot;
        }
        const DistT &image_distance() const noexcept {
            return dist;
        }
        const image_dimensions &img_dim() const noexcept {
            return dims;
        }
        camera<PosT, DirT, DistT> &set_position(const vector3D<PosT> &new_position) {
            pos = new_position;
            this->calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_direction(const vector3D<DirT> &new_direction) {
            dir = new_direction;
            this->calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_rotation(const long double &new_rotation) {
            rot = new_rotation;
            return *this;
        }
        camera<PosT, DirT, DistT> &set_image_distance(const DistT &new_image_distance) {
            dist = new_image_distance;
            return *this;
        }
        camera<PosT, DirT, DistT> &set_image_dimensions(const image_dimensions &new_dimensions) {
            dims = new_dimensions;
            this->calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_position(vector3D<PosT> &&new_position) {
            pos = std::move(new_position);
            this->calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_direction(vector3D<DirT> &&new_direction) {
            dir = std::move(new_direction);
            this->calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_rotation(const long double &&new_rotation) {
            rot = new_rotation;
            return *this;
        }
        camera<PosT, DirT, DistT> &set_image_distance(DistT &&new_image_distance) {
            dist = std::move(new_image_distance);
            return *this;
        }
        camera<PosT, DirT, DistT> &set_image_dimensions(const image_dimensions &&new_dimensions) {
            dims = new_dimensions;
            this->calc_all();
            return *this;
        }
        void recalculate() { // useless method, pretty much
            this->calc_all();
        }
        camera<PosT, DirT, DistT> &invert_rays(bool invert) {
            invert_ray = invert;
            h_factor = -h_factor;
            return *this;
        }
        long double fovh_rad() const noexcept { // horizontal field of view in radians
            return dims.x <= dims.y ? 2*atanl(1.0l/dist) : 2*atanl(dims.x/(dims.y*dist));
        }
        long double fovv_rad() const noexcept { // vertical field of view in radians
            return dims.y <= dims.x ? 2*atanl(1.0l/dist) : 2*atanl(dims.y/(dims.x*dist));
        }
        long double fovd_rad() const noexcept { // diagonal field of view in radians
            long double x_ratio;
            long double y_ratio;
            if (dims.x > dims.y) {
                x_ratio = ((long double) dims.x)/dims.y;
                y_ratio = 1;
            }
            else {
                x_ratio = 1;
                y_ratio = ((long double) dims.y)/dims.x;
            }
            return 2*atanl(sqrtl(x_ratio*x_ratio + y_ratio*y_ratio)/dist);
        }
        long double fovh_deg() const noexcept { // horizontal field of view in degrees
            return rad_to_deg(fovh_rad());
        }
        long double fovv_deg() const noexcept { // vertical field of view in degrees
            return rad_to_deg(fovv_rad());
        }
        long double fovd_deg() const noexcept { // diagonal field of view in degrees
            return rad_to_deg(fovd_rad());
        }
        template <isNumWrapper LenT = long double>
        ray<PosT, DirT, LenT> ray_from_pixel(unsigned int x, unsigned int y) const {
            if (x >= dims.x || y >= dims.y)
                throw std::out_of_range("The pixel attempting to be accessed is not within the image dimensions.\n");
            ray<PosT, DirT, LenT> // ray starts out as if it had been shot out of the camera facing directly down
            ry{this->pos,
                ((((long double) x)/dims.x)*2 - 1)*w_factor, ((((long double) y)/dims.y)*2 - 1)*h_factor, -dist};
            return ry.rodrigues_rotate(perp, angle_with_down).rodrigues_rotate(dir, -rotation_correction - rot);
        } // ^^^ rotate to camera's direction and rotate around camera's z-axis appropriately
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                  isNumWrapper, bool, uint64_t>
        friend class astro_scene;
    };
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double,
              isNumWrapper PosT = long double, isNumWrapper DirT = long double, isNumWrapper LenT = long double,
              isNumWrapper LumT = long double, uint64_t recFreq = 0> requires isConvertible<long double, PosT>
    class star : public body<M, R, T, recFreq> {
        using bod_t = body<M, R, T, recFreq>;
        using src_t = light_src<PosT, DirT, LenT, LumT>;
        using us = unsigned short;
        LumT lum{1};
        bool pt_src = true; // specifies whether the star is approximated as a point source of light in astro_scene
        std::vector<src_t> sources;
        void create_sources(unsigned short num) {
            sources.clear();
            if (num <= 1) {
                sources.push_back({bod_t::curr_pos, lum}); // can't have zero sources - in that case, just use a body
                pt_src = true;
                return;
            }
            if (num % 2) // num should be even
                --num;
            // std::random_device r;
            srand(time(nullptr));
            std::mt19937 mersenne{(unsigned int) rand()};
            std::uniform_real_distribution<long double> dist{-1, 1};
            vector3D<long double> source_pos;
            LumT src_lum = lum/num;
            while (num --> 0) {
                source_pos.set_x(dist(mersenne));
                source_pos.set_y(dist(mersenne));
                source_pos.set_z(dist(mersenne));
                source_pos.set_length(bod_t::radius);
                sources.push_back({bod_t::curr_pos + source_pos, src_lum});
            }
            pt_src = false;
        }
    public:
        explicit star(us num_light_sources = 1) {create_sources(num_light_sources);}
        star(M &&star_mass, R &&star_radius, LumT &&luminosity, us num_light_sources = 1) :
        bod_t{std::move(star_mass), std::move(star_radius)}, lum{std::move(luminosity)}
        {create_sources(num_light_sources);}
        star(const M &star_mass, const R &star_radius, const LumT &luminosity, us num_light_sources = 1) :
        bod_t{star_mass, star_radius}, lum{luminosity} {create_sources(num_light_sources);}
        star(M &&star_mass, R &&star_radius, vector3D<T> &&pos, vector3D<T> &&vel, LumT &&luminosity,
             us num_light_sources = 1):
        bod_t{std::move(star_mass), std::move(star_radius), std::move(pos), std::move(vel)}, lum{std::move(luminosity)}
        {create_sources(num_light_sources);}
        star(const M &star_mass, const R &star_radius, const vector3D<T> &pos, const vector3D<T> &vel,
             const LumT &luminosity, us num_light_sources = 1) :
        bod_t{star_mass, star_radius, pos, vel}, lum{luminosity} {create_sources(num_light_sources);}
        template <isConvertible<M> m, isConvertible<R> r, isConvertible<T> t, uint64_t rF>
        explicit star(const body<m, r, t, rF> &other) : bod_t{other} {}
        template <uint64_t rF>
        explicit star(const body<M, R, T, rF> &other) : bod_t{other} {}
        template <uint64_t rF>
        explicit star(body<M, R, T, rF> &&other) noexcept : bod_t{std::move(other)} {}
        const std::vector<src_t> &light_sources() const {
            return sources;
        }
        bool is_pnt_src() const noexcept {
            return pt_src;
        }
        unsigned short num_sources() {
            return sources.size();
        }
        void recreate_sources(unsigned short new_number = 1) {
            create_sources(new_number);
        }
        vector3D<DirT> normal_at_point(const vector3D<PosT> &pnt) const {
            if (pnt == bod_t::curr_pos)
                return {};
            return (pnt - bod_t::curr_pos).normalise();
        }
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                  isNumWrapper, bool, uint64_t>
        friend class astro_scene;
    };
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double,
              isNumWrapper PosT = long double, isNumWrapper DirT = long double, isNumWrapper DistT = long double,
              isNumWrapper LenT = long double, isNumWrapper LumT = long double, bool modulatedBrightness = true,
              uint64_t recFreq = 0> // stars and bodies must have the same recFreq, or else a body is not a star's parent
    class astro_scene : public bmp {
    private:
        using bod_t = body<M, R, T, recFreq>;
        using ray_t = ray<PosT, DirT, LenT>;
        using src_t = light_src<PosT, DirT, LenT, LumT>;
        using cam_t = camera<PosT, DirT, DistT>;
        using star_t = star<M, R, T, PosT, DirT, LenT, LumT, recFreq>;
        image_dimensions dims{bmp::width, bmp::height};
        unsigned int num_decor_stars = 4*logl(bmp::width*bmp::height);
        long double star_radius = sqrtl(bmp::width*bmp::width)/1000.0l; // in pixels
        std::vector<point> star_points;
        /* indeed, all the std::maps below lead to a larger memory footprint, but, they also allow certain algorithms
         * to be implemented faster */
        std::map<const uint64_t, const bod_t*> tot_map;
        std::map<const uint64_t, const bod_t*> bodies;
        std::map<const uint64_t, const star_t*> stars;
        std::map<const uint64_t, color> body_clrs;
        std::map<const uint64_t, color> star_clrs;
        cam_t cam{dims};
        const cam_t *pcam = &cam;
        bool def_cam = true;
        bool rendered = false;
        bool cleared = false;
        LumT **fdata = nullptr;
        /* the five variables below are to support concurrency in the ray-tracing algorithm - they are only needed, in
         * reality, to ensure brightnesses of pixels are calculated correctly - since the absolute brightness of any
         * pixel depends on the brightness of the brightest pixel, all brightnesses (for all pixels) must be calculated
         * first, hence the need to synchronise the different threads */
        std::mutex max_mutex;
        std::condition_variable cv;
        std::latch *thread_latch;
        LumT max_flux;
        std::vector<LumT> thread_max_fluxes;
        bool have_max = false;
        unsigned int tot_threads;
        long double min_clr_brightness_b = 128; // minimum average BGR value for bodies (in default case)
        long double min_clr_brightness_s = 240; // minimum average BGR value for stars (in default case)
        std::function<color()> col_gen_b = [this](){ // default colour generator for bodies
            std::mt19937 mersenne{(unsigned int) rand()};
            std::uniform_int_distribution dist{0, 255};
            color ret{};
            do {
                ret.b = dist(mersenne);
                ret.g = dist(mersenne);
                ret.r = dist(mersenne);
            } while (avg_bgr(ret) < min_clr_brightness_b);
            return ret;
        }; // std::function objects are being used as these are capturing lambdas (cannot assign to function pointers)
        std::function<color()> col_gen_s = [this](){ // default colour generator for bodies
            std::mt19937 mersenne{(unsigned int) rand()};
            std::uniform_int_distribution dist{0, 255};
            color ret{};
            do {
                ret.b = dist(mersenne);
                ret.g = dist(mersenne);
                ret.r = dist(mersenne);
            } while (avg_bgr(ret) < min_clr_brightness_s);
            return ret;
        };
        void populate_tot_map() {
            tot_map.clear();
            tot_map.insert(bodies.begin(), bodies.end());
            tot_map.insert(stars.begin(), stars.end());
        }
        void create_decor_stars(bool reset) {
            if (!num_decor_stars)
                return;
            sc_clr = colors::white;
            circle c;
            c.set_radius(star_radius);
            std::random_device r;
            if (reset) {
                point p;
                star_points.clear();
                for (unsigned int i = 0; i < num_decor_stars; ++i) {
                    p.x = r() % bmp::width;
                    p.y = r() % bmp::height;
                    c.set_pos(p);
                    bmp::draw_circle(c);
                    star_points.push_back(p);
                }
                return;
            }
            for (const point &p : star_points) {
                c.set_pos(p);
                bmp::draw_circle(c);
            }
        }
        void alloc() {
            if constexpr (!modulatedBrightness)
                return;
            dealloc();
            fdata = new LumT*[bmp::height];
            LumT **ptr = fdata;
            for (unsigned int j = 0; j < bmp::height; ++j, ++ptr)
                *ptr = new LumT[bmp::width];
            if (bmp::data == nullptr)
                bmp::reallocate();
        }
        void dealloc() {
            if (fdata == nullptr)
                return;
            LumT **ptr = fdata;
            for (unsigned int j = 0; j < bmp::height; ++j, ++ptr)
                delete [] *ptr;
            delete [] fdata;
            fdata = nullptr;
        }
        void put_bodies(const std::vector<bod_t*> &vec, bool erase = true) {
            if (erase) {
                bodies.clear();
                stars.clear();
            }
            const star_t *star_ptr = nullptr;
            for (const bod_t *const &bod_ptr : vec) {
                if (bod_ptr == nullptr) // no need to throw an exception, simply ignored
                    continue;
                star_ptr = dynamic_cast<const star_t*>(bod_ptr); // no bad_cast exception thrown with pointers
                if (star_ptr == nullptr) {
                    bodies.emplace(bod_ptr->id, bod_ptr);
                    body_clrs.emplace(bod_ptr->id, col_gen_b());
                }
                else {
                    stars.emplace(star_ptr->id, star_ptr); // all bodies have a unique id, hence why a map is used
                    star_clrs.emplace(bod_ptr->id, col_gen_s());
                }
            }
        }
        void set_default_cam() { // needs work - this is where the default camera will be set
            cam.dims = this->dims;
            // get COM of all bodies
            // get vector perpendicular to the plane the bodies lie in
            // do this by calculating the magnitude of the cross product between pairs of bodies, and then taking the
            // ... average unit vector of all the cross products created
            // place camera on this vector at just enough distance to fit all bodies into FOV (FOV set to arb. value)
            // if (fcam != nullptr)
            //     fcam->dims = this->dims;
        }
        void check_cam() const {
            if (cam.dims != this->dims)
                throw camera_dimensions_error("The dimensions of the camera supplied do not equal that of this "
                                              "astro_scene object.\n");
        }
        void render_subp(unsigned int start_x, unsigned int end_x, unsigned int start_y, unsigned int end_y) {
            typename std::vector<const star_t*>::size_type num_s = stars.size();
            unsigned int x; // inner loop counter
            unsigned int flux_counter;
            color **row_c = bmp::data + start_y;
            LumT **row_f = fdata + start_y;
            color *pix_c; // a single pixel's colour
            LumT *pix_f; // a single pixel's cumulative flux - later converted to brightness - only if modBright true
            const star_t *sptr;
            LumT thread_max_flux{0}; // maximum flux falling on one of the visible points (in the thread) of any object
            typename std::map<const uint64_t, const star_t*>::iterator it;
            auto end_it = stars.end(); // thank goodness for auto ;)
            for (unsigned int y = start_y; y < end_y; ++y, ++row_c, ++row_f) {
                pix_c = *row_c + start_x;
                pix_f = *row_f + start_x;
                for (x = start_x; x < end_x; ++x, ++pix_c, ++pix_f) {
                    ray_t &&cam_ray = pcam->ray_from_pixel(x, y);
                    if (!cam_ray.template intersects<M, R, T, recFreq>(tot_map)) {
                        *pix_c = colors::black; // in case of no intersection, pixel is the black of space
                        *pix_f = LumT{0}; // not truly necessary, since black multiplied by anything will still be black
                        continue;
                    }
                    if (stars.contains(cam_ray.ibod_id)) {
                        *pix_c = star_clrs[cam_ray.ibod_id];
                        *pix_f = LumT{-1}; // only pixels that lie on stars will have a negative cumulative lum.
                        continue;
                    }
                    // next comes the case of the ray having intersected with a body
                    *pix_c = body_clrs[cam_ray.ibod_id]; // this was the easy part ;)
                    *pix_f = LumT{0}; // cumulative flux must start out as zero (W m^-2)
                    vector3D<PosT> &&end_point = cam_ray.end_point();
                    for (it = stars.begin(); it != end_it; ++it) {
                        sptr = it->second;
                        if (sptr->sources.size() == 1) {
                            ray_t &&src_to_ep = sptr->sources[0].cast_ray(end_point);
                            if (!src_to_ep.template intersects<M, R, T, recFreq>(stars, sptr->id) &&
                                !src_to_ep.template intersects<M, R, T, recFreq>(bodies)) // in case of f.p. error
                                continue;
                            if (!src_to_ep.reaches(bodies[cam_ray.ibod_id], end_point)) // point is behind light source
                                continue;
                            *pix_f -= sptr->sources[0].flux_at_point(end_point)*
                                ((src_to_ep.normalise()*(src_to_ep.normal.normalise())));
                            /* value subtracted since the normal and light ray always have an obtuse angle between */
                        }
                        else {
                            for (const src_t &src : sptr->sources) {
                                ray_t &&src_to_ep = src.cast_ray(end_point);
                                if (sptr->normal_at_point(src.pos)*src_to_ep <= 0)
                                    continue;
                                if (!src_to_ep.template intersects<M, R, T, recFreq>(stars, sptr->id) &&
                                    !src_to_ep.template intersects<M, R, T, recFreq>(bodies))
                                    continue;
                                if (!src_to_ep.reaches(bodies[cam_ray.ibod_id], end_point))
                                    continue;
                                *pix_f -= src.flux_at_point(end_point)*
                                ((src_to_ep.normalise()*(src_to_ep.normal.normalise())));
                            }
                        }
                    }
                    if (*pix_f > thread_max_flux)
                        thread_max_flux = *pix_f;
                }
            }
            /* now the entire arrays of colours and flux quantities have to be gone through again to multiply each
             * colour by the normalised flux values (i.e., each flux value is divided by the maximum flux value
             * (to yield a value between 0-1), and then multiplied by the corresponding colour at the given pixel) */
            thread_latch->count_down(); // reduces the latch's internal counter by 1
            thread_latch->wait(); // waits until the internal counter is 0 (so once all threads have reached this point)
            max_mutex.lock(); // ensure that only one thread is accessing the std::vector at a time
            thread_max_fluxes.push_back(thread_max_flux); // each thread adds its maximum brightness to the std::vector
            max_mutex.unlock();
            cv.notify_all(); // notify the running set_max_flux function to continue with the max. calculation
            std::unique_lock<std::mutex> lock{max_mutex}; // std::condition_variable requires a lock
            cv.wait(lock, [this](){return have_max;}); // wait for the set_max_flux function having calculated the max
            row_c = bmp::data + start_y;
            row_f = fdata + start_y;
            for (unsigned int y = start_y; y < end_y; ++y, ++row_c, ++row_f) {
                pix_c = *row_c + start_x;
                pix_f = *row_f + start_x;
                for (x = start_x; x < end_x; ++x, ++pix_c, ++pix_f)
                    if (*pix_f != -1)
                        *pix_c = (*pix_f/max_flux)*(*pix_c);
            }
        }
        void set_max_flux() {
            std::unique_lock<std::mutex> lock{max_mutex};
            cv.wait(lock, [this](){return thread_max_fluxes.size() == tot_threads;});
            max_flux = LumT{0};
            for (const LumT &t_max : thread_max_fluxes) {
                if (t_max > max_flux)
                    max_flux = t_max;
            }
            have_max = true;
            lock.unlock();
            cv.notify_all();
        }
    public:
        astro_scene() :
        bmp{std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp")}
        {alloc(); srand(time(nullptr));}


        explicit astro_scene(const char *source_bmp) :
        bmp{source_bmp, std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp")}
        {alloc(); srand(time(nullptr));}


        astro_scene(unsigned int bmp_width, unsigned int bmp_height) :
        bmp{std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp"),
            colors::black, bmp_width, bmp_height} {alloc(); srand(time(nullptr));}


        astro_scene(unsigned int bmp_width, unsigned int bmp_height, const cam_t &c,
                    const std::vector<bod_t*> &bods) :
                bmp{std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp"),
                    colors::black, bmp_width, bmp_height}, cam{c}
                    {check_cam(); put_bodies(bods); alloc(); srand(time(nullptr)); def_cam = false;}


        astro_scene(unsigned int bmp_width, unsigned int bmp_height, const cam_t &c,
                    const std::vector<bod_t*> &bods, color (*bdy_clr_gen)()) :
                bmp{std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp"),
                    colors::black, bmp_width, bmp_height}, cam{c}, col_gen_b{bdy_clr_gen}
        {check_cam(); put_bodies(bods); alloc(); srand(time(nullptr)); def_cam = false;}


        astro_scene(unsigned int bmp_width, unsigned int bmp_height, const cam_t &c,
                    const std::vector<bod_t*> &bods, color (*bdy_clr_gen)(), color (*star_clr_gen)()) :
                bmp{std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp"),
                    colors::black, bmp_width, bmp_height}, cam{c}, col_gen_b{bdy_clr_gen},
                    col_gen_s{star_clr_gen}
        {check_cam(); put_bodies(bods); alloc(); srand(time(nullptr)); def_cam = false;}


        astro_scene(const String &bmp_path, unsigned int bmp_width, unsigned int bmp_height) :
                bmp{bmp_path, colors::black, bmp_width, bmp_height} {alloc(); srand(time(nullptr));}


        astro_scene(String &&bmp_path, unsigned int bmp_width, unsigned int bmp_height) :
                bmp{std::move(bmp_path), colors::black, bmp_width, bmp_height}
                {alloc(); srand(time(nullptr));}


        astro_scene(const astro_scene &other) : bmp{other} {} // not = default, as I will add stuff in body later


        astro_scene(astro_scene &&other) noexcept : bmp{std::move(other)} {}


        void set_camera(const cam_t &c) { // if an outside camera is being followed, this method does not change that
            cam = c;
            check_cam();
            def_cam = false;
        }
        void set_camera(cam_t &&c) {
            cam = std::move(c);
            check_cam();
            def_cam = false;
        }
        bool set_body_default_min_clr_brightness(const long double &val) noexcept {
            /* sets the minimum brightness for the randomly generated colours for the bodies - has no effect if colours
             * are manually specified by the user */
            if (val < 0 || val > 255)
                return false;
            min_clr_brightness_b = val;
            return true;
        }
        bool set_body_default_min_clr_brightness(long double &&val) noexcept {
            return set_body_default_min_clr_brightness(val);
        }
        bool set_star_default_min_clr_brightness(const long double &val) noexcept {
            /* same as above but for stars */
            if (val < 0 || val > 255)
                return false;
            min_clr_brightness_s = val;
            return true;
        }
        bool set_star_default_min_clr_brightness(long double &&val) noexcept {
            return set_star_default_min_clr_brightness(val);
        }
        template <colorGenerator colGenT>
        void set_body_col_gen(const colGenT &func) noexcept {
            col_gen_b = func;
        }
        template <colorGenerator colGenT>
        void set_star_col_gen(const colGenT &func) noexcept {
            col_gen_s = func;
        }
        const std::function<color()> &get_body_col_gen() const noexcept {
            return col_gen_b;
        }
        const std::function<color()> &get_star_col_gen() const noexcept {
            return col_gen_s;
        }
        color get_clr_by_id(const uint64_t &id) {
            if (body_clrs.contains(id))
                return body_clrs[id];
            if (star_clrs.contains(id))
                return star_clrs[id];
            throw invalid_id_error(id);
        }
        typename std::map<const uint64_t, const bod_t*>::size_type num_bodies() {
            return bodies.size();
        }
        typename std::map<const uint64_t, const star_t*>::size_type num_stars() {
            return stars.size();
        }
        bool set_body_clr(const uint64_t &id, color col) {
            if (!bodies.contains(id))
                return false;
            body_clrs[id] = col;
            return true;
        }
        // bool set_body_clr(unsigned long long &&id, color col) {
        //     return set_body_clr(id, col);
        // }
        bool set_star_clr(const uint64_t &id, color col) {
            if (!stars.contains(id))
                return false;
            star_clrs[id] = col;
            return true;
        }
        bool set_star_clr(uint64_t &&id, color col) {
            return set_star_clr(id, col);
        }
        bool set_num_decor_stars(unsigned int num) {
            if (num*star_radius*star_radius >= bmp::width*bmp::height)
                return false;
            num_decor_stars = num;
            return true;
        }
        bool set_decor_star_rad(long double rad) {
            if (num_decor_stars*rad*rad >= bmp::width*bmp::height || rad < 0)
                return false;
            star_radius = rad;
            return true;
        }
        void use_default_camera() noexcept { // this automatically unfollows an outside camera if it was being followed
            def_cam = true;
            *pcam = cam;
        }
        bool follow_camera(const cam_t *c) noexcept {
            if (c == nullptr)
                return false;
            pcam = c; // as the camera pointed to by pcam can be altered outside, dims are only checked at last moment
            return true;
        }
        bool unfollow_camera() noexcept { // returns true if outside camera was being followed, else false
            return pcam != &cam && (pcam = &cam);
        }
        void generate_body_clrs() {
            std::for_each(bodies.begin(), bodies.end(), [this](const std::pair<const uint64_t,
                                                               const bod_t*> &p) {
                body_clrs[p.first] = col_gen_b(); // unfortunately, cannot iterate over keys only
            });
        }
        void generate_star_clrs() {
            std::for_each(stars.begin(), stars.end(), [this](const std::pair<const uint64_t,
                                                             const star_t*> &p){
                star_clrs[p.first] = col_gen_s();
            });
        }
        void generate_clrs() {
            this->generate_body_clrs();
            this->generate_star_clrs();
        }
        void add_body(const bod_t &bod, color col) {
            try {
                const star_t &s = dynamic_cast<const star_t&>(bod);
                stars.emplace(s.id, &s);
                star_clrs.emplace(s.id, col);
            } catch (const std::bad_cast &e) {
                bodies.emplace(bod.id, &bod);
                body_clrs.emplace(bod.id, col);
            }
        }
        void add_body(const bod_t &bod) {
            this->add_body(bod, col_gen_b());
        }
        void add_star(const star_t &s, color col) {
            stars.emplace(s.id, &s);
            star_clrs.emplace(s.id, col);
        }
        void add_star(const star_t &s) {
            stars.emplace(s.id, &s);
            star_clrs.emplace(s.id, col_gen_s());
        }
        void add_bodies(const std::vector<bod_t*> &new_bodies) {
            put_bodies(new_bodies, false);
        }
        void add_bodies(const std::vector<bod_t*> &&new_bodies) { // no way to move it, so const
            put_bodies(new_bodies, false);
        }
        void add_stars(const std::vector<star_t*> &new_stars) {
            for (const star_t *const &s : new_stars) {
                if (s == nullptr)
                    continue;
                stars.emplace(s->id, s);
                star_clrs.emplace(s->id, col_gen_s());
            }
        }
        void add_stars(const std::vector<star_t*> &&new_stars) {//with these r-value overloads, init. syntax may be used
            this->add_stars(new_stars);
        }
        template <bool progress, bool merge_if_overlapping, int coll, uint64_t mF, uint64_t fF, bool bF>
        void add_system(const system<M, R, T, progress, merge_if_overlapping, coll, mF, fF, bF> &sy) {
            for (const bod_t &b : sy.bods) {
                bodies.emplace(b.id, &b);
                body_clrs.emplace(b.id, col_gen_b());
            }
        }
        void assign_bodies(const std::vector<bod_t*> &new_bodies) {
            put_bodies(new_bodies, true);
        }
        void assign_stars(const std::vector<star_t*> &new_stars) {
            stars.clear();
            for (const star_t *const &s : new_stars) {
                if (s == nullptr)
                    continue;
                stars.emplace(s->id, s);
                star_clrs.emplace(s->id, col_gen_s());
            }
        }
        /* only makes sense to call this method if not following outside camera */
        bool correct_camera_dimensions() noexcept {
            cam.dims = this->dims;
            return pcam == &cam; // returns true only if internal camera is being used
        }
        bool render(bool reset_star_positions = true) {
            if (cleared)
                this->alloc();
            if (rendered) {
                bmp::sc_clr = colors::black;
                bmp::fill_bg();
            }
            if (pcam == &cam)
                if (def_cam)
                    this->set_default_cam();
            if (pcam->dims != this->dims)
                /* as mentioned further above, pcam is only checked as it might have been pointing to an outside cam. */
                /* an exception is raised instead of just silently correcting it because 1. it would be good for the
                 * caller to know that the camera they are providing is of incorrect dimensions and 2. this class
                 * guarantees not to modify a followed camera (hence why pcam points to a const camera) */
                throw camera_dimensions_error("The camera dimensions do not equal astro_scene dimensions.\n");
            rendered = true;
            this->create_decor_stars(reset_star_positions);
            this->populate_tot_map();
            have_max = false;
            thread_max_fluxes.clear();
            /* in case of any error below, single-threaded execution will be reverted to */
            if ((tot_threads = std::thread::hardware_concurrency()) == 0 || tot_threads == 1) {
                tot_threads = 1;
                std::latch once{1};
                thread_latch = &once;
                std::thread max_calculator{&astro_scene::set_max_flux, this};
                this->render_subp(0, bmp::width, 0, bmp::height);
                max_calculator.join();
                return false; // method returns false in case no multithreading was able to be performed
            }
            std::thread max_calculator{&astro_scene::set_max_flux, this};
            /* my approach to the multithreading is that it should be favoured to split the image up into rows, rather
             * than columns, to minimise the jumps around in memory (since rows are stored contiguously, whereas
             * columns are not) - if the number of threads to execute is less than the number of rows (always the case
             * on my home computer), then the image is only split up rows, which are then dealt with separately in each
             * thread - however, if threads > height (could only be the case for a really good processing unit), then
             * the image starts to be split up vertically as well */
            std::latch t_latch{tot_threads};
            thread_latch = &t_latch;
            std::vector<std::pair<long double, long double>> ranges;
            std::vector<std::thread> threads;
            long double bottom_or_left = 0; // most likely bottom
            long double top_or_right = 0; // ... and top
            unsigned int counter = 1;
            if (tot_threads <= bmp::height) { // on most home computers
                ranges.reserve(tot_threads);
                long double y_interval;
                y_interval = ((long double) bmp::height)/tot_threads;
                while (top_or_right <= bmp::height) {
                    top_or_right = y_interval*counter++;
                    ranges.emplace_back(bottom_or_left, top_or_right);
                    bottom_or_left = top_or_right;
                }
                if (ranges.size() > tot_threads) // in case f.p. error above caused the loop to execute an extra time
                    ranges.pop_back();
                ranges.back().second = bmp::height; // again, in case of f.p. error
                try {
                    for (const auto &[y_start, y_end] : ranges)
                        threads.emplace_back(&astro_scene::render_subp, this, 0, bmp::width, y_start, y_end);
                    for (std::thread &t : threads)
                        t.join();
                    max_calculator.join();
                } catch (const std::system_error &e) {
                    tot_threads = 1;
                    std::latch once{1};
                    thread_latch = &once;
                    this->render_subp(0, bmp::width, 0, bmp::height);
                    max_calculator.join();
                    return false;
                }
                return true;
            }
            // for a supercomputer or excellent GPU (such as in gaming setup)
            unsigned int end_beg_row; // less than
            unsigned int num_beg_x_intervals;
            unsigned int num_end_x_intervals;
            long double beg_x_interval;
            if ((end_beg_row = tot_threads % bmp::height)) {
                num_end_x_intervals = tot_threads/bmp::height;
                num_beg_x_intervals = num_end_x_intervals + 1;
            }
            else {
                num_beg_x_intervals = tot_threads/bmp::height;
                end_beg_row = bmp::height;
            }
            ranges.reserve(num_beg_x_intervals);
            beg_x_interval = bmp::width/num_beg_x_intervals;
            bottom_or_left = 0;
            top_or_right = 0;
            while (top_or_right <= bmp::width) {
                top_or_right = beg_x_interval*counter++;
                ranges.emplace_back(bottom_or_left, top_or_right);
                bottom_or_left = top_or_right;
            }
            if (ranges.size() > num_beg_x_intervals)
                ranges.pop_back();
            ranges.back().second = bmp::width;
            std::vector<std::pair<long double, long double>> end_ranges(0); // unlikely not to be used
            if (end_beg_row < bmp::height) {
                long double end_x_interval = bmp::width/num_end_x_intervals;
                bottom_or_left = 0;
                top_or_right = 0;
                counter = 1;
                while (top_or_right <= bmp::width) {
                    top_or_right = end_x_interval*counter++;
                    end_ranges.emplace_back(bottom_or_left, top_or_right);
                    bottom_or_left = top_or_right;
                }
                if (end_ranges.size() > num_end_x_intervals)
                    end_ranges.pop_back();
                end_ranges.back().second = bmp::width;
            }
            try {
                unsigned int y;
                unsigned int next_y;
                for (const auto &[x_start, x_end] : ranges) {
                    for (y = 0, next_y = 1; y < end_beg_row; ++y, ++next_y)
                        threads.emplace_back(&astro_scene::render_subp, this, x_start, x_end, y, next_y);
                }
                for (const auto &[x_start, x_end] : end_ranges) { // simply won't execute if end_ranges is empty
                    for (y = end_beg_row, next_y = end_beg_row + 1; y < bmp::height; ++y, ++next_y)
                        threads.emplace_back(&astro_scene::render_subp, this, x_start, x_end, y, next_y);
                }
                for (std::thread &t : threads)
                    t.join();
                max_calculator.join();
            } catch (const std::system_error &e) {
                tot_threads = bmp::height;
                std::latch again{tot_threads};
                thread_latch = &again;
                threads.clear();
                try { // if the above raises an exception, here try to revert to row-by-row only
                    for (unsigned int y = 0; y < bmp::height; ++y)
                        threads.emplace_back(&astro_scene::render_subp, this, 0, bmp::width, y, y + 1);
                    for (std::thread &t : threads)
                        t.join();
                    max_calculator.join();
                } catch (const std::system_error &e2) {
                    tot_threads = 1;
                    std::latch once{1};
                    thread_latch = &once;
                    this->render_subp(0, bmp::width, 0, bmp::height);
                    max_calculator.join();
                    return false;
                }
            }
            return true; // true if multithreading could be performed
        }
        void clear_bodies() noexcept {
            bodies.clear();
            body_clrs.clear();
        }
        void clear_stars() noexcept {
            stars.clear();
            star_clrs.clear();
        }
        void clear_both() noexcept {
            this->clear_bodies();
            this->clear_stars();
        }
        void clear() override {
            bmp::clear();
            this->dealloc();
            cleared = true;
        }
        void reallocate() override {
            bmp::reallocate();
            this->alloc();
        }
        ~astro_scene() override {
            this->dealloc(); // bmp::~bmp() is automatically called after this
        }
    };
    bool operator==(const image_dimensions &dims1, const image_dimensions &dims2) {
        return dims1.x == dims2.x && dims1.y == dims2.y;
    }
    bool operator!=(const image_dimensions &dims1, const image_dimensions &dims2) {
        return dims1.x != dims2.x || dims1.y != dims2.y;
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU>
    std::ostream &operator<<(std::ostream &os, const ray<PosU, DirU, LenU> &r) {
        return os << '(' << r.pos << ") + " << r.l  << '(' << static_cast<const vector3D<DirU>&>(r) << ')';
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU>
    std::ostream &operator<<(std::ostream &os, const ray<PosU, DirU, LenU> &&r) {
        return os << '(' << r.pos << ") + " << r.l  << '(' << static_cast<const vector3D<DirU>&>(r) << ')';
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU, isNumWrapper LumU>
    std::ostream &operator<<(std::ostream &os, const light_src<PosU, DirU, LenU, LumU> &src) {
        return os << "[gtd::light_src@" << &src << ":pos=(" << src.pos << "),lum=" << src.lum << ']';
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator==(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &r2) {
        return r1.pos == r2.pos && r1.l == r2.l &&
               static_cast<const vector3D<DirU>&>(r1) == static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator==(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &&r2) {
        return r1.pos == r2.pos && r1.l == r2.l &&
               static_cast<const vector3D<DirU>&>(r1) == static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator==(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &r2) {
        return r1.pos == r2.pos && r1.l == r2.l &&
               static_cast<const vector3D<DirU>&>(r1) == static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator==(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &&r2) {
        return r1.pos == r2.pos && r1.l == r2.l &&
               static_cast<const vector3D<DirU>&>(r1) == static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator!=(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &r2) {
        return r1.pos != r2.pos || r1.l != r2.l ||
               static_cast<const vector3D<DirU>&>(r1) != static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator!=(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &&r2) {
        return r1.pos != r2.pos || r1.l != r2.l ||
               static_cast<const vector3D<DirU>&>(r1) != static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator!=(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &r2) {
        return r1.pos != r2.pos || r1.l != r2.l ||
               static_cast<const vector3D<DirU>&>(r1) != static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator!=(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &&r2) {
        return r1.pos != r2.pos || r1.l != r2.l ||
               static_cast<const vector3D<DirU>&>(r1) != static_cast<const vector3D<DirV>&>(r2);
    }
    namespace col_gens {
        template <unsigned char min = 0, unsigned char max = 255>
        color grayscale() {
            static_assert(min <= max, "Min. color value passed must be equal to or less than max. color value.\n");
            static unsigned char val;
            static std::uniform_int_distribution<unsigned short> dist{(unsigned short) min,(unsigned short) max};//[a,b]
            static time_t seed_val = 1;
            static time_t now;
            now = time(nullptr);
            seed_val = now == seed_val ? ++seed_val : now;
            srand(seed_val);
            std::mt19937 mt{(unsigned int) rand()};
            val = dist(mt);
            return {val, val, val};
        }
        template <color col>
        color monoscale() {
            static std::uniform_real_distribution<long double> dist{0, 1}; // [a, b)
            static time_t seed_val = 1;
            static time_t now;
            now = time(nullptr);
            seed_val = now == seed_val ? ++seed_val : now;
            srand(seed_val);
            std::mt19937 mt{(unsigned int) rand()};
            return dist(mt)*col;
        }
    }
    typedef ray<long double, long double, long double> ray_t;
    typedef light_src<long double, long double, long double, long double> src_t;
    typedef star<long double, long double, long double, long double, long double, long double, long double, 0> star_t;
    typedef camera<long double, long double, long double> cam;
    typedef astro_scene<long double, long double, long double, long double, long double, long double, long double,
                        long double, true, 0> asc_0f;
    typedef astro_scene<long double, long double, long double, long double, long double, long double, long double,
            long double, true, 1> asc_1f;
    typedef astro_scene<long double, long double, long double, long double, long double, long double, long double,
            long double, true, 10> asc_10f;
    typedef astro_scene<long double, long double, long double, long double, long double, long double, long double,
            long double, true, 100> asc_100f;
    typedef astro_scene<long double, long double, long double, long double, long double, long double, long double,
            long double, true, 1000> asc_1000f;
    typedef astro_scene<long double, long double, long double, long double, long double, long double, long double,
                        long double, false, 0> asc_n;
}
#endif
