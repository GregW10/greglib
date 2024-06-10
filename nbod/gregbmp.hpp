//
// Created by mario on 10/10/2022.
//

#ifndef GREGBMP_H
#define GREGBMP_H

#ifndef __cplusplus
#error "The gregbmp.hpp header file is a C++ header file only."
#endif

#include <iostream>
#include <fstream>
#include <cmath>

#include "gregstr.hpp"
#include "gregmat.hpp"

namespace gtd {
#pragma pack(push, 1)
    struct color {
        unsigned char b;
        unsigned char g;
        unsigned char r;
    };
#pragma pack(pop)
    inline std::ostream &operator<<(std::ostream &os, const color &col);
    template <isNumWrapper T>
    inline color operator*(const T &scalar, const color &col);
    template <isNumWrapper T>
    inline color operator*(const color &col, const T &scalar);
    /* template <isNumWrapper T>
    inline color operator*(const T &&scalar, const color &col);
    template <isNumWrapper T>
    inline color operator*(const color &&col, const T &scalar);
    template <isNumWrapper T>
    inline color operator*(const T &scalar, const color &&col);
    template <isNumWrapper T>
    inline color operator*(const color &col, const T &&scalar);
    template <isNumWrapper T>
    inline color operator*(const T &&scalar, const color &&col);
    template <isNumWrapper T>
    inline color operator*(const color &&col, const T &&scalar); */
    namespace colors {
        color black{0, 0, 0};
        color white{255, 255, 255};
        color blue{255, 0, 0};
        color green{0, 255, 0};
        color red{0, 0, 255};
        color pink{180, 105, 255};
        color cerise{99, 49, 222};
        color fuchsia{255, 0, 255};
        color neon_pink{240, 16, 255};
        color pink_orange{128, 152, 248};
        color purple{128, 0, 128};
        color salmon{114, 128, 250};
        color watermelon_pink{131, 115, 227};
        color orange{0, 165, 255};
        color gold{0, 215, 255};
        color yellow{0, 255, 255};
        color lavender{250, 230, 230};
        color indigo{130, 0, 75};
        color violet{238, 130, 238};
        color lime_green{50, 205, 50};
        color forest_green{34, 139, 34};
        color dark_green{0, 100, 0};
        color aqua{255, 255, 0};
        color sky_blue{235, 206, 135};
        color royal_blue{225, 105, 65};
        color navy{128, 0, 0};
        color wheat{179, 222, 245};
        color tan{140, 180, 210};
        color rosy_brown{143, 143, 188};
        color peru{63, 133, 205};
        color chocolate{30, 105, 210};
        color brown{42, 42, 165};
        color maroon{0, 0, 128};
        color snow{250, 250, 255};
        color honey_dew{240, 255, 240};
        color azure{255, 255, 240};
        color ghost_white{255, 248, 248};
        color beige{220, 245, 245};
        color ivory{240, 255, 255};
        color gainsboro{220, 220, 220};
        color silver{192, 192, 192};
        color gray{128, 128, 128};
        color slate_gray{144, 128, 112};
    }
    class zero_size_bmp : public std::invalid_argument {
    public:
        zero_size_bmp() : std::invalid_argument{"The resulting BMP would have a size of zero. Width and height must be "
                                                "non-zero.\n"} {}
    };
    class invalid_bmp_format : public std::exception {
        String msg;
    public:
        invalid_bmp_format() : msg{"BMP is not in the correct format.\n"} {}
        explicit invalid_bmp_format(const char *str) : msg{str} {}
        const char *what() const noexcept override {
            return msg.c_str();
        }
    };
    class point {
    public:
        long double x{0}; // no good reason to make these private
        long double y{0};
        point() = default;
        point(const long double &x_pos, const long double &y_pos) : x{x_pos}, y{y_pos} {}
        point(const long double &&x_pos, const long double &&y_pos) : x{x_pos}, y{y_pos} {}
        point(const point &other) : x{other.x}, y{other.y} {}
        static long double distance_between(const point &p1, const point &p2) {
            return sqrtl((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
        }
        friend inline bool operator==(const point &p1, const point &p2);
    };
    class rectangle;
    class shape {
    protected:
        point pos{};
    public:
        shape() = default;
        shape(const long double &x_pos, const long double &y_pos) : pos{x_pos, y_pos} {}
        shape(const long double &&x_pos, const long double &&y_pos) : pos{x_pos, y_pos} {}
        explicit shape(const point &position) : pos{position} {}
        explicit shape(const point &&position) : pos{position} {}
        shape(const shape &other) = default;
        virtual bool contains(const point &p) const noexcept = 0;
        virtual bool contains(const long double &x, const long double &y) const noexcept = 0;
        virtual bool contains(const long double &&x, const long double &&y) const noexcept = 0;
        virtual rectangle get_rect_bounds() const noexcept = 0;
        const point &get_pos() const noexcept {
            return pos;
        }
        virtual void set_pos(const point &new_pos) {
            pos = new_pos;
        }
        virtual void set_pos(const point &&new_pos) {
            pos = new_pos;
        }
        virtual void set_pos(const long double &x, const long double &y) {
            pos.x = x;
            pos.y = y;
        }
        virtual void set_pos(const long double &&x, const long double &&y) {
            set_pos(x, y);
        }
        friend class bmp;
    };
    class rectangle : public shape {
        long double width{};
        long double height{};
        void check_dim() {
            if (width < 0) {
                width = 0;
            }
            if (height < 0) {
                height = 0;
            }
        }
    public:
        rectangle() = default;
        rectangle(const long double &x_pos, const long double &y_pos) : shape{x_pos, y_pos} {}
        rectangle(const long double &&x_pos, const long double &&y_pos) : shape{x_pos, y_pos} {}
        rectangle(const long double &x_pos, const long double &y_pos, const long double &w, const long double &h) :
        shape{x_pos, y_pos}, width{w}, height{h} {check_dim();}
        rectangle(const long double &&x_pos, const long double &&y_pos, const long double &&w, const long double &&h) :
        shape{x_pos, y_pos}, width{w}, height{w} {check_dim();}
        rectangle(const point &position, const long double &w, const long double &h) :
        shape{position}, width{w}, height{h} {check_dim();}
        rectangle(const point &&position, const long double &w, const long double &h) :
        shape{position}, width{w}, height{h} {check_dim();}
        rectangle(const rectangle &other) = default;
        virtual void set_width(const long double &w) noexcept {
            if (w < 0)
                return;
            width = w;
        }
        virtual void set_width(const long double &&w) noexcept {
            set_width(w);
        }
        virtual void set_height(const long double &h) noexcept {
            if (h < 0)
                return;
            height = h;
        }
        virtual void set_height(const long double &&h) noexcept {
            set_height(h);
        }
        long double get_width() const noexcept {
            return width;
        }
        long double get_height() const noexcept {
            return height;
        }
        bool contains(const long double &x, const long double &y) const noexcept override {
            if (x >= shape::pos.x && x < shape::pos.x + width && y >= shape::pos.y && y < shape::pos.y + height)
                return true;
            return false;
        }
        bool contains(const long double &&x, const long double &&y) const noexcept override {
            return contains(x, y);
        }
        bool contains(const point &p) const noexcept override {
            if (p.x >= shape::pos.x && p.x < shape::pos.x + width && p.y >= shape::pos.y && p.y < shape::pos.y + height)
                return true;
            return false;
        }
        rectangle get_rect_bounds() const noexcept override {
            return {*this};
        }
        friend class bmp;
    };
    class square : public rectangle {
    public:
        square() = default;
        square(const long double &x_pos, const long double &y_pos) : rectangle{x_pos, y_pos} {}
        square(const long double &&x_pos, const long double &&y_pos) : rectangle{x_pos, y_pos} {}
        square(const long double &x_pos, const long double &y_pos, const long double &side_length) :
        rectangle{x_pos, y_pos, side_length, side_length} {}
        square(const long double &&x_pos, const long double &&y_pos, const long double &&side_length) :
        rectangle{x_pos, y_pos, side_length, side_length} {}
        square(const point &position, const long double &side_length) :
        rectangle{position, side_length, side_length} {}
        square(const point &&position, const long double &&side_length) :
        rectangle{position, side_length, side_length} {}
        void set_width(const long double &new_side_length) noexcept override {
            rectangle::set_width(new_side_length);
            rectangle::set_height(new_side_length);
        }
        void set_height(const long double &new_side_length) noexcept override {
            rectangle::set_width(new_side_length);
            rectangle::set_height(new_side_length);
        }
        void set_width(const long double &&new_side_length) noexcept override {
            set_width(new_side_length);
        }
        void set_height(const long double &&new_side_length) noexcept override {
            set_height(new_side_length);
        }
        friend class bmp;
    };
    class pixel;
    class circle : public shape {
        long double radius{};
    public:
        circle() = default;
        circle(const long double &x_pos, const long double &y_pos) : shape(x_pos, y_pos) {}
        circle(const long double &&x_pos, const long double &&y_pos) : shape(x_pos, y_pos) {}
        explicit circle(const point &p) : shape(p) {}
        explicit circle(const point &&p) : shape(p) {}
        circle(const long double &x_pos, const long double &y_pos, const long double &r) :
        shape(x_pos, y_pos), radius{r} {}
        circle(const long double &&x_pos, const long double &&y_pos, const long double &&r) :
        shape(x_pos, y_pos), radius{r} {}
        explicit circle(const point &p, const long double &r) :
        shape(p), radius{r} {}
        explicit circle(const point &&p, const long double &&r) :
        shape(p), radius{r} {}
        circle(const circle &other) = default;
        bool set_radius(const long double &r) noexcept {
            if (r < 0)
                return false;
            radius = r;
            return true;
        }
        bool set_radius(const long double &&r) noexcept {
            return set_radius(r);
        }
        bool contains(const long double &x, const long double &y) const noexcept override {
            return sqrtl((x - shape::pos.x)*(x - shape::pos.x) + (y - shape::pos.y)*(y - shape::pos.y)) <= radius;
        }
        bool contains(const long double &&x, const long double &&y) const noexcept override {
            return contains(x, y);
        }
        bool contains(const point &p) const noexcept override {
            return sqrtl((p.x - shape::pos.x)*(p.x - shape::pos.x)+(p.y - shape::pos.y)*(p.y - shape::pos.y)) <= radius;
        }
        bool fully_contains(const pixel &p) const noexcept;
        bool partially_contains(const pixel &p) const noexcept;
        rectangle get_rect_bounds() const noexcept override {
            return {shape::pos.x - radius, shape::pos.y - radius, 2*radius, 2*radius};
        }
        friend class bmp;
    };
    class pixel : public square { // essentially a square with side length of 1
    private:
        color col{};
        void correct_pos() {
            if (shape::pos.x < 0) { // pixel can't have a negative position
                shape::pos.x = -shape::pos.x;
            }
            if (shape::pos.y < 0) {
                shape::pos.y = -shape::pos.y;
            }
            shape::pos.x = floor(shape::pos.x);
            shape::pos.y = floor(shape::pos.y);
        }
    public:
        pixel() = default;
        pixel(const long double &x_pos, const long double &y_pos) : square{x_pos, y_pos, 1.0l} {correct_pos();}
        pixel(const long double &&x_pos, const long double &&y_pos) : square{x_pos, y_pos, 1.0l} {correct_pos();}
        pixel(const long double &x_pos, const long double &y_pos, const color &color) :
        square{x_pos, y_pos, 1.0l}, col{color} {correct_pos();}
        pixel(const long double &&x_pos, const long double &&y_pos, const color &&color) :
        square{x_pos, y_pos, 1.0l}, col{color} {correct_pos();}
        pixel(const pixel &other) = default;
        void set_color(const color &clr) {
            col = clr;
        }
        const color &get_color() const noexcept {
            return col;
        }
        void set_pos(const point &p) noexcept override {
            shape::pos = p;
            correct_pos();
        }
        void set_pos(const point &&p) noexcept override {
            set_pos(p);
        }
        void set_pos(const long double &x, const long double &y) noexcept override {
            shape::pos.x = x;
            shape::pos.y = y;
            correct_pos();
        }
        void set_pos(const long double &&x, const long double &&y) noexcept override {
            set_pos(x, y);
        }
        pixel &calc_shade(const circle &c, const color &bg_color = color{}) noexcept {
            if (!c.partially_contains(*this)) {
                col = bg_color; // the whole pixel is outside the circle, so its colour must be the background colour
                return *this;
            }
            if (c.fully_contains(*this)) {
                return *this; // the whole pixel is inside the circle, so its colour remains unchanged
            }
            point p;
            unsigned short count = 0;
            for (unsigned char i = 1; i <= 16; ++i) { // creates a 16x16 grid of points on the pixel - the proportion
                for (unsigned char j = 1; j <= 16; ++j) { // that the circle contains sets the brightness of the colour
                    p.x = shape::pos.x + i/17.0l;
                    p.y = shape::pos.y + j/17.0l;
                    if (c.contains(p))
                        ++count;
                }
            }
            if (count == 256) // case for the pixel being entirely within circle, so fully coloured
                return *this;
            col.r = bg_color.r + (unsigned char)
                    roundl(col.r >= bg_color.r ? (col.r - bg_color.r)*(count/256.0l) : (bg_color.r - col.r)*((256.0l - count)/256.0l));
            col.g = bg_color.g + (unsigned char)
                    roundl(col.g >= bg_color.g ? (col.g - bg_color.g)*(count/256.0l) : (bg_color.g - col.g)*((256.0l - count)/256.0l));
            col.b = bg_color.b + (unsigned char)
                    roundl(col.b >= bg_color.b ? (col.b - bg_color.b)*(count/256.0l) : (bg_color.b - col.b)*((256.0l - count)/256.0l));
            return *this;
        }
        friend inline bool operator==(const pixel &p1, const pixel &p2);
        friend class circle;
        friend class bmp;
    };
    bool circle::fully_contains(const pixel &p) const noexcept {
        return contains(p.pos) && contains(p.pos.x, p.pos.y + 1) && contains(p.pos.x + 1, p.pos.y) &&
               contains(p.pos.x + 1, p.pos.y + 1);
    }
    bool circle::partially_contains(const pixel &p) const noexcept {
        return (contains(p.pos) || contains(p.pos.x, p.pos.y + 1) || contains(p.pos.x + 1, p.pos.y) ||
               contains(p.pos.x + 1, p.pos.y + 1)) ||
               (/* radius <= sqrtl(0.5l) && */p.pos.x <= shape::pos.x && p.pos.x + 1 >= shape::pos.x &&
               p.pos.y <= shape::pos.y && p.pos.y + 1 >= shape::pos.y);
    }
    class bmp {
    protected:
        color **data = nullptr;
        String path;
        unsigned int width = def_width;
        unsigned int height = def_height;
        color sc_clr = colors::black; // selected colour
    private:
        static constexpr unsigned int def_width = 2400;
        static constexpr unsigned int def_height = 1600;
        unsigned int stroke_t = 10;
        unsigned char padding[3] = {0, 0, 0};
        unsigned char (*calc_pad)(const unsigned int&) =
                [](const unsigned int &w){return (unsigned char) ((w*3) % 4 == 0 ? 0 : 4 - ((w*3) % 4));};
        void fill_background() {
            data = new color*[height];
            color **ptr = data;
            color *col;
            for (unsigned int j = 0; j < height; ++j, ++ptr) {
                *ptr = new color[width];
                col = *ptr;
                for (unsigned int i = 0; i < width; ++i, ++col) {
                    *col = sc_clr;
                }
            }
        }
        void check_size() const {
            if (width == 0 || height == 0) {
                throw zero_size_bmp();
            }
        }
        void dealloc() {
            if (data == nullptr)
                return;
            color **ptr = data;
            for (unsigned int j = 0; j < height; ++j, ++ptr)
                delete [] *ptr;
            delete [] data;
            data = nullptr;
        }
        void alloc() {
            dealloc();
            data = new color*[height];
            color **ptr = data;
            for (unsigned int j = 0; j < height; ++j, ++ptr)
                *ptr = new color[width];
        }
        unsigned char *get_hdr() const noexcept {
            unsigned int size = file_size();
            static unsigned char header[14] = {'B', 'M', (unsigned char) (size), (unsigned char) (size >> 8),
                                               (unsigned char) (size >> 16), (unsigned char) (size >> 24), 0, 0, 0, 0,
                                               54, 0, 0, 0};
            return header;
        }
        unsigned char *get_info_hdr() const noexcept {
            static unsigned char info_hdr[40] = {40, 0, 0, 0, (unsigned char) width, (unsigned char) (width >> 8),
                                                 (unsigned char) (width >> 16), (unsigned char) (width >> 24),
                                                 (unsigned char) height, (unsigned char) (height >> 8), (unsigned char)
                                                 (height >> 16), (unsigned char) (height >> 24), 1, 0, 24, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            return info_hdr;
        }
        void read_bmp(const char *source) {
            if (source == nullptr || !*source)
                throw std::invalid_argument("A nullptr or empty string cannot be passed as a source BMP path.\n");
            std::ifstream in(source, std::ios_base::binary | std::ios_base::in);
            if (!in.good()) {
                String msg = "There was an error opening the file \"";
                msg.append_back(source).append_back("\".\n");
                throw std::invalid_argument(msg.c_str());
            }
            unsigned int offset;
            unsigned int info_hdr_size;
            unsigned short bpp;
            in.seekg(10, std::ios_base::beg);
            in.read((char *) &offset, 4);
            in.read((char *) &info_hdr_size, 4);
            in.read((char *) &width, 4);
            in.read((char *) &height, 4);
            in.seekg(2, std::ios_base::cur);
            in.read((char *) &bpp, 2);
            if (bpp != 24)
                throw invalid_bmp_format("The specified BMP must have a pixel depth of 24 bpp.\n");
            unsigned char pad = calc_pad(width);
            if (!in.good()) {
                bad:
                String msg = "There was an error reading from the file \"";
                msg.append_back(source).append_back("\".\n");
                throw std::invalid_argument(msg.c_str());
            }
            if (offset != 54 || info_hdr_size != 40)
                throw invalid_bmp_format();
            check_size();
            in.seekg(54, std::ios_base::beg);
            alloc();
            color **ptr = data;
            std::streamsize triple_width = 3*width;
            for (unsigned int j = 0; j < height; ++j) {
                in.read((char *) *ptr++, triple_width);
                in.seekg(pad, std::ios_base::cur);
            }
            if (!in.good())
                goto bad;
            in.close();
        }
        void copy_pixels(const color *const *source) {
            if (source == nullptr)
                return;
            for (unsigned int j = 0; j < height; ++j)
                for (unsigned int i = 0; i < width; ++i)
                    data[j][i] = source[j][i];
        }
    protected:
        explicit bmp(String &&str) : path{std::move(str)} {fill_background();}
        explicit bmp(const char *source_bmp, String &&output_path) : path{output_path} {
            read_bmp(source_bmp);
        }
    public:
        bmp() : path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_" + get_date_and_time() + ".bmp")}
        {fill_background();}
        explicit bmp(const char *source_bmp) :
        path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_" + get_date_and_time() + ".bmp")} {
            read_bmp(source_bmp);
        }
        bmp(unsigned int bmp_width, unsigned int bmp_height) : width{bmp_width}, height{bmp_height},
        path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_" + get_date_and_time() + ".bmp")} {
            check_size();
            fill_background();
        }
        explicit bmp(const color &background_color) : sc_clr{background_color},
        path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_" + get_date_and_time() + ".bmp")}
        {fill_background();}
        bmp(color background_color, unsigned int bmp_width, unsigned int bmp_height) : sc_clr{background_color},
        width{bmp_width}, height{bmp_height},
        path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_" + get_date_and_time() + ".bmp")} {
            check_size();
            fill_background();
        }
        bmp(const String &bmp_path, color background_color, unsigned int bmp_width, unsigned int bmp_height) :
        path{bmp_path}, sc_clr{background_color}, width{bmp_width},
            height{bmp_height} {
            check_size();
            fill_background();
        }
        bmp(String &&bmp_path, color background_color, unsigned int bmp_width, unsigned int bmp_height) :
        path{std::move(bmp_path)}, sc_clr{background_color}, width{bmp_width}, height{bmp_height} {
            check_size();
            fill_background();
        }
        bmp(const bmp &other) : width{other.width}, height{other.height}, path{other.path}, sc_clr{other.sc_clr},
        stroke_t{other.stroke_t} {
            copy_pixels(other.data);
        }
        bmp(bmp &&other) noexcept : width{other.width}, height{other.height}, path{std::move(other.path)},
        sc_clr{other.sc_clr}, stroke_t{other.stroke_t} {
            data = other.data;
            other.data = nullptr;
        }
        static bmp create_bmp(const char *source_bmp) {
            return bmp(source_bmp); // can't use copy-list-initialisation as I have marked this ctor explicit
        }
        unsigned int get_height() const noexcept {
            return height;
        }
        unsigned int get_width() const noexcept {
            return width;
        }
        void set_color(const color &col) noexcept {
            sc_clr = col;
        }
        void set_color(const color &&col) noexcept {
            sc_clr = col;
        }
        void set_color(unsigned char r, unsigned char g, unsigned char b) {
            sc_clr.r = r;
            sc_clr.g = g;
            sc_clr.b = b;
        }
        bool set_pixel(unsigned int x, unsigned int y, const color &col) {
            if (x >= width || y >= height)
                return false;
            data[y][x] = col;
            return true;
        }
        bool set_pixel(unsigned int x, unsigned int y, const color &&col) {
            return set_pixel(x, y, col);
        }
        bool set_pixel(unsigned int x, unsigned int y) {
            return set_pixel(x, y, sc_clr);
        }
        color get_color() {
            return sc_clr;
        }
        color get_color(unsigned int x, unsigned int y) {
            if (x >= width || y >= height)
                throw std::out_of_range("The requested pixel is out of range.\n");
            return data[y][x];
        }
        bool fill_rect(unsigned int x, unsigned int y, unsigned int w, unsigned int h) {
            if (x >= width || y >= height || !w || !h)
                return false;
            if (x + w - 1 >= width)
                w = width - x;
            if (y + h - 1 >= height)
                h = height - y;
            unsigned int x_end = x + w;
            unsigned int y_end = y + h;
            for (unsigned int j = y; j < y_end; ++j)
                for (unsigned int i = x; i < x_end; ++i)
                    data[j][i] = sc_clr;
            return true;
        }
        void fill_bg() { // convenience method
            fill_rect(0, 0, width, height);
        }
        unsigned int get_rgb(unsigned int x, unsigned int y) {
            if (x >= width || y >= height)
                return -1;
            return (data[y][x].r << 16) + (data[y][x].g << 8) + data[y][x].b;
        }
        unsigned char *as_ycbcr(unsigned int x, unsigned int y) const noexcept {
            if (x >= width || y >= height)
                return nullptr;
            /* thanks to: https://web.archive.org/web/20180423091842/http://www.equasys.de/colorconversion.html for
             * the conversion equations */
            static unsigned char ret[3];
            unsigned char &r = data[y][x].r;
            unsigned char &g = data[y][x].g;
            unsigned char &b = data[y][x].b;
            ret[0] = (unsigned char) roundl(0.299l*r + 0.587l*g + 0.114l*b);
            ret[1] = (unsigned char) roundl(128.0l - 0.169l*r - 0.331l*g + 0.5l*b);
            ret[2] = (unsigned char) roundl(128.0l + 0.5l*r - 0.419l*g - 0.081l*b);
            return ret;
        }
        bool contains(const pixel &p) const noexcept {
            return p.pos.x >= 0 && p.pos.y >= 0 && p.pos.x < width && p.pos.y < height;
        }
        void draw_circle(const circle &c, bool fill = true) {
            if (c.radius == 0 || c.pos.x + c.radius < 0 || c.pos.y + c.radius < 0 || c.pos.x - c.radius >= width ||
                c.pos.y - c.radius >= height) { // case for where the entire circle would be hidden
                return;
            }
            pixel p;
            if (fill) {
                long double top = c.pos.y + c.radius;
                if (top > height - 1)
                    top = height;
                unsigned int left;
                unsigned int right;
                long double root;
                long double y_btm = c.pos.y - c.radius;
                if (y_btm < 0)
                    y_btm = 0;
                long double y_counter = y_btm;
                long double center_bound;
                long double left_ld;
                long double right_ld;
                unsigned int count = 0;
                bool skip_right;
                while (y_counter < top) { // first loop fills in all the totally covered pixels
                    skip_right = false;
                    root = sqrtl(c.radius * c.radius - (c.pos.y - y_counter) * (c.pos.y - y_counter));
                    left_ld = roundl(c.pos.x - root);
                    right_ld = roundl(c.pos.x + root);
                    if (left_ld >= width || right_ld - left_ld <= 0 || right_ld < 0)
                        goto end;
                    left = (unsigned int) (left_ld < 0 ? 0 : left_ld);
                    right = (unsigned int) (right_ld >= width ? width - 1 : right_ld);
                    center_bound = c.pos.x;
                    if (c.pos.x >= width) {
                        center_bound = width - 1;
                        skip_right = true;
                    }
                    while (left <= center_bound) { // cast to uns. int takes care of correct conversion to axes coords
                        p.set_pos(left, (unsigned int) y_counter);
                        if (c.fully_contains(p) && this->contains(p))
                            data[(unsigned int) y_counter][left] = sc_clr;
                        ++left;
                    }
                    if (skip_right)
                        goto end;
                    if (right >= width)
                        right = width - 1;
                    if (center_bound < 0)
                        center_bound = 0;
                    while (right > center_bound) {
                        p.set_pos(right, (unsigned int) y_counter);
                        if (c.fully_contains(p) && this->contains(p))
                            data[(unsigned int) y_counter][right] = sc_clr;
                        --right;
                    }
                    end:
                    ++count;
                    y_counter = y_btm + count; // to avoid accumulated f.p. error
                }
            }
            pixel prev;
            prev.pos.x = 0.5; // this way prev is guaranteed to be different from p in the first iteration
            prev.pos.y = 0.5;
            std::vector<pixel> pixels(8); // pixels surrounding the one being adjusted
            long double theta = 0;
            long double d_theta = 0.000244140625/c.radius; // a small enough angle step to include all outer pixels
            long double two_pi = 2*PI; // PI is defined in gregmat.h
            color background_color{};
            unsigned short tot_r;
            unsigned short tot_g;
            unsigned short tot_b;
            unsigned char bg_count;
            while (theta < two_pi) {//this loop fills in missing pixels and/or corrects existing ones with anti-aliasing
                p.set_pos((unsigned int) (c.pos.x + cosl(theta)*c.radius),
                          (unsigned int) (c.pos.y + sinl(theta)*c.radius));
                if (p == prev)
                    goto bye;
                p.set_color(sc_clr);
                tot_r = tot_g = tot_b = 0;
                pixels[0].set_pos(p.pos.x, p.pos.y + 1);
                pixels[1].set_pos(p.pos.x - 1, p.pos.y + 1);
                pixels[2].set_pos(p.pos.x - 1, p.pos.y);
                pixels[3].set_pos(p.pos.x - 1, p.pos.y - 1);
                pixels[4].set_pos(p.pos.x, p.pos.y - 1);
                pixels[5].set_pos(p.pos.x + 1, p.pos.y - 1);
                pixels[6].set_pos(p.pos.x + 1, p.pos.y);
                pixels[7].set_pos(p.pos.x + 1, p.pos.y + 1);
                bg_count = 0;
                for (const pixel &pix : pixels) {
                    if (!c.partially_contains(pix) && contains(pix)) {
                        tot_r += data[(unsigned int) pix.pos.y][(unsigned int) pix.pos.x].r;
                        tot_g += data[(unsigned int) pix.pos.y][(unsigned int) pix.pos.x].g;
                        tot_b += data[(unsigned int) pix.pos.y][(unsigned int) pix.pos.x].b;
                        ++bg_count;
                    }
                }
                if (!bg_count)
                    background_color = sc_clr;
                else {
                    background_color.r = (unsigned char) (tot_r/((long double) bg_count) + 0.5l);
                    background_color.g = (unsigned char) (tot_g/((long double) bg_count) + 0.5l);
                    background_color.b = (unsigned char) (tot_b/((long double) bg_count) + 0.5l);
                }
                p.calc_shade(c, background_color);
                if (contains(p))
                    data[(unsigned int) roundl(p.pos.y - 0.5l)][(unsigned int) roundl(p.pos.x - 0.5l)] = p.get_color();
                bye:
                theta += d_theta; // f.p. rounding errors are not an issue here
                prev = p;
            }
        }
        void make_mono(unsigned char cutoff, color dark_color = colors::black, color bright_color = colors::white) {
            for (unsigned int j = 0; j < height; ++j) {
                for (unsigned int i = 0; i < width; ++i) {
                    if ((data[j][i].r + data[j][i].g + data[j][i].b)/3.0l > cutoff) {
                        data[j][i] = bright_color;
                    }
                    else {
                        data[j][i] = dark_color;
                    }
                }
            }
        }
        void make_grayscale() {
            unsigned char avg;
            for (unsigned int j = 0; j < height; ++j) {
                for (unsigned int i = 0; i < width; ++i) {
                    avg = (unsigned char) ((data[j][i].r + data[j][i].g + data[j][i].b)/3.0l);
                    data[j][i].r = data[j][i].g = data[j][i].b = avg;
                }
            }
        }
        virtual void clear() {
            dealloc();
        }
        bool set_path(const char *str) {
            if (str == nullptr || *str == 0) {
                return false;
            }
            path = str;
            return true;
        }
        const char *get_path() const noexcept {
            return path.c_str();
        }
        unsigned int file_size() const noexcept { // has to be unsigned int, otherwise the BMP format wouldn't allow it
            return 54 + 3*width*height + height*calc_pad(width);
        }
        virtual void reallocate() {
            if (data != nullptr) {
                return;
            }
            alloc();
        }
        matrix<unsigned char> as_matrix(const unsigned int &x_index, const unsigned int &y_index) const {
            if (x_index >= width || y_index >= height)
                throw std::out_of_range("The requested pixel is out of range.\n");
            const color &ref = data[y_index][x_index];
            return matrix<unsigned char>(3, 1) << ref.r << ref.g << ref.b;
        }
        matrix<unsigned char> as_matrix(const unsigned int &&x_index, const unsigned int &&y_index) const {
            return as_matrix(x_index, y_index);
        }
        static matrix<unsigned char> as_matrix(const color &col) {
            return matrix<unsigned char>(3, 1) << col.r << col.g << col.b;
        }
        bool read(const char *source_bmp) {
            if (source_bmp == nullptr || *source_bmp == 0) {
                return false;
            }
            read_bmp(source_bmp);
            return true;
        }
        bool write(const char *path_to_bmp) {
            if (path_to_bmp == nullptr || *path_to_bmp == 0) {
                return false;
            }
            std::ofstream out(path_to_bmp, std::ios_base::trunc | std::ios_base::binary);
            if (!out.good()) {
                return false;
            }
            unsigned char pad = calc_pad(width);
            out.write((char *) get_hdr(), 14);
            out.write((char *) get_info_hdr(), 40);
            color **ptr = data;
            std::streamsize triple_width = 3*width;
            for (unsigned int j = 0; j < height; ++j) {
                out.write((char *) *ptr++, triple_width);
                out.write((char *) padding, pad);
            }
            out.close();
            return true;
        }
        bool write() {
            return write(path.c_str());
        }
        int open_image() {
#ifndef _WIN32
            return std::system(("open " + path).c_str());
#else
            return std::system(path.c_str());
#endif
        }
        bmp &operator=(const bmp &other) {
            if (&other == this) {
                return *this;
            }
            bool need_to_alloc = false;
            if (other.data == nullptr || this->width != other.width || this->height != other.height) {
                this->clear();
                need_to_alloc = true;
            }
            this->width = other.width;
            this->height = other.height;
            this->sc_clr = other.sc_clr;
            this->path = other.path;
            this->stroke_t = other.stroke_t;
            if (other.data != nullptr) {
                if (need_to_alloc) {
                    alloc();
                }
                copy_pixels(other.data);
            }
            return *this;
        }
        bmp &operator=(bmp &&other) noexcept {
            if (&other == this) { // could happen!
                return *this;
            }
            this->data = other.data;
            other.data = nullptr;
            this->width = other.width;
            this->height = other.height;
            this->sc_clr = other.sc_clr;
            this->path = std::move(other.path);
            this->stroke_t = other.stroke_t;
            return *this;
        }
        color *operator[](unsigned int index) {
            if (index >= height) {
                throw std::out_of_range("The requested pixel is out of range.\n");
            }
            return data[index];
        }
        virtual ~bmp() {
            dealloc();
        }
        friend inline std::ostream &operator<<(std::ostream &os, const bmp &image);
    };
    inline std::ostream &operator<<(std::ostream &os, const bmp &image) {
        os << "[gtd::bmp@" << &image << ":width=" << image.width << ",height=" << image.height
           << ",size_allocated=";
        if (image.data == nullptr) {
            return os << "0]";
        }
        return os << image.width*image.height*3 << ']';
    }
    template <isNumWrapper T>
    inline color operator*(const T &scalar, const color &col) {
        return {(unsigned char) (scalar*col.b), (unsigned char) (scalar*col.g), (unsigned char) (scalar*col.r)};
    }
    template <isNumWrapper T>
    inline color operator*(const color &col, const T &scalar) {
        return {(unsigned char) (scalar*col.b), (unsigned char) (scalar*col.g), (unsigned char) (scalar*col.r)};
    }
    /* template <isNumWrapper T>
    inline color operator*(const T &&scalar, const color &col) {
        return {(unsigned char) (scalar*col.b), (unsigned char) (scalar*col.g), (unsigned char) (scalar*col.r)};
    }
    template <isNumWrapper T>
    inline color operator*(const color &&col, const T &scalar) {
        return {(unsigned char) (scalar*col.b), (unsigned char) (scalar*col.g), (unsigned char) (scalar*col.r)};
    }
    template <isNumWrapper T>
    inline color operator*(const T &scalar, const color &&col) {
        return {(unsigned char) (scalar*col.b), (unsigned char) (scalar*col.g), (unsigned char) (scalar*col.r)};
    }
    template <isNumWrapper T>
    inline color operator*(const color &col, const T &&scalar) {
        return {(unsigned char) (scalar*col.b), (unsigned char) (scalar*col.g), (unsigned char) (scalar*col.r)};
    }
    template <isNumWrapper T>
    inline color operator*(const T &&scalar, const color &&col) {
        return {(unsigned char) (scalar*col.b), (unsigned char) (scalar*col.g), (unsigned char) (scalar*col.r)};
    }
    template <isNumWrapper T>
    inline color operator*(const color &&col, const T &&scalar) {
        return {(unsigned char) (scalar*col.b), (unsigned char) (scalar*col.g), (unsigned char) (scalar*col.r)};
    } */
    inline std::ostream &operator<<(std::ostream &os, const color &col) {
        return os << "[gtd::color@" << &col << ":r=" << +col.r << ",g=" << +col.g << ",b=" << +col.b << ']';
    }
    inline bool operator==(const color &c1, const color &c2) {
        return c1.r == c2.r && c1.g == c2.g && c1.b == c2.b;
    }
    inline bool operator==(const point &p1, const point &p2) {
        return p1.x == p2.x && p2.y == p2.y;
    }
    inline bool operator==(const pixel &p1, const pixel &p2) {
        return p1.pos.x == p2.pos.x && p1.pos.y == p2.pos.y && p1.col == p2.col;
    }
    long double avg_bgr(color col) {
        return ((long double) (col.b + col.g + col.r))/3;
    }
}
#endif
