#ifndef GREGMAT_H
#define GREGMAT_H

#ifndef __cplusplus
#error "The gregmat.hpp header file is a C++ header file only."
#endif

#include <sstream>
#include "gregalg.hpp"

namespace gtd {
    class empty_matrix_error : public std::invalid_argument {
    public:
        empty_matrix_error() : std::invalid_argument("A matrix cannot have a size of zero.\n") {}
        explicit empty_matrix_error(const char *message) : std::invalid_argument(message) {}
    };
    class invalid_axis_error : public std::invalid_argument {
    public:
        invalid_axis_error() : std::invalid_argument("Invalid axis.\n") {}
        explicit invalid_axis_error(const char *message) : std::invalid_argument(message) {}
    };
    class invalid_matrix_format : public std::invalid_argument {
    public:
        invalid_matrix_format() : std::invalid_argument("The requested matrix has an invalid format.\n") {}
        explicit invalid_matrix_format(const char *message) : std::invalid_argument(message) {}
    };
    template <isNumWrapper T, size_t rows1, size_t rows2, size_t columns1, size_t columns2> // matrix multiplication for
    std::vector<std::vector<T>> matmul(T mat1[rows1][columns1], T mat2[rows2][columns2]) { // matrices represented as
        std::vector<std::vector<T>> ret_matrix; //                                            2D arrays
        if constexpr (columns1 != rows2 || rows1 == 0 || columns1 == 0 || rows2 == 0 || columns2 == 0)
            return ret_matrix;
        std::vector<T> sub;
        size_t total;
        for (size_t i = 0; i < rows1; ++i) {
            sub.clear();
            for (size_t j = 0; j < columns2; ++j) {
                total = 0;
                for (size_t k = 0; k < columns1; ++k)
                    total += mat1[i][k]*mat2[k][j];
                sub.push_back(total); // single element in final matrix
            }
            ret_matrix.push_back(sub); // row in final matrix
        }
        return ret_matrix;
    }
    template <isNumWrapper T, isNumWrapper U>
    auto matmul(const std::vector<std::vector<T>> &mat1, // matrix multiplication for matrices repr. as std::vectors
                const std::vector<std::vector<U>> &mat2) ->
    std::vector<std::vector<decltype(std::declval<T>()*std::declval<U>())>> {
        std::vector<std::vector<T>> ret_matrix;
        if (mat1.empty() || mat2.empty() || mat1[0].empty() || mat2[0].empty())
            return ret_matrix;
        const size_t rows1 = mat1.size();
        const size_t columns1 = mat1[0].size();
        const size_t rows2 = mat2.size();
        const size_t columns2 = mat2[0].size();
        for (const auto &[e1, e2] : zip_cref(mat1, mat2)) { // probing for valid matrix
            if (e1.size() != columns1 || e1.size() != rows2 || e2.size() != columns2) {
                return ret_matrix;
            }
        }
        std::vector<T> sub;
        size_t total;
        for (size_t i = 0; i < rows1; ++i) {
            sub.clear();
            for (size_t j = 0; j < columns2; ++j) {
                total = 0;
                for (size_t k = 0; k < columns1; ++k)
                    total += mat1[i][k]*mat2[k][j];
                sub.push_back(total);
            }
            ret_matrix.push_back(sub);
        }
        return ret_matrix;
    }
    template <isNumWrapper>
    class vector;
    template <isNumWrapper>
    class vector2D;
    template <isNumWrapper>
    class vector3D;
    template <isNumWrapper T>
    class matrix {
    private:
        std::vector<std::vector<T>> mat;
        size_t next_row = 0; // used for determining the next element at which to insert << a value
        size_t next_col = 0;
        bool auto_stop = true; // determines whether insertion stops once the last element has been set
        bool stop = false; // used to signal the insertion to stop
        inline bool is_square() const noexcept {
            return mat.size() == mat[0].size();
        }
        void check_validity() const {
            size_t column_num;
            if (mat.empty() || (column_num = mat[0].size()) == 0) {
                throw empty_matrix_error();
            }
            for (const std::vector<T> &row : mat) {
                if (row.size() != column_num) {
                    throw invalid_matrix_format("A matrix must have an equal number of elements in every row.");
                }
            }
        }
    public:
        static constexpr uint64_t nopos = -1;
        matrix() { // constructs a 2x2 identity matrix
            mat.push_back(std::vector<T>{T{1}, T{0}});
            mat.push_back(std::vector<T>{T{0}, T{1}});
            next_row = 0;
            next_col = 0;
        }
        matrix(uint64_t rows, uint64_t columns) { // constructs a zero matrix with specified size
            std::vector<T> row;
            for (uint64_t i = 0; i < rows; ++i) {
                row.clear();
                for (uint64_t j = 0; j < columns; ++j) {
                    row.push_back(T{});
                }
                mat.push_back(row);
            }
        }
        matrix(const matrix<T> &other) : mat{other.mat}, next_row{other.next_row}, next_col{other.next_col},
        auto_stop{other.auto_stop}, stop{other.stop} {}
        matrix(matrix<T> &&other) : mat{std::move(other.mat)}, next_row{other.next_row}, next_col{other.next_col},
        auto_stop{other.auto_stop}, stop{other.stop} {}
        template <typename U> requires isConvertible<U, T>
        matrix(const matrix<U> &other) : next_row{other.next_row}, next_col{other.next_col},
                                         auto_stop{other.auto_stop}, stop{other.stop} {
            std::vector<T> sub; // have to iterate through values in case U and T are not the same type
            for (const std::vector<U> &row : other) {
                sub.clear();
                for (const U &elem : row) {
                    sub.push_back(elem);
                }
                mat.push_back(sub);
            }
        }
        template <typename U> requires isConvertible<U, T>
        matrix(matrix<U> &&other) : next_row{other.next_row}, next_col{other.next_col},
                                    auto_stop{other.auto_stop}, stop{other.stop} {
            std::vector<T> sub;
            for (const std::vector<U> &row : other) {
                sub.clear();
                for (const U &elem : row) {
                    sub.push_back(elem);
                }
                mat.push_back(sub);
            }
        } // using the below two constructors is very easy thanks to C++'s initializer lists
        matrix(std::vector<std::vector<T>> &&matrix_elements) : mat(std::move(matrix_elements)) {
            check_validity();
        }
        matrix(const std::vector<std::vector<T>> &matrix_elements) : mat(matrix_elements) {
            check_validity();
        }
        matrix<T> &set_dimensions(size_t num_rows = 2, size_t num_cols = 2) {
            if (num_rows == 0 || num_cols == 0) {
                throw empty_matrix_error();
            }
            size_t rows = mat.size();
            size_t cols = mat[0].size(); // better than repeatedly calling the function (2 functions, actually)
            if (rows == num_rows && cols == num_cols) {
                return *this;
            }
            if (num_rows < rows) {
                mat.erase(mat.begin() + num_rows, mat.end());
            }
            else {
                size_t diff = num_rows - rows;
                std::vector<T> sub;
                for (size_t i = 0; i < diff; ++i) {
                    sub.clear();
                    for (size_t j = 0; j < cols; ++j) {
                        sub.push_back(T{0});
                    }
                    mat.push_back(sub);
                }
            }
            if (num_cols < cols) {
                for (std::vector<T> &row : mat) {
                    row.erase(row.begin() + num_cols, row.end());
                }
                return *this;
            }
            size_t diff = num_cols - cols;
            for (std::vector<T> &row : mat) {
                for (size_t i = 0; i < diff; ++i) {
                    row.push_back(T{0});
                }
            }
            return *this;
        }
        matrix<T> &transpose() {
            std::vector<std::vector<T>> new_mat;
            std::vector<T> new_row;
            size_t current_rows = mat.size();
            size_t current_cols = mat[0].size();
            for (size_t i = 0; i < current_cols; ++i) {
                new_row.clear();
                for (size_t j = 0; j < current_rows; ++j) {
                    new_row.push_back(mat[j][i]);
                }
                new_mat.push_back(new_row);
            }
            mat.assign(new_mat.begin(), new_mat.end());
            return *this;
        }
        matrix<T> get_transpose() const {
            return matrix<T>(*this).transpose();
        }
        matrix<T> get_sub_matrix(size_t start_row_index, size_t start_col_index,
                                 size_t num_rows = nopos, size_t num_cols = nopos) const {
            if ((start_row_index == nopos && mat.size() == nopos) || // would never happen anyway
            (start_col_index == nopos && mat[0].size() == nopos) || start_row_index + 1 > mat.size() ||
            start_col_index + 1 > mat[0].size() || num_rows == 0 || num_cols == 0) {
                throw empty_matrix_error("The resulting matrix would have a size of zero.");
            }
            if (num_rows > mat.size() - start_row_index) {
                num_rows = mat.size() - start_row_index;
            }
            if (num_cols > mat[0].size() - start_col_index) {
                num_cols = mat[0].size() - start_col_index;
            }
            matrix<T> ret_mat(num_rows, num_cols);
            size_t end_row = start_row_index + num_rows;
            size_t end_col = start_col_index + num_cols;
            for (size_t i = start_row_index, r_i = 0; i < end_row; ++i, ++r_i) {
                for (size_t j = start_col_index, r_j = 0; j < end_col; ++j, ++r_j) {
                    ret_mat.at(r_i, r_j) = mat[i][j];
                }
            }
            return ret_mat;
        }
        T determinant() const {
            if (!is_square()) {
                throw invalid_matrix_format("Only square matrices have determinants.");
            }
            if (mat.size() == 1) {
                return mat[0][0];
            }
            if (mat.size() == 2) { // base case
                return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
            }
            size_t rows = mat.size();
            size_t cols = mat[0].size();
            T ret{0};
            bool positive = true;
            std::vector<matrix<T>> matrices; // used to store sub-matrices whose determinant will be calculated
            matrix<T> sub(rows - 1, cols - 1);
            for (size_t i = 0; i < cols; ++i) { // outer loop for elements in the first row
                for (size_t j = 1, j_s = 0; j < rows; ++j, ++j_s) {
                    for (size_t k = 0, k_s = 0; k < cols; ++k, ++k_s) {
                        if (i != k) {
                            sub.at(j_s, k_s) = mat[j][k];
                            continue;
                        }
                        --k_s;
                    }
                }
                ret += positive ? mat[0][i]*sub.determinant() : -mat[0][i]*sub.determinant();
                positive = !positive;
            }
            return ret;
        }
        inline bool has_inverse() {
            return this->determinant() != 0;
        }
        template <typename U> requires isConvertible<U, T>
        matrix<T> &append_row(const std::vector<U> &row) {
            if (row.size() != mat[0].size()) {
                throw invalid_matrix_format("The row to be appended does not contain the same number of columns "
                                            "(entries)\n as the number of columns in the matrix.");
            }
            std::vector<T> row_to_add;
            for (const U &elem : row) { // since T and U might be of different types
                row_to_add.push_back(elem);
            }
            mat.push_back(row_to_add);
            return *this;
        }
        template <typename U> requires isConvertible<U, T>
        matrix<T> &append_col(const std::vector<U> &col) {
            if (col.size() != mat.size()) {
                throw invalid_matrix_format("The row to be appended does not contain the same number of columns "
                                            "(entries)\n as the number of columns in the matrix.");
            }
            std::vector<U> &col_ref = const_cast<std::vector<U>&>(col);
            for (auto &[row, col_el] : zip_ref(mat, col_ref)) {
                row.push_back(col_el);
            }
            return *this;
        }
        matrix<T> &pop_row() {
            if (mat.size() == 1) {
                return *this;
            }
            mat.pop_back();
            if (next_row == mat.size() - 1) {
                --next_row;
            }
            return *this;
        }
        matrix<T> &pop_col() {
            if (mat[0].size() == 1) {
                return *this;
            }
            for (std::vector<T> row : mat) {
                row.pop_back();
            }
            if (next_col == mat[0].size() - 1) {
                --next_col;
            }
            return *this;
        }
        size_t num_rows() {
            return mat.size();
        }
        size_t num_cols() {
            return mat[0].size();
        }
        matrix<T> copy() const {
            return matrix<T>(*this);
        }
        T &at(size_t row_index, size_t col_index) {
            if (row_index >= mat.size() || col_index >= mat[0].size()) {
                throw std::invalid_argument("Index out of range.");
            }
            return mat[row_index][col_index];
        }
        const T &at(size_t row_index, size_t col_index) const {
            if (row_index >= mat.size() || col_index >= mat[0].size()) {
                throw std::invalid_argument("Index out of range.");
            }
            return mat[row_index][col_index];
        }
        matrix<T> &make_zero() noexcept {
            for (std::vector<T> &row : mat) {
                for (T &elem : row) {
                    elem = T{0};
                }
            }
            next_row = 0;
            next_col = 0;
            return *this;
        }
        matrix<T> &make_identity() {
            if (!is_square()) {
                throw invalid_matrix_format("Only square matrices can be identity matrices.");
            }
            size_t outer_count = 0;
            size_t inner_count;
            for (std::vector<T> &row : mat) {
                inner_count = 0;
                for (T &elem : row) {
                    elem = inner_count++ == outer_count ? T{1} : T{0};
                }
                ++outer_count;
            }
            return *this;
        }
        bool is_zero() const {
            for (const std::vector<T> &row : mat) {
                for (const T &elem : row) {
                    if (elem != T{0})
                        return false;
                }
            }
            return true;
        }
        bool is_identity() const {
            size_t rows = mat.size();
            size_t cols = mat[0].size();
            size_t one_index = 0;
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    if (i == j) {
                        if (mat[i][j] != T{1}) {
                            return false;
                        }
                    }
                    else {
                        if (mat[i][j] != T{0}) {
                            return false;
                        }
                    }
                }
                ++one_index;
            }
            return true;
        }
        auto begin() const noexcept {
            return mat.cbegin();
        }
        auto end() const noexcept {
            return mat.cend();
        }
        auto begin() noexcept {
            return mat.begin();
        }
        auto end() noexcept {
            return mat.end();
        }
        auto cbegin() const noexcept {
            return mat.cbegin();
        }
        auto cend() const noexcept {
            return mat.cend();
        }
        matrix<T> &transform_entries(void (*func)(T &val)) {
            for (std::vector<T> &row : mat) {
                for (T &elem : row) {
                    func(elem);
                }
            }
            return *this;
        }
        matrix<T> &discard_entries(bool (*unary_predicate)(const T &val)) {
            for (std::vector<T> &row : mat) {
                for (T &elem : row) {
                    if (unary_predicate(elem)) {
                        elem = T{0};
                    }
                }
            }
            return *this;
        }
        matrix<T> &keep_entries(bool (*unary_predicate)(const T &val)) {
            for (std::vector<T> &row : mat) {
                for (T &elem : row) {
                    if (!unary_predicate(elem)) {
                        elem = T{0};
                    }
                }
            }
            return *this;
        }
        matrix<T> &set_insertion_start(size_t row_index, size_t col_index) {
            if (row_index >= mat.size() || col_index >= mat[0].size()) {
                throw std::invalid_argument("Index out of range.");
            }
            next_row = row_index;
            next_col = col_index;
            stop = false;
            return *this;
        }
        matrix<T> &rewind_insertion_start() noexcept {
            next_row = 0;
            next_col = 0;
            stop = false;
            return *this;
        }
        matrix<T> &set_auto_stop(bool auto_stop_insertion = true) noexcept {
            auto_stop = auto_stop_insertion;
            return *this;
        }
        String str() const { // returns a String representation of the matrix
            std::ostringstream out;
            std::streampos before;
            std::streamoff chars_written;
            std::streamoff max_chars_written = 0;
            size_t diff;
            for (const std::vector<T> &row : mat) {
                for (const T &elem : row) {
                    before = out.tellp();
                    out << +elem << ' ';
                    chars_written = out.tellp() - before;
                    if (chars_written > max_chars_written) {
                        max_chars_written = chars_written;
                    }
                }
            }
            out.str("");
            out.clear();
            for (const std::vector<T> &row : mat) {
                out << "| ";
                for (const T &elem : row) {
                    before = out.tellp();
                    out << +elem << ' ';
                    chars_written = out.tellp() - before;
                    diff = (size_t) (max_chars_written - chars_written);
                    for (size_t i = 0; i < diff; ++i) {
                        out << ' ';
                    }
                }
                out << "|\n";
            }
            return {out.str().c_str()};
        }
        static inline matrix<long double> get_2D_rotation_matrix(const long double &angle_rad = PI) {
            if (angle_rad == PI/2) { // returns an exact matrix for the following 4 cases (no f.p. rounding errors)
                return matrix<long double>(2, 2) << 0 << -1 << 1 << 0;
            }
            if (angle_rad == PI) {
                return matrix<long double>(2, 2) << -1 << 0 << 0 << -1;
            }
            if (angle_rad == 3*PI/4) {
                return matrix<long double>(2, 2) << 0 << 1 << -1 << 0;
            }
            if (angle_rad == 2*PI) {
                return matrix<long double>(2, 2).make_identity();
            }
            return matrix<long double>(2, 2) << cosl(angle_rad) << -sinl(angle_rad)
                                             << sinl(angle_rad) <<  cosl(angle_rad);
        }
        static inline matrix<long double> get_2D_rotation_matrix(const long double &&angle_rad = _PI_) {
            return get_2D_rotation_matrix(angle_rad);
        }
        static inline matrix<long double> get_2D_scale_matrix(const long double &&scale) {
            return matrix<long double>(2, 2) << scale << 0 << 0 << scale;
        }
        static inline matrix<long double> get_2D_scale_matrix(const long double &scale) {
            return matrix<long double>(2, 2) << scale << 0 << 0 << scale;
        }
        static inline matrix<long double> get_3D_rotation_matrix(const long double &angle_rad = PI,
                                                                 char about = 'z') {
            if (!(about >= 'x' && about <= 'z') && !(about >= 'X' && about <= 'Z')) {
                throw invalid_axis_error("Invalid axis specified. Axes are: 'x', 'y' and 'z'.");
            }
            if (about < 91) { // make 'about' lower-case in case the char is upper-case
                about += 32;
            }
            if (about == 'x') {
                if (angle_rad == PI/2) { // again, returns an exact matrix for the following 4 cases (no rounding error)
                    return matrix<long double>(3, 3) << 1 << 0 <<  0
                                                     << 0 << 0 << -1
                                                     << 0 << 1 <<  0;
                }
                if (angle_rad == PI) {
                    return matrix<long double>(3, 3) << 1 <<  0 <<  0
                                                     << 0 << -1 <<  0
                                                     << 0 <<  0 << -1;
                }
                if (angle_rad == 3*PI/4) {
                    return matrix<long double>(3, 3) << 1 <<  0 << 0
                                                     << 0 <<  0 << 1
                                                     << 0 << -1 << 0;
                }
                if (angle_rad == 2*PI) {
                    return matrix<long double>(3, 3).make_identity();
                }
                return matrix<long double>(3, 3) << 1 << 0               << 0
                                                 << 0 << cosl(angle_rad) << -sinl(angle_rad)
                                                 << 0 << sinl(angle_rad) <<  cosl(angle_rad);
            }
            else if (about == 'y') {
                if (angle_rad == PI/2) {
                    return matrix<long double>(3, 3) <<  0 << 0 << 1
                                                     <<  0 << 1 << 0
                                                     << -1 << 0 << 0;
                }
                if (angle_rad == PI) {
                    return matrix<long double>(3, 3) << -1 << 0 <<  0
                                                     <<  0 << 1 <<  0
                                                     <<  0 << 0 << -1;
                }
                if (angle_rad == 3*PI/4) {
                    return matrix<long double>(3, 3) << 0 << 0 << -1
                                                     << 0 << 1 <<  0
                                                     << 1 << 0 <<  0;
                }
                if (angle_rad == 2*PI) {
                    return matrix<long double>(3, 3).make_identity();
                }
                return matrix<long double>(3, 3) <<  cosl(angle_rad) << 0 << sinl(angle_rad)
                                                 <<  0               << 1 << 0
                                                 << -sinl(angle_rad) << 0 << cosl(angle_rad);
            }
            if (angle_rad == PI/2) { // returns an exact matrix for the following 4 cases (no f.p. rounding errors)
                return matrix<long double>(3, 3) << 0 << -1 << 0
                                                 << 1 <<  0 << 0
                                                 << 0 <<  0 << 1;
            }
            if (angle_rad == PI) {
                return matrix<long double>(3, 3) << -1 <<  0 << 0
                                                 <<  0 << -1 << 0
                                                 <<  0 <<  0 << 1;
            }
            if (angle_rad == 3*PI/4) {
                return matrix<long double>(3, 3) <<  0 << 1 << 0
                                                 << -1 << 0 << 0
                                                 <<  0 << 0 << 1;
            }
            if (angle_rad == 2*PI) {
                return matrix<long double>(3, 3).make_identity();
            }
            return matrix<long double>(3, 3) << cosl(angle_rad) << -sinl(angle_rad) << 0
                                             << sinl(angle_rad) <<  cosl(angle_rad) << 0
                                             << 0               <<  0               << 1;
        }
        static inline matrix<T> get_3D_rotation_matrix(const long double &&angle_rad = _PI_, char about = 'z') {
            return get_3D_rotation_matrix(angle_rad, about);
        }
        std::vector<T> operator[](size_t index) const { // using [] returns copy of row, so cannot make changes to mat.
            if (index >= mat.size()) {
                throw std::invalid_argument("Index out of range.");
            }
            return {mat[index]};
        }
        template <typename U> requires isConvertible<U, T>
        matrix<T> &operator<<(U value) { // inserts a single value into the first element that is zero in the matrix
            if (auto_stop && stop) { // allows operations such as: matrix << 4 << 2 << 90 << 65 << 24;
                return *this;
            }
            mat[next_row][next_col] = value;
            next_col = next_col == mat[0].size() - 1 ? ++next_row, 0 : next_col + 1;
            if (next_row == mat.size()) {
                if (auto_stop) {
                    stop = true;
                }
                next_row = 0;
            }
            return *this;
        }
        template <typename U> requires isConvertible<U, T>
        matrix<T> &operator=(const matrix<U> &other) {
            if (&other == this) {
                return *this;
            }
            if constexpr (std::same_as<U, T>) {
                this->mat = other.mat;
            }
            else {
                mat.clear();
                std::vector<T> sub;
                for (const std::vector<U> &row : other) {
                    sub.clear();
                    for (const U &elem : row) {
                        sub.push_back(elem);
                    }
                    mat.push_back(sub);
                }
            }
            this->next_row = other.next_row;
            this->next_col = other.next_col;
            this->auto_stop = other.auto_stop;
            this->stop = other.stop;
            return *this;
        }
        template <typename U> requires isConvertible<U, T>
        matrix<T> &operator=(matrix<U> &&other) {
            if (&other == this) {
                return *this;
            }
            if constexpr (std::same_as<U, T>) {
                this->mat = std::move(other.mat);
            }
            else {
                mat.clear();
                std::vector<T> sub;
                for (const std::vector<U> &row : other) {
                    sub.clear();
                    for (const U &elem : row) {
                        sub.push_back(elem);
                    }
                    mat.push_back(sub);
                }
            }
            this->next_row = other.next_row;
            this->next_col = other.next_col;
            this->auto_stop = other.auto_stop;
            this->stop = other.stop;
            return *this;
        }
        template <isNumWrapper U>
        friend std::ostream &operator<<(std::ostream &out, const matrix<U> &mat);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(U s, const matrix<V> &m2) -> matrix<decltype(std::declval<U>()*std::declval<V>())>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m1, const matrix<V> &m2) ->
        matrix<decltype(std::declval<U>()*std::declval<V>() + std::declval<U>()*std::declval<V>())>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const matrix<U> &m1, const matrix<V> &m2) ->
        matrix<decltype(std::declval<U>() + std::declval<V>())>;
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator==(const matrix<U> &m1, const matrix<V> &m2);
        template <isNumWrapper U, isNumWrapper V>
        friend bool operator!=(const matrix<U> &m1, const matrix<V> &m2);
        template <isNumWrapper U>
        friend class matrix;
        template <isNumWrapper U>
        friend class vector;
        template <isNumWrapper U>
        friend class vector2D;
        template <isNumWrapper U>
        friend class vector3D;
    };
    template <isNumWrapper U>
    std::ostream &operator<<(std::ostream &out, const matrix<U> &mat) {
        return out << mat.str();
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &m1, const matrix<V> &m2) ->
    matrix<decltype(std::declval<U>()*std::declval<V>() + std::declval<U>()*std::declval<V>())> {
        using T = decltype(std::declval<U>()*std::declval<V>() + std::declval<U>()*std::declval<V>());
        std::vector<std::vector<T>> r;
        const size_t rows1 = m1.mat.size();
        const size_t columns1 = m1.mat[0].size();
        const size_t rows2 = m2.mat.size();
        const size_t columns2 = m2.mat[0].size();
        for (const auto &[e1, e2] : zip_cref(m1.mat, m2.mat)) { // probing for valid matrix
            if (e1.size() != columns1 || e1.size() != rows2 || e2.size() != columns2) {
                throw std::invalid_argument("For matrix multiplication to be possible, the number of columns in the\n"
                                            "first matrix must equal the number of rows in the second.");
            }
        }
        std::vector<T> sub;
        T total;
        for (size_t i = 0; i < rows1; ++i) {
            sub.clear();
            for (size_t j = 0; j < columns2; ++j) {
                total = 0;
                for (size_t k = 0; k < columns1; ++k) {
                    total += m1.mat[i][k]*m2.mat[k][j]; // faster than calling operator[]
                }
                sub.push_back(total);
            }
            r.push_back(sub);
        }
        return matrix<T>(std::move(r));
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(U s, const matrix<V> &m) -> matrix<decltype(std::declval<U>()*std::declval<V>())> {
        size_t rows = m.mat.size();
        size_t cols = m.mat[0].size();
        matrix<decltype(std::declval<U>()*std::declval<V>())> ret_mat(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                ret_mat.mat[i][j] = s*m.mat[i][j]; // faster than calling << (although it would have looked prettier)
            }
        }
        return ret_mat;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const matrix<U> &m1, const matrix<V> &m2) ->
    matrix<decltype(std::declval<U>() + std::declval<V>())> {
        size_t rows = m1.mat.size();
        size_t cols = m1.mat[0].size();
        using T = decltype(std::declval<U>() + std::declval<V>());
        if (rows != m2.mat.size() || cols != m2.mat[0].size()) {
            throw std::invalid_argument("Matrix dimensions must be equal for matrix addition to be possible.");
        }
        matrix<T> ret_mat(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                ret_mat.mat[i][j] = m1.mat[i][j] + m2.mat[i][j];
            }
        }
        return ret_mat;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator==(const matrix<U> &m1, const matrix<V> &m2) {
        if (m1.mat.size() != m2.mat.size() || m1.mat[0].size() != m2.mat[0].size()) {
            return false;
        }
        for (const auto &[v1, v2] : zip_cref(m1.mat, m2.mat)) {
            for (const auto &[el1, el2] : zip_cref(v1, v2)) {
                if (el1 != el2) {
                    return false;
                }
            }
        }
        return true;
    }
    template <isNumWrapper U, isNumWrapper V>
    bool operator!=(const matrix<U> &m1, const matrix<V> &m2) {
        return !(m1 == m2);
    }
}
#endif
