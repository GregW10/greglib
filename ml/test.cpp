#include <iostream>
#include <concepts>
#include <sstream>
#include <iomanip>

namespace gml {
    template <typename T>
    concept Numeric = requires (T value) {
        T{1};
    };
    class crap {
        int g{};
    public:
        crap() = default;
        crap(int a) : g{a} {}
        template <Numeric T>
        friend class tensor;
    };
    template <Numeric T>
    class tensor {
    protected:
        crap c{};
    public:
        class shape {
        public:
            shape() = default;
            friend std::ostream &operator<<(std::ostream &out, const shape &s);
            template <Numeric U, Numeric V>
            friend bool operator==(const typename tensor<U>::shape &s1, const typename tensor<V>::shape &s2) {
                return true;
            }
        };
        tensor() {
            std::cout << "tensor ctor" << std::endl;
        }
        void other() {
            std::cout << this->c.g << std::endl;
        }
    };
    template <Numeric T>
    class matrix : public tensor<T> {
    public:
        void something() {
            std::cout << this->c.g << std::endl;
        }
    };
    template <Numeric U>
    std::ostream &operator<<(std::ostream &out, const typename tensor<U>::shape &s) {
        return out << "PRINTING AT LAST!!!\n" << std::endl;
    }
}

int main() {
    std::ostringstream oss;
    std::cout << "Beginning: " << oss.tellp() << std::endl;
    oss << "Hello!";
    std::cout << "Now: " << oss.tellp() << std::endl;
    std::cout << "Current contents: " << oss.str() << std::endl;
    oss.seekp(0);
    std::cout << "After seeking: " << oss.tellp() << std::endl;
    std::cout << "Current contents: " << oss.str() << std::endl;
    oss << "Shh";
    std::cout << "Now: " << oss.tellp() << std::endl;
    std::cout << "Current contents: " << oss.str() << std::endl;
    std::cout << "Default width: " << std::cout.width() << ", default precision: " << std::cout.precision() << std::endl;
    std::cout << std::setw(3) << std::setprecision(5);
    std::cout << 98 << ", " << 987 << ", " << 9876 << std::endl;
    std::cout << 987654321 << std::endl;
    std::cout << 1234546.123456l << std::endl;
    gml::matrix<int> m;
    m.other();
    return 0;
}
