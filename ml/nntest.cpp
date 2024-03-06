#include "gregffnn.hpp"

enum shite {
    zero, one
};

void func(shite v) {
    std::cout << v << std::endl;
}

int main() {

    gml::ffnn<long double> nn;

    long double val = 1000000000;

    gml::activations::sigmoid(val);

    std::cout << val << std::endl;

    val = -10000000000;

    gml::activations::sigmoid(val);

    std::cout << val << std::endl;

    func(zero);
    func(one);

    return 0;
}
