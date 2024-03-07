#include "gregffnn.hpp"

enum shite {
    zero, one
};

void func(shite v) {
    std::cout << v << std::endl;
}

int main() {

    gml::ffnn<long double>::layer layer{128, 64, gml::activations::sigmoid, gml::activations::sigmoid};

    gml::vector<long double> vec{128ull};

    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<long double> dist{-1, 1};

    uint64_t counter = 128;

    while (counter --> 0)
        vec[counter] = dist(rng);

    gml::vector<long double> tv = vec;

    tv.for_each(sinl);

    gml::matrix<long double> mat{16, 128};

    for (long double &value : mat)
        value = dist(rng);

    tv.apply(mat);

    std::cout << "vec:\n" << vec << "\ntv:\n" << tv << std::endl;

    gml::ffnn<long double> ffnn{};

    ffnn.emplace_layer(128, 64, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);

    ffnn.emplace_layer(64, 32, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);

    ffnn.emplace_layer(32, 16, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);

    auto l = ffnn.forward_pass(vec, tv);

    std::cout << "l:\n" << l << std::endl;

    ffnn.backward_pass();

    return 0;
}
