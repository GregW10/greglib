#include "gregffnn.hpp"

#define LR (1.0l/powl(2, 24))

#define NUM_EXAMPLES 10'000
#define BATCH_SIZE 10'000
#define NUM_PASSES 1'000

static_assert(BATCH_SIZE <= NUM_EXAMPLES, "Batch size greater than number of examples.\n");

#define INPUT_DIM 10
#define OUTPUT_DIM 10

int main() {
    gml::vector<long double> input((uint64_t) INPUT_DIM);
    gml::vector<long double> output((uint64_t) OUTPUT_DIM);
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<long double> dist{-1, 1};
    gml::matrix<long double> mat{OUTPUT_DIM, INPUT_DIM};
    for (long double &value : mat)
        value = dist(rng);
    uint64_t counter = 0;
    std::vector<std::pair<gml::vector<long double>, gml::vector<long double>>> pairs;
    pairs.reserve(NUM_EXAMPLES);
    while (counter++ < NUM_EXAMPLES) {
        uint64_t _i = INPUT_DIM;
        while (_i --> 0)
            input[_i] = dist(rng); // initialise input to random numbers
        output = input;
        output.apply(mat); // apply linear transformation
        output.for_each(sinl); // follow with non-linearity
        pairs.emplace_back(input, output);
        // std::cout << "Sin(" << pairs.back().first[0] << ") = " << pairs.back().second[0] << std::endl;
        std::cout << "\rGenerated " << counter << '/' << NUM_EXAMPLES << " input/output pairs.";
    }
    putchar('\n');
    gml::ffnn<long double> ffnn{};
    ffnn.emplace_back(INPUT_DIM, 128, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);
    ffnn.emplace_back(128, 128, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);
    ffnn.emplace_back(128, 128, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);
    // ffnn.emplace_back(256, 128, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);
    ffnn.emplace_back(128, OUTPUT_DIM, nullptr, nullptr, gml::GLOROT_UNIFORM);
    long double loss;
    std::vector<long double> losses;
    losses.reserve(NUM_EXAMPLES);
    uint64_t passes = 0;
    uint64_t batches;
    while (passes++ < NUM_PASSES) {
        counter = 0;
        while (counter < NUM_EXAMPLES) {
            batches = 0;
            while (batches++ < BATCH_SIZE) {
                auto &[_in, _out] = pairs[counter++];
                ffnn.forward_pass(_in);
                ffnn.backward_pass(_out);
            }
            // std::cout << "\rIteration " << counter << '/' << NUM_EXAMPLES;
            losses.push_back(ffnn.update_params(LR));
        }
        // putchar('\n');
        std::cout << "\rPass " << passes << '/' << NUM_PASSES << " completed";
    }
    putchar('\n');
    std::ofstream out{"losses.dat", std::ios_base::out | std::ios_base::binary};
    if (!out.good()) {
        std::cerr << "Couldn't write data file.\n";
        return 1;
    }
    out.write((char *) losses.data(), sizeof(long double)*losses.size());
    out.close();
    std::cout << "Starting loss: " << losses.front() << "\nFinal loss: " << losses.back() << std::endl;
    // for (const auto &layer : ffnn)
    //     std::cout << layer.biases() << std::endl;
    long double val = -1;
    long double step = 0.1;
    counter = 20;
    // while (counter --> 0) {
    //     input[0] = val;
    //     std::cout << "sin(" << val << "): " << ffnn.fpass_eval(input)[0] << std::endl;
    //     val += step;
    // }
    int _i = INPUT_DIM;
    while (_i --> 0)
        input[_i] = dist(rng); // initialise input to random numbers
    output = input;
    std::cout << "Input:\n" << input << std::endl;
    output.apply(mat);
    output.for_each(sinl);
    std::cout << "True output:\n" << output << std::endl;
    std::cout << "Predicted output:\n" << ffnn.fpass_eval(input) << std::endl;
    return 0;
}
