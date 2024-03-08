#include "gregffnn.hpp"
#include "../misc/gregparse.hpp"

#define LR (1.0l/powl(2, 10))

#define NUM_EXAMPLES 4
#define BATCH_SIZE 4
#define NUM_PASSES ((uint64_t) 1'000'000)

static_assert(BATCH_SIZE <= NUM_EXAMPLES, "Batch size greater than number of examples.\n");

#define INPUT_DIM 2
#define OUTPUT_DIM 1

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    uint64_t num_passes = parser.get_arg("--num", NUM_PASSES);
    gml::vector<long double> input((uint64_t) INPUT_DIM);
    gml::vector<long double> output((uint64_t) OUTPUT_DIM);
    uint64_t counter = 0;
    std::vector<std::pair<gml::vector<long double>, gml::vector<long double>>> pairs =
            {{{{0}, {0}}, {{0}}},
             {{{0}, {1}}, {{1}}},
             {{{1}, {0}}, {{1}}},
             {{{1}, {1}}, {{0}}}};
    gml::ffnn<long double> ffnn{};
    ffnn.emplace_back(INPUT_DIM, 8, gml::activations::softsign<long double>, gml::activations::softsign_d<long double>, gml::GLOROT_UNIFORM);
    // ffnn.emplace_back(4, 16, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);
    // ffnn.emplace_back(16, 64, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);
    // ffnn.emplace_back(64, 16, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);
    // ffnn.emplace_back(16, 4, gml::activations::sigmoid<long double>, gml::activations::sigmoid_d<long double>, gml::GLOROT_UNIFORM);
    ffnn.emplace_back(8, OUTPUT_DIM, nullptr, nullptr, gml::GLOROT_UNIFORM);
    long double loss;
    std::vector<long double> losses;
    losses.reserve(NUM_EXAMPLES);
    uint64_t passes = 0;
    uint64_t batches;
    while (passes++ < num_passes) {
        counter = 0;
        while (counter < NUM_EXAMPLES) {
            batches = 0;
            while (batches++ < BATCH_SIZE) {
                auto &[_in, _out] = pairs[counter++];
                ffnn.forward_pass(_in);
                loss = ffnn.backward_pass(_out);
            }
            ffnn.update_params(LR);
            // std::cout << "\rIteration " << counter << '/' << NUM_EXAMPLES;
            losses.push_back(loss);
        }
        // putchar('\n');
        std::cout << "\rPass " << passes << '/' << num_passes << " completed";
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
    std::cout << "XOR " << pairs[0].first[0] << ", " << pairs[0].first[1] << ": " << ffnn.fpass_eval(pairs[0].first)[0] << std::endl;
    std::cout << "XOR " << pairs[1].first[0] << ", " << pairs[1].first[1] << ": " << ffnn.fpass_eval(pairs[1].first)[0] << std::endl;
    std::cout << "XOR " << pairs[2].first[0] << ", " << pairs[2].first[1] << ": " << ffnn.fpass_eval(pairs[2].first)[0] << std::endl;
    std::cout << "XOR " << pairs[3].first[0] << ", " << pairs[3].first[1] << ": " << ffnn.fpass_eval(pairs[3].first)[0] << std::endl;
    return 0;
}
