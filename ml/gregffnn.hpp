#ifndef GREGFFNN_HPP
#define GREGFFNN_HPP

#include "gregvct.hpp"
#include <cmath>
#include <random>

namespace gml {
    namespace activations {
        template <Numeric T>
        void sigmoid(T &val) {
            val = 1/(1 + expl(-val)); // safe even in case of underflow or overflow
        }
        template <Numeric T>
        T softsign(T &val) {
            val = val/(1 + (val < 0 ? -val : val));
        }
        template <Numeric T>
        T softroot(T &val) {
            val = val/(1 + sqrtl(val < 0 ? -val : val)); // my own invention
        }
    }
    namespace losses {
        template <Numeric T>
        T sqloss(const vector<T> &_y, const vector<T> &_f) {
            auto mag = (_y - _f).magnitude(); // consider doing this manually to avoid the construction of another vec.
            return mag*mag;
        }
    }
    template <Numeric T>
    class ffnn {
        /* Feed-Forward Neural Network. */
    public:
        enum initialiser {
            ZEROS, ONES, GLOROT_UNIFORM, GLOROT_NORMAL
        };
        class layer {
            uint64_t DN; // dimensionality of output vector of layer - NOT INCLUDING BIAS
            uint64_t DNm1; // dimensionality of input vector to layer - INCLUDING BIAS
            matrix<T> W{DN, DNm1}; // matrix of weights for layer - includes biases
            vector<T> _hNm1{DNm1}; // input vector to the layer (\(\mathbf{h}^{N-1}\) in my notes)
            vector<T> _gN{DN}; // intermediate vector (after matmul.) (\(\mathbf{g}^N\) in my notes)
            vector<T> _hN{DN}; // output vector of the layer (\(\mathbf{h}^N\) in my notes)
            void (*_sigma)(T&){}; // activation function - can be `nullptr` to signify no activation
            void (*_sigmap)(T&){}; // derivative of activation function - can be `nullptr` to signify no activation
            vector<T> dLdhN{DN}; // gradient of loss w.r.t. outputs of layer - only this is passed down through layers
            matrix<T> dLdW{DN, DNm1}; // gradient of loss w.r.t. weights, stored so can updated later
            layer *prev{}; // pointer to previous layer
            layer *next{}; // pointer to next layer
            void init_gluni() requires (std::is_floating_point_v<T>) { // initialise weights to a uniform Glorot dist.
                std::mt19937_64 rng{std::random_device{}()};
                T bound = sqrtl(6/(DNm1 + DN));
                std::uniform_real_distribution<long double> dist{-bound, bound};
                T *ptr = W.data;
            }
            void init_glono() { // initialise weights to a normal Glorot distribution

            }
        public:
            layer(uint64_t _input_dim, uint64_t _output_dim, initialiser _init = GLOROT_UNIFORM) :
            DNm1{_input_dim}, DN{_output_dim} {
                switch (_init) {
                    case ZEROS:
                        return;
                    case ONES:
                        W.assign(1);
                        break;
                    case GLOROT_UNIFORM:
                        init_gluni();
                        break;
                    case GLOROT_NORMAL:
                        init_glono();
                        break;
                }
            }
            friend class ffnn;
        };
    private:
        std::vector<layer> layers{};
    };
}

#endif
