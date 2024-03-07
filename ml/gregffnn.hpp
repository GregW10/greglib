#ifndef GREGFFNN_HPP
#define GREGFFNN_HPP

#include "gregvct.hpp"
#include <cmath>
#include <random>

namespace gml {
    namespace exceptions {
        class invalid_nn_sequence : public std::logic_error {
        public:
            invalid_nn_sequence() :std::logic_error{"Error: invalid sequence of steps in neural network processes.\n"}{}
            explicit invalid_nn_sequence(const char *msg) : std::logic_error{msg} {}
        };
    }
    namespace activations {
        template <Numeric T>
        void sigmoid(T &val) { // logistic sigmoid
            val = 1/(1 + expl(-val)); // safe even in case of underflow or overflow
        }
        template <Numeric T>
        void sigmoid_d(T &val) { // logistic sigmoid derivative
            T sig = val;
            sigmoid(sig);
            val = sig*(1 - sig);
        }
        template <Numeric T>
        void softsign(T &val) {
            val /= (1 + (val < 0 ? -val : val));
        }
        template <Numeric T>
        void softsign_d(T &val) { // softsign derivative
            T denom_root = 1 + (val < 0 ? -val : val);
            val = 1/denom_root*denom_root;
        }
        template <Numeric T>
        void softroot(T &val) {
            val = val/(1 + sqrtl(val < 0 ? -val : val)); // my own invention
        }
    }
    namespace losses {
        template <Numeric T>
        std::pair<T, vector<T>> sqloss(const vector<T> &_y, const vector<T> &_f) {
            /* Square loss. Returns the loss itself and the derivative of the loss w.r.t. output vector `_f`. */
            vector<T> diff = _f - _y; // this way round because of also returning derivative of loss
            return {diff.mag_sq(), 2*diff.transpose()};
        }
    }
    enum initialiser {
        ZEROS, ONES, GLOROT_UNIFORM, GLOROT_NORMAL
    };
    template <Numeric T>
    class ffnn {
        /* Feed-Forward Neural Network. */
    public:
        class layer {
            uint64_t DN; // dimensionality of output vector of layer
            uint64_t DNm1; // dimensionality of input vector to layer
            matrix<T> W{DN, DNm1}; // matrix of weights for layer (\(W^n\) in my notes)
            uint64_t wvol = W.vol; // total number of weights in `W` matrix
            vector<T> _b{DN}; // bias vector for layer (\(\mathbf{b}^n\) in my notes)
            vector<T> _hNm1{DNm1}; // input vector to the layer (\(\mathbf{h}^{N-1}\) in my notes)
            vector<T> _gN{DN}; // intermediate vector (after matmul. and adding bias) (\(\mathbf{g}^N\) in my notes)
            vector<T> _hN{DN}; // output vector of the layer (\(\mathbf{h}^N\) in my notes)
            void (*_sigma)(T&){}; // activation function - can be `nullptr` to signify no activation
            void (*_sigmap)(T&){}; // derivative of activation function - can be `nullptr` to signify no activation
            // vector<T> dLdhN{DN}; // gradient of loss w.r.t. outputs of layer - only this is passed down through layers
            matrix<T> dLdW{DN, DNm1}; // gradient of loss w.r.t. weights, stored so can be updated later
            vector<T> dLdb{DN}; // gradient of loss w.r.t. biases, stored for updating NN later
            layer *prev{}; // pointer to previous layer
            layer *next{}; // pointer to next layer
            bool fpropped = false; // boolean indicating whether forward propagation has taken place
            bool bpropped = false; // to indicate whether backpropagation has taken place
            void init_gluni() requires (std::is_floating_point_v<T>) { // initialise weights to a uniform Glorot dist.
                std::mt19937_64 rng{std::random_device{}()};
                T bound = sqrtl(6.0l/(DNm1 + DN));
                std::uniform_real_distribution<long double> dist{-bound, bound};
                T *ptr = W.data;
                uint64_t counter = wvol;
                while (counter --> 0)
                    *ptr++ = dist(rng);
            }
            void init_glono() { // initialise weights to a normal Glorot distribution
                // will not implement this until I have derived it myself (from Glorot's paper which gave the above)
            }
        public:
            layer(uint64_t _input_dim,
                  uint64_t _output_dim,
                  void (*act_func)(T&) = nullptr,
                  void (*act_func_derivative)(T&) = nullptr,
                  initialiser _init = GLOROT_UNIFORM) :
            DN{_output_dim}, DNm1{_input_dim}, _sigma{act_func}, _sigmap{act_func_derivative} {
                if (act_func)
                    if (!act_func_derivative)
                        throw std::invalid_argument{"Error: if the activation function is provided, its derivative "
                                                    "must be provided too.\n"};
                if (act_func_derivative)
                    if (!act_func)
                        throw std::invalid_argument{"Error: if the derivative of the activation function is provided, "
                                                    "so must be the activation function itself.\n"};
                switch (_init) {
                    case ZEROS:
                        return; // because `W` has already been initialised to zeros
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
            const matrix<T> &weights() const noexcept {
                return this->W;
            }
            const vector<T> &biases() const noexcept {
                return this->_b;
            }
            vector<T> &fprop(vector<T> &_input) {
                if (bpropped)
                    throw exceptions::invalid_nn_sequence{"Error: fprop cannot be called after bprop without updating "
                                                          "the weights and biases first.\n"};
                /* Forward-propagate the input vector `_input`. This will change `_input` itself and return a
                 * reference to it. */
                this->_hNm1 = _input; // have to save for backpropagation later
                _input.apply_cv(W); // multiply `W` by `_input` and reassign to `_input`-does dimension checking already
                _input += _b; // add biases to `_input`
                this->_gN = _input; // store pre-activated output for backpropagation later
                if (_sigma)
                    _input.for_each(_sigma);
                this->_hN = _input; // store final output of layer for backpropagation
                fpropped = true; // to indicate that the forward propagation has taken place
                return _input;
            }
            vector<T> bprop(vector<T> &_dLdhN) {
                if (!fpropped)
                    throw exceptions::invalid_nn_sequence{"Error: bprop must only be called after fprop.\n"};
                if (this->_sigmap)
                    this->_gN.for_each(this->_sigmap);
                vector<T> dLdgN = tens_ops::hadamard(_dLdhN, this->_gN);
                this->dLdb = dLdgN;
                this->dLdW = dLdgN*this->_hNm1.transpose();
                bpropped = true;
                return dLdgN.transpose().apply_rv(this->W); //returns derivative of loss w.r.t. output of previous layer
            }
            layer &update(long double _lr) { // update weights and biases, `_lr` is the learning rate
                if (!bpropped)
                    throw exceptions::invalid_nn_sequence{"Error: weights and biases can only be updated after the "
                                                          "forwards and backwards passes.\n"};
                this->W -= _lr*this->dLdW;
                this->_b -= _lr*this->dLdb;
                fpropped = false;
                bpropped = false;
                return *this;
            }
            friend class ffnn;
        };
    private:
        std::vector<layer> layers{};
        std::pair<T, vector<T>> (*_lossf)(const vector<T>&, const vector<T>&) = losses::sqloss;
        std::pair<T, vector<T>> _loss{};
        bool fpassed = false;
        bool bpassed = false;
    public:
        using vec_size_t = typename std::vector<layer>::size_type;
        ffnn() = default;
        explicit ffnn(std::pair<T, vector<T>> (*loss_function)(const vector<T>&, const vector<T>&)) {
            if (!loss_function)
                throw std::invalid_argument{"Error: loss function passed cannot be nullptr.\n"};
            _lossf = loss_function;
        }
        template <typename ...Args>
        const layer &emplace_layer(Args&&... args) { // preferred method for adding a layer as it constructs it in place
            return this->layers.emplace_back(std::forward<Args>(args)...); // check this
        }
        const layer &add_layer(const layer &_layer) {
            this->push_back(_layer);
            return this->layers.back();
        }
        /* const layer &add_layer(layer &&_layer) { // commented out until I have a move ctor in `layer`
            this->layers.emplace_back(std::move(_layer));
            return this->layers.back();
        } */
        vec_size_t num_layers() const noexcept {
            return this->layers.size();
        }
        bool empty() const noexcept {
            return this->layers.empty();
        }
        vector<T> forward_pass(const vector<T> &input, const vector<T> &true_output) {
            if (this->layers.empty())
                throw exceptions::invalid_nn_sequence{"Error: the forward pass cannot be computed if no layers are "
                                                      "present in the network.\n"};
            vector<T> _input = input;
            for (layer &_layer : this->layers)
                _layer.fprop(_input);
            _loss = _lossf(true_output, _input); // `_input` has actually been transformed into the output of the nn
            fpassed = true;
            return _input;
        }
        void backward_pass() {
            if (!fpassed)
                throw exceptions::invalid_nn_sequence{"Error: the backward pass can only be performed once the forward "
                                                      "pass has taken place.\n"};
            vector<T> dLdhn = _loss.second; // gradient of loss w.r.t. output of a layer
            for (layer &_layer : this->layers)
                dLdhn = _layer.bprop(dLdhn);
            bpassed = true;
        }
        void update_params(long double _lr) {
            if (!bpassed)
                throw exceptions::invalid_nn_sequence{"Error: parameters can only be updated once the forwards and "
                                                      "backwards passes have been completed.\n"};
            for (layer &_layer : this->layers)
                _layer.update(_lr);
            fpassed = false;
            bpassed = false;
        }
        const layer &operator[](uint64_t index) const {
            return this->layers[index];
        }
    };
}
#endif
