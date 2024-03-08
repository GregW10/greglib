#ifndef GREGFFNN_HPP
#define GREGFFNN_HPP

#include "gregvct.hpp"
#include <cmath>
#include <random>
#include <deque>

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
            val = 1/(denom_root*denom_root);
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
            return {diff.mag_sq(), 2*diff};
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
            matrix<T> dLdW{DN, DNm1}; // cumulative gradient of loss w.r.t. weights, stored so can be updated later
            vector<T> dLdb{DN}; // cumulative gradient of loss w.r.t. biases, stored for updating NN later
            vector<T> *_trvec{}; // "travelling" vector - is a pointer to the vector passing through the entire NN
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
                if (act_func) {
                    if (!act_func_derivative)
                        throw std::invalid_argument{"Error: if the activation function is provided, its derivative "
                                                    "must be provided too.\n"};
                }
                else if (act_func_derivative)
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
            const vector<T> &input() const noexcept {
                /* Returns the input to the layer. Only makes sense to call this once fprop has taken place. */
                return this->_hNm1;
            }
            void input(vector<T> &_input) const noexcept {
                /* Sets the input to the entire network and, thus, can only be called on the first layer of the NN. */
                if (this->prev)
                    throw exceptions::invalid_nn_sequence{"Error: input can only be set for first layer of network.\n"};
                if (fpropped)
                    throw exceptions::invalid_nn_sequence{"Error: input can only be set before forward propagation.\n"};
                _trvec = _input;
            }
            void grad_wrt_foutput(vector<T> &_grad) {
                if (this->next)
                    throw exceptions::invalid_nn_sequence{"Error: the gradient of the loss w.r.t. the final output of "
                                                          "the network can only be set for the final layer.\n"};
                if (!fpropped || bpropped)
                    throw exceptions::invalid_nn_sequence{"Error: the gradient of the loss w.r.t. the final output can "
                                                          "only be set after forward propagation and "
                                                          "before backpropagation.\n"};
                _trvec = _grad;
            }
            vector<T> &fprop() {
                // printf("Prev: %p, this: %p, next: %p\n", prev, this, next);
                if (fpropped)
                    throw exceptions::invalid_nn_sequence{"Error: fprop must be called before bprop and at most once "
                                                          "at one time.\n"};
                if (!_trvec)
                    throw exceptions::invalid_nn_sequence{"Error: fprop cannot be called without setting first layer's "
                                                          "input first.\n"};
                this->_hNm1 = *_trvec;
                _trvec->apply_cv(W);
                *_trvec += _b;
                this->_gN = *_trvec;
                if (_sigma)
                    _trvec->for_each(_sigma);
                this->_hN = *_trvec;
                vector<T> *ptr = _trvec;
                _trvec = nullptr;
                if (next)
                    next->_trvec = ptr;
                fpropped = true;
                bpropped = false;
                return *ptr;
            }
            vector<T> &bprop() {
                /* Back-propagate gradients throughout the layer. The pointer `_trvec` will already be pointing to the
                 * rate of change of the loss w.r.t. the outputs of the current layer. Thus, the current layer is
                 * responsible for computing the gradient of the loss w.r.t. the current layer's weights and biases, and
                 * also for setting the vector pointed to by `_trvec` to the gradient of the loss w.r.t. the outputs of
                 * the previous layer. */
                if (!fpropped)
                    throw exceptions::invalid_nn_sequence{"Error: bprop must only be called after fprop.\n"};
                if (!_trvec)
                    throw exceptions::invalid_nn_sequence{"Error: bprop cannot be called without setting final layer's "
                                                          "gradient w.r.t. final layer output first.\n"};
                vector<T> *ptr = _trvec;
                if (_sigmap) {
                    this->_gN.for_each(_sigmap); // convert `_gN` into the grad. of sigma(`_gN`) w.r.t. `_gN`
                    vector<T> &&dLdgN = tens_ops::hadamard(*_trvec, this->_gN);
                    this->dLdW += dLdgN*_hNm1.transpose();
                    this->dLdb += dLdgN;
                    if (prev) {
                        *_trvec = dLdgN.transpose().apply_rv(W).transpose(); // set to grad. w.r.t. prev. layer output
                        prev->_trvec = ptr;
                        this->_trvec = nullptr;
                    }
                } else {
                    /* Case for layer without activation function, so the pre-activated output is equal to the
                     * activated output, and, therefore, the grad. w.r.t. the pre-activated output is equal to the grad.
                     * w.r.t. the final output of the layer. */
                    this->dLdW += (*_trvec)*_hNm1.transpose();
                    this->dLdb += *_trvec;
                    if (prev) {
                        _trvec->transpose().apply_rv(W).transpose(); // make into grad. w.r.t. prev. layer output
                        prev->_trvec = ptr;
                        this->_trvec = nullptr;
                    }
                }
                fpropped = false;
                bpropped = true;
                return *ptr;
            }
            vector<T> &fprop_eval(vector<T> &_input) const {
                /* Forward-propagate without recording anything. */
                _input.apply_cv(W).operator+=(_b);
                if (_sigma)
                    _input.for_each(_sigma);
                return _input;
            }
            /* vector<T> &fprop(vector<T> &_input) {
                // Forward-propagate the input vector `_input`. This will change `_input` itself and return a
                // reference to it.
                if (fpropped)
                    throw exceptions::invalid_nn_sequence{"Error: fprop must be called before bprop and at most once "
                                                          "at a time.\n"};
                this->_hNm1 = _input; // have to save for backpropagation later
                _input.apply_cv(W); // multiply `W` by `_input` and reassign to `_input`-does dimension checking already
                _input += _b; // add biases to `_input`
                this->_gN = _input; // store pre-activated output for backpropagation later
                if (_sigma)
                    _input.for_each(_sigma);
                this->_hN = _input; // store final output of layer for backpropagation
                fpropped = true; // to indicate that the forward propagation has taken place
                bpropped = false;
                return _input;
            }
            vector<T> bprop(vector<T> &_dLdhN) {
                if (!fpropped)
                    throw exceptions::invalid_nn_sequence{"Error: bprop must only be called after fprop.\n"};
                if (this->_sigmap)
                    this->_gN.for_each(this->_sigmap);
                vector<T> dLdgN = tens_ops::hadamard(_dLdhN, this->_gN);
                this->dLdb += dLdgN;
                this->dLdW += dLdgN*this->_hNm1.transpose();
                fpropped = false;
                bpropped = true;
                return dLdgN.transpose().apply_rv(this->W).transpose(); //returns derivative of loss w.r.t. output of previous layer
            } */
            layer &update(long double _lr, uint64_t batch_size = 1) { // update weights/biases, `_lr` is learning rate
                if (!bpropped)
                    throw exceptions::invalid_nn_sequence{"Error: weights and biases can only be updated after the "
                                                          "forwards and backwards passes.\n"};
                if (!batch_size)
                    throw std::invalid_argument{"Error: batch size cannot be zero.\n"};
                auto _frac = _lr/((long double) batch_size);
                this->W -= _frac*this->dLdW;
                this->_b -= _frac*this->dLdb;
                this->dLdW.to_zeros(); // must reset accumulated gradients to zero
                this->dLdb.to_zeros();
                fpropped = false;
                bpropped = false;
                return *this;
            }
            friend class ffnn;
        };
    private:
        std::deque<layer> layers{};
        std::pair<T, vector<T>> (*_lossf)(const vector<T>&, const vector<T>&) = losses::sqloss;
        T _loss{}; // cumulative loss
        vector<T> _foutput{}; // to store the final output of the network after a forward pass
        uint64_t _bsize = 0; // counter for number of forwards/backwards passes performed ( == num. batches)
        bool fpassed = false;
        // bool bpassed = false;
    public:
        using deq_size_t = typename std::deque<layer>::size_type;
        ffnn() = default;
        explicit ffnn(std::pair<T, vector<T>> (*loss_function)(const vector<T>&, const vector<T>&)) {
            if (!loss_function)
                throw std::invalid_argument{"Error: loss function passed cannot be nullptr.\n"};
            _lossf = loss_function;
        }
        template <typename ...Args>
        const layer &emplace_back(Args&&... args) { // preferred method for adding a layer as it constructs it in place
            this->layers.emplace_back(std::forward<Args>(args)...);
            if (this->layers.size() > 1) {
                auto _it = --this->layers.end();
                layer *_last = &(*_it); // pointer to last layer
                layer *_penultimate = &(*--_it); // pointer to penultimate layer
                _penultimate->next = _last; // make penultimate layer's `next` pointer point to last layer
                _last->prev = _penultimate; // make last layer's `prev` pointer point to penultimate layer
                /*
                layer *ptr = this->layers.data() + this->layers.size() - 2; // point to penultimate layer
                ptr->next = ptr + 1;
                ++ptr;
                ptr->prev = ptr - 1;
                return *ptr; */
                return *_last;
            }
            return this->layers.back();
        }
        /* const layer &emplace_back(const layer &_layer) {
            this->push_back(_layer);
            if (this->layers.size() > 1) {
                auto _it = --this->layers.end();
                layer *_last = &_it--; // pointer to last layer
                layer *_penultimate = &_it; // pointer to penultimate layer
                _penultimate->next = _last; // make penultimate layer's `next` pointer point to last layer
                _last->prev = _penultimate; // make last layer's `prev` pointer point to penultimate layer
                return *_last;
            }
            return this->layers.back();
        } */
        /* const layer &add_layer(layer &&_layer) { // commented out until I have a move ctor in `layer`
            this->layers.emplace_back(std::move(_layer));
            return this->layers.back();
        } */
        deq_size_t num_layers() const noexcept {
            return this->layers.size();
        }
        bool empty() const noexcept {
            return this->layers.empty();
        }
        const vector<T> &forward_pass(const vector<T> &input /* , const vector<T> &true_output */) {
            if (this->layers.empty())
                throw exceptions::invalid_nn_sequence{"Error: the forward pass cannot be computed if no layers are "
                                                      "present in the network.\n"};
            if (fpassed)
                throw exceptions::invalid_nn_sequence{"Error: another forward pass cannot be performed without "
                                                      "performing a backwards pass first.\n"};
            _foutput = input;
            this->layers.front()._trvec = &_foutput; // `_foutput` starts out the same as `input`, but becomes the outp.
            for (layer &_layer : this->layers)
                _layer.fprop();
            // for (layer &_layer : this->layers)
            //     _layer.fprop(_input);
            // _loss = _lossf(true_output, _input); // `_input` has actually been transformed into the output of the nn
            fpassed = true;
            return _foutput;
        }
        T backward_pass(const vector<T> &_observed_output) {
            /* Performs the backwards pass through the entire network. Returns the loss. */
            if (!fpassed)
                throw exceptions::invalid_nn_sequence{"Error: the backward pass can only be performed once the forward "
                                                      "pass has taken place.\n"};
            auto [L, dLdhn] = _lossf(_observed_output, _foutput); // loss and grad. w.r.t. final output of NN
            // vector<T> dLdhn = _loss.second; // gradient of loss w.r.t. output of final layer
            auto _rit = this->layers.rbegin(); // have to traverse layers in reverse order, starting from final layer
            auto _rend = this->layers.rend();
            _rit->_trvec = &dLdhn;
            while (_rit != _rend)
                _rit++->bprop();
            // for (layer &_layer : this->layers)
            //     dLdhn = _layer.bprop(dLdhn);
            // bpassed = true;
            _loss += L;
            fpassed = false;
            ++_bsize;
            return L;
        }
        T update_params(long double _lr) {
            /* Updates parameters in entire NN based upon cumulative gradients computed. Returns mean loss. */
            if (!_bsize || fpassed)
                throw exceptions::invalid_nn_sequence{"Error: parameters can only be updated once at least one "
                                                      "forwards and backwards pass have been completed, with an equal "
                                                      "number of both.\n"};
            for (layer &_layer : this->layers)
                _layer.update(_lr, _bsize);
            // fpassed = false;
            // bpassed = false;
            T mean_loss = _loss/_bsize;
            _loss = 0;
            _bsize = 0;
            return mean_loss;
        }
        vector<T> fpass_eval(const vector<T> &_input) const {
            vector<T> retvec = _input;
            for (const layer &_layer : this->layers)
                _layer.fprop_eval(retvec);
            return retvec;
        }
        typename std::deque<layer>::const_iterator begin() const {
            return this->layers.begin();
        }
        typename std::deque<layer>::const_iterator end() const {
            return this->layers.end();
        }
        const layer &operator[](uint64_t index) const {
            return this->layers[index];
        }
    };
}
#endif
