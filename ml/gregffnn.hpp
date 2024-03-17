#ifndef GREGFFNN_HPP
#define GREGFFNN_HPP

#include "gregvct.hpp"
#include <cmath>
#include <random>
#include <deque>
#include <map>

namespace gml {
    namespace exceptions {
        class invalid_nn_sequence : public std::logic_error {
        public:
            invalid_nn_sequence() :std::logic_error{"Error: invalid sequence of steps in neural network processes.\n"}{}
            explicit invalid_nn_sequence(const char *msg) : std::logic_error{msg} {}
        };
        class invalid_nnw_format : public std::logic_error {
        public:
            invalid_nnw_format() :std::logic_error{"Error: invalid .nnw file format.\n"} {}
            explicit invalid_nnw_format(const char *msg) : std::logic_error{msg} {}
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
        template <Numeric T>
        std::pair<void (*)(T&), void (*)(T&)> get_func_by_id(uint32_t _id) {
            /* Returns the activation function and its derivative with ID `_id`, throws exception if it doesn't exist.*/
            /* An activation function and its derivative are assigned the same ID, as they are always used together. */
            static const std::map<uint32_t, std::pair<void (*)(T&), void (*)(T&)>> _funcs = {
                    {0, {nullptr, nullptr}}, // an ID of zero means no activation function
                    {1, {sigmoid<T>, sigmoid_d<T>}},
                    {2, {softsign<T>, softsign_d<T>}}
            };
            return _funcs.at(_id);
        }
        template <Numeric T>
        uint32_t get_id_by_func(void (*_f)(T&)) {
            /* Returns the ID of activation function `_f`. */
            static const std::map<void (*)(T&), uint32_t> _funcs = {
                    {nullptr, 0},
                    {sigmoid<T>, 1},
                    {softsign<T>, 2}
            };
            return _funcs.at(_f);
        }
        template <Numeric T>
        uint32_t get_id_by_dfunc(void (*_f)(T&)) {
            /* Returns the ID of the activation function derivative `_f`. */
            static const std::map<void (*)(T&), uint32_t> _funcs = {
                    {nullptr, 0},
                    {sigmoid_d<T>, 1},
                    {softsign_d<T>, 2}
            };
            return _funcs.at(_f);
        }
    }
    namespace losses {
        template <Numeric T, Numeric U/*,bool WithGrad = true*/> requires (std::is_convertible<subComT<T, U>, T>::value)
        std::pair<T, vector<T>> sqloss(const vector<U> &_y, const vector<T> &_f) {
            /* Square loss. Returns the loss itself and the derivative of the loss w.r.t. output vector `_f`. */
            vector<T> diff = _f - _y; // this way round because of also returning derivative of loss
            return {diff.mag_sq(), 2*diff};
        }
        // template <Numeric T, Numeric U> requires (std::is_convertible<subComT<T, U>, T>::value)
        // T sqloss<T, U, false>(const vector<U> &_y, const vector<T> &_f) {
        //     /* Square loss. Returns only the loss itself. */
        //     return (_y - _f).mag_sq(); // write manually instead to avoid creating unnecessary `vector<T>` object
        // }
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
            layer(T *_dataW, int64_t rankW, uint64_t *_sW, uint64_t _volumeW, bool copy_shapeW,
                  T *_data_b, int64_t rank_b, uint64_t *_s_b, uint64_t _volume_b, bool copy_shape_b,
                  uint64_t _actfunc_id) :
            DN{*_sW}, DNm1{*(_sW + 1)}, W{_dataW, rankW, _sW, _volumeW, copy_shapeW},
            _b{_data_b, rank_b, _s_b, _volume_b, copy_shape_b, true} {
                try {
                    auto [f1, f2] = activations::get_func_by_id<T>(_actfunc_id);
                    _sigma = f1;
                    _sigmap = f2;
                } catch(const std::out_of_range&) {/* `_sigma` and `_sigmap` already `nullptr`, no action required */}
            }
        public:
#pragma pack(push, 1)
            struct info_chunk {
                uint32_t _fid; // activation function ID
                uint64_t _idim; // input dimension of layer
                uint64_t _odim; // output dimension of layer
            };
            struct layer_chunk {
                T *_weights;
                T *_biases;
            };
#pragma pack(pop)
            class ctor_helper {
                T *_dataW{};
                int64_t rankW{};
                uint64_t *_sW{};
                uint64_t _volumeW{};
                bool copy_shapeW{};
                T *_data_b{};
                int64_t rank_b{};
                uint64_t *_s_b{};
                uint64_t _volume_b{};
                bool copy_shape_b{};
                uint64_t _actfunc_id{};
                ctor_helper(T *_dataW_p,
                            int64_t rankW_p,
                            uint64_t *_sW_p,
                            uint64_t _volumeW_p,
                            bool copy_shapeW_p,
                            T *_data_b_p,
                            int64_t rank_b_p,
                            uint64_t *_s_b_p,
                            uint64_t _volume_b_p,
                            bool copy_shape_b_p,
                            uint64_t _actfunc_id_p) :_dataW{_dataW_p}, rankW{rankW_p}, _sW{_sW_p}, _volumeW{_volumeW_p},
                            copy_shapeW{copy_shapeW_p}, _data_b{_data_b_p}, rank_b{rank_b_p}, _s_b{_s_b_p},
                            _volume_b{_volume_b_p}, copy_shape_b{copy_shape_b_p}, _actfunc_id{_actfunc_id_p} {}
            public:
                friend class ffnn;
            };
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
            layer(const ctor_helper &_ch) :
            DN{*_ch._sW}, DNm1{*(_ch._sW + 1)}, W{_ch._dataW, _ch.rankW, _ch._sW, _ch._volumeW, _ch.copy_shapeW},
            _b{_ch._data_b, _ch.rank_b, _ch._s_b, _ch._volume_b, _ch.copy_shape_b, true} {
                try {
                    auto [f1, f2] = activations::get_func_by_id<T>(_ch._actfunc_id);
                    _sigma = f1;
                    _sigmap = f2;
                } catch(const std::out_of_range&) {/* `_sigma` and `_sigmap` already `nullptr`, no action required */}
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
                // _input = matcvecmul(W, _input) + _b;
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
            // template <Numeric>
            friend class ffnn;
        };
        using ctorh_t = typename layer::ctor_helper;
    private:
        std::deque<layer> layers{};
        std::pair<T, vector<T>> (*_lossf)(const vector<T>&, const vector<T>&) = losses::sqloss<T, T>;
        T _loss{}; // cumulative loss
        vector<T> _foutput{}; // to store the final output of the network after a forward pass
        uint64_t _bsize = 0; // counter for number of forwards/backwards passes performed ( == num. batches)
        bool fpassed = false;
        // bool bpassed = false;
        void load_nnw(const char *path) {
            struct stat buff{};
            if (stat(path, &buff) == -1)
                throw std::ios_base::failure{"Error: could not obtain .nnw file info.\n"};
            if (!S_ISREG(buff.st_mode))
                throw std::ios_base::failure{"Error: .nnw file is not a regular file.\n"};
            if (buff.st_size < sizeof(nnw_header))
                throw exceptions::invalid_nnw_format{"Error: invalid .nnw format - insufficient file size.\n"};
            std::ifstream in{path, std::ios_base::in | std::ios_base::binary};
            if (!in.good())
                throw std::ios_base::failure{"Error: could not open .nnw file for reading.\n"};
            nnw_header _hdr;
            in.read((char *) &_hdr, sizeof(nnw_header));
            if (_hdr.hd[0] != 'N' || _hdr.hd[1] != 'N' || _hdr.hd[2] != 'W')
                throw exceptions::invalid_nnw_format{"Error: invalid .nnw format - invalid file signature.\n"};
            if (_hdr._sizeof_T != sizeof(T))
                throw exceptions::invalid_nnw_format{"Error: invalid .nnw format - reported \"sizeof(T)\" does not "
                                                     "match actual \"sizeof(T)\".\n"};
            if (!_hdr._numl) {
                if (buff.st_size > sizeof(nnw_header))
                    throw exceptions::invalid_nnw_format{"Error: invalid .nnw format - file too large for storing "
                                                         "network with zero layers.\n"};
                in.close();
                return; // an empty .nnw file is a no-op, it is fine
            }
            uint64_t _ichunks_size = _hdr._numl*sizeof(ichunk_t); // size of all layer info chunks
            uint64_t _cum_size = sizeof(nnw_header) + _ichunks_size; // cumulative size
            if (buff.st_size < _cum_size)
                throw exceptions::invalid_nnw_format{"Error: invalid .nnw format - file size insufficient for reported "
                                                     "number of layers.\n"};
            ichunk_t *_ichunks = new ichunk_t[_hdr._numl];
            ichunk_t *iptr = _ichunks;
            in.read((char *) _ichunks, _ichunks_size);
            uint64_t _num_elems = 0;
            uint64_t counter = _hdr._numl;
            while (counter --> 0) {
                _num_elems += iptr->_odim*iptr->_idim + iptr->_odim;
                ++iptr;
            }
            _cum_size += _num_elems*sizeof(T); // by here, `_cum_size` equals the correct total file size
            if (buff.st_size != _cum_size)
                throw exceptions::invalid_nnw_format{"Error: invalid .nnw format - invalid file size.\n"};
            T *_dataW;
            T *_data_b;
            uint64_t _wvol;
            uint64_t _wshape[2];
            uint64_t _bshape[2];
            iptr = _ichunks;
            counter = 0;
            typename std::deque<layer>::iterator _it;
            layer *_last;
            layer *_pn;
            while (counter++ < _hdr._numl) {
                _dataW = new T[(_wvol = iptr->_odim*iptr->_idim)];
                _data_b = new T[iptr->_odim];
                in.read((char *) _dataW, _wvol*sizeof(T));
                in.read((char *) _data_b, iptr->_odim*sizeof(T));
                _wshape[0] = iptr->_odim;
                _wshape[1] = iptr->_idim;
                _bshape[0] = iptr->_odim;
                _bshape[1] = 1;
                /* this->layers.emplace_back(_dataW, 2, _wshape, _wvol, true, _data_b, 2, _bshape, iptr->_odim, true,
                                          iptr->_fid); */
                this->layers.emplace_back(ctorh_t{_dataW, 2, _wshape, _wvol, true, _data_b, 2, _bshape,
                                                  iptr->_odim, true, iptr->_fid});
                if (counter >= 2) {
                    _last = &(*(--(_it = this->layers.end())));
                    --_it;
                    _pn = &(*_it);
                    _last->prev = _pn;
                    _pn->next = _last;
                }
                ++iptr;
            }
            delete [] _ichunks;
            in.close();
        }
    public:
        using ichunk_t = typename layer::info_chunk;
        using lchunk_t = typename layer::layer_chunk;
#pragma pack(push, 1)
        struct nnw_header {
            const char hd[3] = {'N', 'N', 'W'};
            const unsigned char _sizeof_T = sizeof(T);
            uint64_t _numl{}; // number of layers
        };
#pragma pack(pop)
        using deq_size_t = typename std::deque<layer>::size_type;
        ffnn() = default;
        explicit ffnn(std::pair<T, vector<T>> (*loss_function)(const vector<T>&, const vector<T>&)) {
            if (!loss_function)
                throw std::invalid_argument{"Error: loss function passed cannot be nullptr.\n"};
            _lossf = loss_function;
        }
        explicit ffnn(const char *nnw_path) {
            if (!nnw_path)
                throw std::invalid_argument{"Error: path to .nnw file cannot be `nullptr`.\n"};
            this->load_nnw(nnw_path);
        }
        void load_model(const char *nnw_path) {
            if (!nnw_path)
                return;
            this->layers.clear();
            this->load_nnw(nnw_path);
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
        std::pair<T, vector<T>> (*loss_func())(const vector<T>&, const vector<T>&) {
            return this->_lossf;
        }
        std::pair<T, vector<T>> (*loss_func(std::pair<T, vector<T>> (*_new_lossf)(const vector<T>&, const vector<T>&)))
                                (const vector<T>&, const vector<T>&) {
            /* Sets the loss function to `_new_lossf` and returns the previous loss function. */
            if (!_new_lossf)
                return nullptr;
            auto _oldfunc = this->_lossf;
            this->_lossf = _new_lossf;
            return _oldfunc;
        }
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
            // std::cout << "Retvec:\n" << retvec << std::endl;
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
        uint64_t nnw_fsize() const noexcept {
            /* Returns the size a .nnw file would need to be to represent the entire network. */
            uint64_t _tot_params = 0;
            for (const layer &_layer : this->layers)
                _tot_params += _layer.wvol + _layer.DN;
            return 12 + this->layers.size()*sizeof(ichunk_t) + sizeof(T)*_tot_params;
        }
        [[nodiscard("Returns pointer to dynamically allocated memory.\n")]] char* to_nnw() const {
            std::unique_ptr<char[]> _uptr{gen::now_str("NN_weights_", ".nnw")}; // in case exception is thrown below
            this->to_nnw(_uptr.get());
            // delete [] path;
            // return _pos;
            return _uptr.release(); // releases ownership and returns pointer to path, caller must free the memory
        }
        const char* to_nnw(const char *path) const {
            /* Writes the entire network's weights, biases and activation function IDs to a .nnw file. */
            if (!path)
                throw std::invalid_argument{"Error: path to .nnw passed cannot be nullptr.\n"};
            std::ofstream out{path, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc};
            if (!out.good())
                throw std::ios_base::failure{"Error: could not open .nnw file.\n"};
            nnw_header _hdr;
            _hdr._numl = this->layers.size();
            out.write((char *) &_hdr, sizeof(nnw_header));
            if (!_hdr._numl) {
                // std::ofstream::pos_type _pos = out.tellp();
                out.close();
                // return _pos; // not much point in writing an empty .nnw file, but I don't see any reason to prohibit it
                return path;
            }
            ichunk_t *_ichunks = new ichunk_t[_hdr._numl]; // have chosen to allocate to avoid more I/O calls
            ichunk_t *iptr = _ichunks;
            for (const layer &_layer : this->layers) {
                iptr->_idim = _layer.DNm1;
                iptr->_odim = _layer.DN;
                if (!_layer._sigma) {
                    iptr++->_fid = 0;
                    continue;
                }
                try {
                    iptr->_fid = activations::get_id_by_func<T>(_layer._sigma);
                } catch(const std::out_of_range&) {
                    iptr->_fid = -1; // 4294967295 due to underflow - represents unknown function (not one of mine)
                }
                ++iptr;
            }
            out.write((char*) _ichunks, _hdr._numl*sizeof(ichunk_t));
            delete [] _ichunks;
            for (const layer &_layer : this->layers) { // can't avoid the second loop
                out.write((char *) _layer.W.data, _layer.wvol*sizeof(T)); // write layer weights
                out.write((char *) _layer._b.data, _layer.DN*sizeof(T)); // write layer biases
            }
            // std::ofstream::pos_type _pos = out.tellp();
            out.close();
            // return _pos;
            return path;
        }
        const layer &operator[](uint64_t index) const {
            return this->layers[index];
        }
    };
}
#endif
