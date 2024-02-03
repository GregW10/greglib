#ifndef GREG8TREE_H
#define GREG8TREE_H

#ifndef __cplusplus
#error "greg8tree.hpp is a C++ header file only.\n"
#endif

#include "gregbod.hpp"

namespace gtd {
    class barnes_hut_error : public nbody_error {
    public:
        using nbody_error::nbody_error;
        barnes_hut_error() :
        nbody_error{"An error occurred with the creation and/or usage of a gtd::bh_cube<> object.\n"} {}
    };
    class invalid_body_placement : public barnes_hut_error {
        char msg[100]{};
    public:
        invalid_body_placement() : msg{"Error: body does not fall within the current Barnes-Hut cube.\n"} {}
        invalid_body_placement(const char *message) {
            strcpy_c(this->msg, message);
        }
        invalid_body_placement(uint64_t id) {
            snprintf(msg, 100,"Error: body with ID = %" PRIu64" is not contained within the current Barnes-Hut cube.\n",
                     id);
        }
        const char *what() const noexcept override {
            return msg;
        }
    };
    class indistinguishable_bodies : public barnes_hut_error {
        char msg[176]{};
    public:
        indistinguishable_bodies() : msg{"Error: body positions cannot be distinguished from one another.\n"} {}
        indistinguishable_bodies(const char *message) {
            strcpy_c(this->msg, message);
        }
        indistinguishable_bodies(uint64_t id1, uint64_t id2) {
            snprintf(msg, 176, "Error: body with ID = %" PRIu64" and body with ID = %" PRIu64" have positions that are "
                               "either equal or so close together that they cannot be told apart.\n", id1, id2);
        }
        const char *what() const noexcept override {
            return msg;
        }
    };
    template <isNumWrapper, isNumWrapper, isNumWrapper, uint64_t>
    class bh_cube;
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t rF>
    std::ostream &operator<<(std::ostream&, const bh_cube<M, R, T, rF>&);
    // template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, uint64_t, uint64_t, bool>
    // class system;
    template <isNumWrapper, isNumWrapper, isNumWrapper, uint64_t>
    class bh_tree;
    template <isNumWrapper M = long double, isNumWrapper R = long double,
              isNumWrapper T = long double, uint64_t rF = 0>
    class bh_cube {
        using vec = vector3D<T>;
        using bod_t = body<M, R, T, rF>;
        using cube_t = bh_cube<M, R, T, rF>;
        vec _btm{}; // corner of cube with smallest x, y & z coordinates
        vec _top{}; // corner of cube with largest x, y & z coordinates
        vec _com{}; // centre of mass of all bodies contained within cube
        M _mass{}; // total mass of all bodies contained within cube
#ifndef RUNNING_COM
#ifdef RUNNING_MR_SUM
        // using mr_t = decltype(std::declval<M>()*std::declval<T>());
        vec mr_sum{};
#endif
#endif
        /* Only a leaf node (or cube) will have a _bod that is NOT a nullptr. */
        const bod_t *_bod; // pointer to body stored - only stores if leaf node
        cube_t *parent = nullptr; // parent cube
        /* Perhaps unconventionally, I am indexing the octants as follows: 0 -> (-,-,-), 1 -> (+,-,-),
         * 2 -> (-,+,-), 3 -> (+,+,-), 4 -> (-,-,+), 5 -> (+,-,+), 6 -> (-,+,+), 7 -> (+,+,+), */
        cube_t *_sub[8]{}; // 8 pointers to child cubes
        std::pair<vec, vec> corners[8]; // bottom and top corners of child cubes
#ifdef BH_REC_NUM
        uint64_t _num = 1; // number of bodies in cube - impossible to create empty cube, so initialised to 1
#endif
        void set_corners() {
            /* Method to set the bounds of the 8 child cubes (even if they are never created). */
            T half_x = (_btm.x + _top.x)/2.0l;
            T half_y = (_btm.y + _top.y)/2.0l;
            T half_z = (_btm.z + _top.z)/2.0l;
            std::pair<vec, vec> *ptr = corners;
            ptr->first = _btm; // octant 0 bottom
            ptr->second.x = half_x; // octant 0 top
            ptr->second.y = half_y;
            ptr++->second.z = half_z;
            ptr->first.x = half_x; // octant 1 bottom
            ptr->first.y = _btm.y;
            ptr->first.z = _btm.z;
            ptr->second.x = _top.x; // octant 1 top
            ptr->second.y = half_y;
            ptr++->second.z = half_z;
            ptr->first.x = _btm.x; // octant 2 bottom
            ptr->first.y = half_y;
            ptr->first.z = _btm.z;
            ptr->second.x = half_x; // octant 2 top
            ptr->second.y = _top.y;
            ptr++->second.z = half_z;
            ptr->first.x = half_x; // octant 3 bottom
            ptr->first.y = half_y;
            ptr->first.z = _btm.z;
            ptr->second.x = _top.x; // octant 3 top
            ptr->second.y = _top.y;
            ptr++->second.z = half_z;
            ptr->first.x = _btm.x; // octant 4 bottom
            ptr->first.y = _btm.y;
            ptr->first.z = half_z;
            ptr->second.x = half_x; // octant 4 top
            ptr->second.y = half_y;
            ptr++->second.z = _top.z;
            ptr->first.x = half_x; // octant 5 bottom
            ptr->first.y = _btm.y;
            ptr->first.z = half_z;
            ptr->second.x = _top.x; // octant 5 top
            ptr->second.y = half_y;
            ptr++->second.z = _top.z;
            ptr->first.x = _btm.x; // octant 6 bottom
            ptr->first.y = half_y;
            ptr->first.z = half_z;
            ptr->second.x = half_x; // octant 6 top
            ptr->second.y = _top.y;
            ptr++->second.z = _top.z;
            ptr->first.x = half_x; // octant 7 bottom
            ptr->first.y = half_y;
            ptr->first.z = half_z;
            ptr->second = _top; // octant 7 top
        }
        void check_first_body() {
            if (_bod == nullptr)
                throw std::invalid_argument{"Error: pointer to gtd::body<M, R, T, rF> cannot be nullptr.\n"};
            if (!_bod->curr_pos.between(_btm, _top))
                throw invalid_body_placement{_bod->id};
            this->_mass = this->_bod->mass_;
#ifndef RUNNING_COM
#ifdef RUNNING_MR_SUM
            this->mr_sum = this->_mass*this->_bod->curr_pos;
#endif
#else
            this->_com = this->_bod->curr_pos;
#endif
        }
        void dispatch(const bod_t *_body) {
            std::pair<vec, vec> *ptr;
            cube_t *cptr = this;
            unsigned char i;
            // for (unsigned char i = 0; i < 8; ++i, ++ptr) {
#ifdef RUNNING_COM
            M prev_mass;
#endif
            while (true) {
                ptr = cptr->corners;
                for (i = 0; i < 8;) {
                    if (_body->curr_pos.between(ptr->first, ptr->second)) {
                        // if (_body->id == 0 && i == 7)
                        //     std::cout << "first: " << ptr->first << ", second: " << ptr->second << std::endl;
                        // if (_sub[i] == nullptr) {
                        //     _sub[i] = new bh_cube{this, ptr->first, ptr->second, _body};
                        // } else {
                        //     _sub[i]->add_body(_body);
                        // }
                        // return true;
                        // cptr->tot_mass += body->mass_;
                        if (cptr->_sub[i] == nullptr) {
                            cptr->_sub[i] = new bh_cube{cptr, ptr->first, ptr->second, _body};
                            // std::cout << "Created: " << cptr->_sub[i] << std::endl;
                            if (cptr->_bod != nullptr) { // means *cptr was a leaf node
                                if (_body->curr_pos == cptr->_bod->curr_pos)
                                    throw indistinguishable_bodies{_body->id, cptr->_bod->id};
                                _body = cptr->_bod;
                                cptr->_bod = nullptr;
                                ptr -= i;
                                i = 0;
                                continue;
                            }
                            return;
                        }
                        else {
                            cptr = cptr->_sub[i];
#ifdef RUNNING_COM
                            prev_mass = cptr->_mass;
#endif
                            cptr->_mass += _body->mass_;
#ifndef RUNNING_COM
#ifdef RUNNING_MR_SUM
                            cptr->mr_sum += _body->mass_*_body->curr_pos;
#endif
#else
                            cptr->_com = (cptr->_com*prev_mass + _body->mass_*_body->curr_pos)/cptr->_mass;
#endif
#ifdef BH_REC_NUM
                            ++cptr->_num;
#endif
                            break;
                        }
                    }
                    // std::cout << "i: " << +i << std::endl;
                    ++i;
                    ++ptr;
                }
                if (i == 8) [[unlikely]] {
                    // std::cout << "bottom: " << cptr->_btm << ", top: " << cptr->_top << std::endl;
                    throw invalid_body_placement{_body->id}; // in case of f.p. error
                }
            }
            // throw invalid_body_placement{_body->id}; // in case of f.p. error
        }
#ifndef RUNNING_COM
    public:
        void calc_com() noexcept {
            if (this->_bod != nullptr) {
                this->_com = this->_bod->curr_pos;
                return;
            }
            cube_t *cptr = this;
            cube_t **sptr = this->_sub;
            unsigned char i = 0;
#ifdef RUNNING_MR_SUM
            this->_com = this->mr_sum/this->_mass; // calculate com of this cube outside loop
            while (true) {
                // cptr->_com = cptr->mr_sum/cptr->_mass;
                // sptr = cptr->_sub;
                while (i < 8) {
                    if (*sptr != nullptr) {
                        if ((*sptr)->_bod == nullptr) {
                            cptr = *sptr;
                            sptr = cptr->_sub;
                            i = 0;
                            continue;
                        } else {
                            (*sptr)->_com = (*sptr)->mr_sum/(*sptr)->_mass;
                        }
                    }
                    ++i;
                    ++sptr;
                }
                if (cptr->parent == nullptr)
                    break;
                cptr->_com = cptr->mr_sum/cptr->_mass;
                i = 1;
                sptr = cptr->parent->_sub;
                while (*sptr++ != cptr) ++i;
                cptr = cptr->parent;
            }
#else
            vec sub_mr{};
            while (true) {
                while (i++ < 8) {
                    if (*sptr != nullptr) {
                        if ((*sptr)->_bod == nullptr) {
                            cptr = *sptr;
                            sptr = cptr->_sub;
                            i = 0;
                            continue;
                        } else {
                            (*sptr)->_com = (*sptr)->_bod->curr_pos;
                        }
                    }
                    ++sptr;
                }
                while (i --> 0) {
                    if (*--sptr != nullptr) {
                        sub_mr += (*sptr)->_mass*(*sptr)->_com;
                    }
                }
                cptr->_com = sub_mr/cptr->_mass;
                if (cptr->parent == nullptr)
                    return;
                i = 1;
                sptr = cptr->parent->_sub;
                while (*sptr++ != cptr) ++i;
                cptr = cptr->parent;
                sub_mr.x = T{}; // better not to call make_zero in case compiler doesn't optimise function call away
                sub_mr.y = T{};
                sub_mr.z = T{};
            };
#endif
        }
#endif
    public:
        class nn_iterator; // forward declarations required for in-class friend declarations
        class cell_iterator;
        class body_iterator {
        protected:
            cube_t *_cptr{}; // pointer to current leaf node/cube
            cube_t **_sptr{}; // pointer to pointer to current node within parent node's array
            unsigned char index{}; // position of _cptr within parent node's array
        public:
            body_iterator() = default;
            body_iterator(cube_t *_cube) : _cptr{_cube} {}
            void erase(bool adjust_com = true) { // for the future: make it return an iterator to the next body
                if (_cptr->parent == nullptr) {
                    delete _cptr;
                    return;
                }
                M mass = _cptr->_bod->mass_; // valid body iterators will always be pointing to some body
                M mass_before;
                vec pos = _cptr->_bod->curr_pos;
                vec _bmr = mass*pos;
                // cube_t *_cpy = _cptr;
                // cube_t *_cpy_sub;
                // bool assign_body_to_updated_cube = false; // really went to town with var. names here eh
#ifndef BH_REC_NUM
                uint64_t _num;
#define BH_BODY_NUMBER(...) _num
#define BH_CALC_BODY_NUMBER(cube_token) \
                _num = 0; \
                _sptr = cube_token->_sub; \
                index = 8; \
                while (index --> 0) { \
                    if (*_sptr != nullptr) { \
                        if ((*_sptr)->_bod == nullptr) { /* means a sibling cube has children, so > 1 body in it */ \
                            _num = 3; /* could be any value that is NOT 2 */ \
                            break; \
                        } \
                        ++_num; \
                    } \
                    ++_sptr; \
                }
#else
#define BH_BODY_NUMBER(...) __VA_ARGS__
#define BH_CALC_BODY_NUMBER(...)
#endif
                if (_cptr->parent->parent == nullptr) { // I'd rather deal with this explicitly to avoid everyth. below
                    BH_CALC_BODY_NUMBER(_cptr->parent)
                    if (BH_BODY_NUMBER(_cptr->parent->_num) == 2) { // case for when parent is root and entire tree has
                        _sptr = _cptr->parent->_sub; // only two bodies
                        index = 0;
                        while (index++ < 8) {
                            if (*_sptr != nullptr && *_sptr != _cptr)
                                break;
                            ++_sptr;
                        }
                        _cptr->parent->_bod = (*_sptr)->_bod;
                        delete _cptr;
                        _cptr = (*_sptr)->parent;
                        delete *_sptr;
                        _sptr = _cptr->_sub;
                        index = 8;
                        while (index --> 0) *_sptr++ = nullptr;
                        if (adjust_com) {
                            _cptr->_mass = _cptr->_bod->mass_;
                            _cptr->_com = _cptr->_bod->curr_pos;
                        }
                        return;
                    }
                    _sptr = _cptr->parent->_sub;
                    while (*_sptr != _cptr) ++_sptr;
                    _cptr = _cptr->parent;
                    delete *_sptr;
                    *_sptr = nullptr;
                    if (adjust_com) {
                        mass_before = _cptr->_mass;
                        _cptr->_mass -= mass;
                        _cptr->_com = (_cptr->_com*mass_before - mass*pos)/_cptr->_mass;
                    }
                    return;
                }
                cube_t *to_update = _cptr->parent;
                cube_t *to_delete = _cptr;
                BH_CALC_BODY_NUMBER(to_update)
                if (to_update->parent != nullptr && BH_BODY_NUMBER(to_update->_num) == 2) {
                    _sptr = to_update->_sub;
                    while (true) {
                        if (*_sptr != nullptr && *_sptr != to_delete)
                            break;
                        ++_sptr;
                    }
                    const bod_t *to_assign_to_updated_cube = (*_sptr)->_bod;
                    delete *_sptr;
                    *_sptr = nullptr;
                    cube_t *to_assign_to = to_update;
                    to_update = to_update->parent;
#ifndef BH_REC_NUM
                    _sptr = to_update->_sub;
                    index = 0;
                    while (index++ < 8) {
                        if (*_sptr != nullptr && *_sptr != to_assign_to)
                            goto after_loop;
                        ++_sptr;
                    }
                    // BH_CALC_BODY_NUMBER(to_update)
                    while (to_update->parent != nullptr) {
#else
                    while (to_update->_num == 2 && to_update->parent != nullptr) {
#endif
                        to_delete = to_assign_to;
                        to_assign_to = to_update;
                        to_update = to_update->parent;
#ifndef BH_REC_NUM
                        _sptr = to_update->_sub;
                        index = 0;
                        while (index++ < 8) {
                            if (*_sptr != nullptr && *_sptr != to_assign_to)
                                goto after_loop;
                            ++_sptr;
                        }
#endif
                    }
                    after_loop:
                    // to_delete = to_update;
                    // to_update = to_update->parent;
                    // BH_CALC_BODY_NUMBER
                    /*
                    do {
                        to_assign_to = to_update;
                        to_update = to_update->parent;
                        // assign_body_to_updated_cube = true;
#ifndef BH_REC_NUM
                        if (to_update->parent == nullptr)
                            break;
                        BH_CALC_BODY_NUMBER(to_update->parent)
                    } while (BH_BODY_NUMBER(blah_blah_blah) == 2);
#else
                    } while (to_update->parent != nullptr && BH_BODY_NUMBER == 2);
#endif
                     */
                    to_assign_to->_bod = to_assign_to_updated_cube;
                    to_assign_to->_mass = to_assign_to_updated_cube->mass_;
#ifdef RUNNING_COM
                    to_assign_to->_com = to_assign_to_updated_cube->curr_pos;
#else
#ifdef RUNNING_MR_SUM
                    to_assign_to->mr_sum = to_assign_to_updated_cube->mass_*to_assign_to_updated_cube->curr_pos;
#endif
#endif
                    _sptr = to_assign_to->_sub;
                }
#undef BH_BODY_NUMBER
#undef BH_CALC_BODY_NUMBER
                else {
                    _sptr = to_update->_sub;
                }
                // --to_update->_num;
                while (*_sptr != to_delete) ++_sptr;
                *_sptr = nullptr;
                to_delete->parent = nullptr; // setting the parent to nullptr and then deleting it will destroy all the
                delete to_delete; // cells below it recursively (if present), destroying completely its part of the tree
                // if (to_assign_to_updated_cube != nullptr) {
                //     to_update->_bod = to_assign_to_updated_cube;
                //     // _sptr = to_update->_sub;
                //     // index = 8;
                //     // while (index --> 0) *_sptr++ = nullptr; // likely just as fast as first finding index
                // }
                if (!adjust_com) { // could be for when a tree is going to be destroyed or when COMs are not imp. anymo.
#ifdef BH_REC_NUM
                    do {
                        --to_update->_num;
                    } while ((to_update = to_update->parent) != nullptr);
#endif
                    return;
                }
                do {
                    // _cpy_sub = _cpy;
                    // _cpy = _cpy->parent;
                    // index = 0;
                    // _sptr = _cpy->_sub;
                    // while (*_sptr != _cpy_sub) {
                    //     ++index;
                    //     ++_sptr;
                    // }
                    // delete _cpy_sub;
                    // _cpy->_sub[index] = nullptr;
                    mass_before = to_update->_mass;
                    to_update->_mass -= mass;
                    to_update->_com = (to_update->_com*mass_before - _bmr)/to_update->_mass;
#ifdef RUNNING_MR_SUM
                    to_update->mr_sum -= _bmr;
#endif
#ifdef BH_REC_NUM
                    --to_update->_num;
#endif
                } while ((to_update = to_update->parent) != nullptr);
            }
            const bod_t &operator*() const {
                return *_cptr->_bod;
            }
            const bod_t *operator->() const {
                return _cptr->_bod;
            }
            bod_t &operator*() {
                return const_cast<bod_t&>(*(_cptr->_bod));
            }
            bod_t *operator->() {
                return const_cast<bod_t*>(_cptr->_bod);
            }
            /* Whilst it is not my taste to define overloads for binary operators (such as == and !=) as member
             * functions (nor is it the convention), template argument deduction requires it, given that
             * gtd::bh_cube<>::iterator is a nested class. */
            bool operator==(const body_iterator &other) const noexcept {
                return this->_cptr == other._cptr;
            }
            bool operator!=(const body_iterator &other) const noexcept {
                return this->_cptr != other._cptr;
            }
            friend class nn_iterator;
            friend class cell_iterator;
        };
        class iterator : public body_iterator {
            cube_t *_root = nullptr;
            using body_iterator::_cptr;
            using body_iterator::_sptr;
            using body_iterator::index;
            // const cube_t *_cptr; // pointer to current leaf node/cube
            // const cube_t *const *_sptr{}; // pointer to pointer to current node within parent node's array
            // unsigned char index{}; // position of _cptr within parent node's array
            void find_next_leaf() {
                index = 0;
                _sptr = _cptr->_sub;
                while (index < 8) {
                    if (*_sptr != nullptr) {
                        _cptr = *_sptr;
                        if (_cptr->_bod != nullptr) { // means _cptr is pointing to leaf node
                            return;
                        }
                        _sptr = _cptr->_sub;
                        index = 0;
                        continue;
                    }
                    ++index;
                    ++_sptr;
                }
            }
        public:
            iterator() = default; // iterator is unusable until assigned to a cube
            iterator(cube_t *root, bool is_end = false) : _root{root} { // doesn't have to be root
                /* The constructor initially points _cptr to the leaf node (containing a non-NULL pointer to a body)
                 * with the lowest octant coordinates (according to how I have labelled them - see further up). */
                if (root == nullptr) [[unlikely]]
                    throw std::invalid_argument{"Error: nullptr passed as root cube.\n"};
                if (is_end) {
                    if (root->_bod != nullptr) {
                        _cptr = root->parent;
                        return;
                    }
                    _cptr = root;
                    return;
                }
                _cptr = root;
                this->find_next_leaf();
            }
            iterator(const iterator &other) = default;
            iterator &operator++() {
                if (_cptr == _root->parent) // means iterator is already pointing to end
                    return *this;
                if (_cptr == _root) { // the "end" for a non-leaf root is the root itself, else its parent
                    if (_root->_bod != nullptr)
                        _cptr = _root->parent;
                    return *this;
                }
                if (index == 7) {
                    loop_start:
                    while (true) {
                        if ((_cptr = _cptr->parent) == _root)
                            return *this;
                        index = 1;
                        _sptr = _cptr->parent->_sub;
                        while (*_sptr++ != _cptr) ++index;
                        while (index < 8 && *_sptr == nullptr) {
                            ++_sptr;
                            ++index;
                        }
                        if (index != 8) {
                            _cptr = *_sptr;
                            break;
                        }
                    }
                    // while (true) {
                    //     index = 1;
                    //     _sptr = _cptr->parent->_sub;
                    //     while (*_sptr++ != _cptr) ++index;
                    //     if (index != 8) {
                    //         _cptr = *_sptr;
                    //         break;
                    //     }
                    // }
                } else {
                    while (++index < 8 && *++_sptr == nullptr);
                    if (index == 8) {
                        goto loop_start;
                    }
                    _cptr = *_sptr;
                }
                if (_cptr->_bod != nullptr)
                    return *this;
                this->find_next_leaf();
                return *this;
            }
            iterator operator++(int) {
                iterator _copy = *this;
                this->operator++();
                return _copy;
            }
            explicit operator bool() const noexcept {
                return _cptr->_bod != nullptr; // covers all cases
            }
            iterator &operator=(const iterator &other) = default;
            iterator &operator=(cube_t *cube) {
                if (cube == nullptr) [[unlikely]]
                    throw std::invalid_argument{"Error: nullptr passed as cube.\n"};
                _cptr = _root = cube;
                // _root = cube;
                this->find_next_leaf();
                return *this;
            }
            friend class nn_iterator;
            friend class cell_iterator;
            template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prog, bool merge, int coll, uint64_t mF,
                    uint64_t fF, bool bin>
            friend class system;
        };
        class nn_iterator : public body_iterator {// nearest-neighbour iterator: iterates over the bodies nearest to one
            cube_t *_bc; // pointer to the bh_cube containing body around which bodies will be found
            using body_iterator::_cptr;
            using body_iterator::_sptr;
            using body_iterator::index;
            // const cube_t *_cptr{}; // pointer to current neighbouring body-containing cube
            // int64_t _h = 0; // height above the _bc cube in the tree (negative means below _bc cube - hence signed int)
            // uint64_t _mh = 0; // max. height reached above _bc
            // const cube_t *const *_sptr{}; // pointer to array of pointers that _cptr resides in (in parent's array)
            // unsigned char index = 0; // index of NEXT pointer to cube in parent array
            bool touches(const cube_t *_cube) const noexcept { // tests whether the _bc cube is touching the _cptr cube
                return !(_bc->_top.x < _cube->_btm.x || _bc->_btm.x > _cube->_top.x ||
                         _bc->_top.y < _cube->_btm.y || _bc->_btm.y > _cube->_top.y ||
                         _bc->_top.z < _cube->_btm.z || _bc->_btm.z > _cube->_top.z);
            }
            void find_next_neighbour() noexcept {
                _sptr = _cptr->parent->_sub + index;
                while (true) {
                    while (index < 8) {
                        if (*_sptr != nullptr && *_sptr != _bc && this->touches(*_sptr)) {
                            if ((*_sptr)->_bod != nullptr) {
                                _cptr = *_sptr;
                                ++index;
                                return;
                            }
                            _cptr = *_sptr;
                            _sptr = _cptr->_sub;
                            index = 0;
                            continue;
                        }
                        ++index;
                        ++_sptr;
                    }
                    if (_cptr->parent == nullptr) {
                        _cptr = _bc;
                        return;
                    }
                    _sptr = _cptr->parent->_sub;
                    index = 1;
                    while (*_sptr++ != _cptr) ++index;
                    _cptr = _cptr->parent;
                }
            }
        public:
            nn_iterator(const body_iterator &body_it, bool is_end = false) : body_iterator{body_it._cptr}, _bc{_cptr} {
                if (_bc == nullptr || _bc->_bod == nullptr)
                    throw std::invalid_argument{"Error: body iterator does not contain a pointer to a valid "
                                                "body-containing gtd::bh_cube<> object.\n"};
                if (_bc->parent == nullptr || is_end) {
                    // _cptr = _bc; // define _bc as end
                    return; // _cptr is pointing to end (no nearest neighbours to iterate over)
                }
                while (_cptr->parent->parent != nullptr) // find children of root
                    _cptr = _cptr->parent;
                this->find_next_neighbour();
            }
            // long long height() const noexcept {
            //     return _h;
            // }
            bool has_body() const noexcept {
                return _cptr != _bc;
            }
            explicit operator bool() const noexcept {
                return _cptr != _bc;
            }
            nn_iterator &operator++() {
                if (_cptr == _bc) [[unlikely]] // means _cptr is pointing to end
                    return *this;
                this->find_next_neighbour();
                return *this;
            }
            nn_iterator operator++(int) {
                nn_iterator cpy = *this;
                this->operator++();
                return cpy;
            }
            // bod_t &operator*() {
            //     return const_cast<bod_t&>(*(this->_cptr->_bod));
            // }
            // bod_t *operator->() {
            //     return const_cast<bod_t*>(this->_cptr->_bod);
            // }
            friend class cell_iterator;
            template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prog, bool merge, int coll, uint64_t mF,
                    uint64_t fF, bool bin>
            friend class system;
        };
        class cell_iterator { // does not inherit from body_iterator as does not necessarily iterate over bodies
            const cube_t *_bc = nullptr; // pointer to the leaf in which the body is contained
            const cube_t *_root = nullptr;
            const cube_t *_cptr = nullptr; // is const in cell_iterator since a cell_iterator cannot alter the cube
            const cube_t *const *_sptr{};
            unsigned char index = 0; // these variables have the same uses as in the body_iterator class
            vec _bpos{}; // position of the body on which force will be calculated
            long double _theta; // maximum opening angle from body to cell for cell not to be recursively subdivided
            void find_first_cell() noexcept {
                _sptr = _root->_sub;
                while (index < 8) {
                    if (*_sptr != nullptr && *_sptr != _bc) {
                        if ((*_sptr)->_bod != nullptr || this->is_source_cell(*_sptr)) {
                            _cptr = *_sptr++; // _sptr is left pointing to next cell to investigate
                            ++index;
                            return;
                        }
                        _sptr = (*_sptr)->_sub;
                        index = 0;
                        continue;
                    }
                    ++index;
                    ++_sptr;
                }
                _cptr = _bc; // would only occur in the case of single body tree
            }
            void find_next_cell() noexcept {
                while (true) {
                    while (index < 8) {
                        if (*_sptr != nullptr && *_sptr != _bc) {
                            if ((*_sptr)->_bod != nullptr || this->is_source_cell(*_sptr)) {
                                _cptr = *_sptr++; // _sptr is left pointing to next cell to investigate
                                ++index;
                                return;
                            }
                            // _cptr = *_sptr;
                            _sptr = (*_sptr)->_sub;
                            index = 0;
                            continue;
                        }
                        ++_sptr;
                        ++index;
                    }
                    if ((_cptr = _cptr->parent)->parent == nullptr)
                        break;
                    index = 1;
                    _sptr = _cptr->parent->_sub;
                    while (*_sptr++ != _cptr) ++index;
                }
                _cptr = _bc; // have reached end
            }
            bool is_source_cell(const cube_t *cell) { // banking on this being inlined by the compiler (if not, macro!)
                return (cell->_top.x - cell->_btm.x)/vec_ops::distance(this->_bpos, cell->_com) < _theta;
            }
            void error_logging() {
                if (_root == nullptr)
                    throw std::invalid_argument{"Error: root cube passed is nullptr.\n"};
                if (_root->parent != nullptr)
                    throw std::invalid_argument{"Error: root cube supplied is not root, as parent is not nullptr.\n"};
                if (_theta >= PI || _theta < 0)
                    throw std::invalid_argument{"Error: opening angle cannot be negative or equal to or larger than PI "
                                                "radians.\n"};
            }
        public:
            cell_iterator(long double opening_angle) : _theta{opening_angle} {
                if (_theta >= PI || _theta < 0)
                    throw std::invalid_argument{"Error: opening angle cannot be negative or equal to or larger than PI "
                                                "radians.\n"};
            }
            cell_iterator(cube_t *root, long double opening_angle) : _root{root}, _theta{opening_angle}{
                /* This constructor will construct a cell_iterator object that needs to be assigned to a body_iterator
                 * before use (otherwise -> undefined behaviour (most likely seg. fault due to deref. nullptr)). */
                this->error_logging();
            }
            cell_iterator(cube_t *root, const vec &particle_pos, long double opening_angle) :
            _root{root}, _cptr{root}, _bpos{particle_pos}, _theta{opening_angle} {
                this->error_logging();
                if ((_bc = root->find_leaf(particle_pos)) == nullptr)
                    throw invalid_body_placement{"Error: particle does not lie within root cell.\n"};
                this->find_first_cell();
            }
            cell_iterator(cube_t *root, const body_iterator &body_it, long double opening_angle) :
            _root{root}, _cptr{body_it._cptr}, _bc{body_it._cptr}, _theta{opening_angle} {
                this->error_logging();
                if (_cptr->_bod == nullptr)
                    throw std::invalid_argument{"Error: iterator passed does not point to a body-containing cell.\n"};
                if (!(_bpos = _cptr->_bod->curr_pos).between(_root->_btm, _root->_top))
                    throw invalid_body_placement{"Error: particle does not lie within root cell.\n"};
                this->find_first_cell();
            }
            cell_iterator(const body_iterator &body_it, long double opening_angle) :
            cell_iterator{body_it._cptr->root(), body_it, opening_angle} {}
            void set_root(cube_t *root) { // could have been an assignment operator overload but this is more explicit
                if (root == nullptr || root->parent != nullptr)
                    throw std::invalid_argument{"Error: root cannot be nullptr and parent must be nullptr.\n"};
                _root = root;
            }
            cell_iterator &operator++() {
                if (_cptr == _bc)
                    return *this;
                this->find_next_cell();
                return *this;
            }
            cell_iterator operator++(int) {
                cell_iterator _cpy = *this;
                this->operator++();
                return _cpy;
            }
            const cube_t &operator*() const noexcept {
                return *_cptr;
            }
            const cube_t *operator->() const noexcept {
                return _cptr;
            }
            explicit operator bool() const noexcept {
                return _cptr != _bc;
            }
            cell_iterator &operator=(const cell_iterator &other) = default;
            cell_iterator &operator=(const body_iterator &body_it) {
                if (body_it._cptr->_bod == nullptr)
                    throw std::invalid_argument{"Error: body iterator does not point to body.\n"};
                if (!body_it._cptr->_bod->curr_pos.between(_root->_btm, _root->_top))
                    throw invalid_body_placement{body_it->id};
                this->_cptr = this->_bc = body_it._cptr;
                this->_bpos = body_it._cptr->_bod->curr_pos;
                this->find_first_cell();
                return *this;
            }
            template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prog, bool merge, int coll, uint64_t mF,
                    uint64_t fF, bool bin>
            friend class system;
        };
        // bh_cube() = default;
        // bh_cube(cube_t *parent_cube, const vec &btm, const vec &top) : parent{parent_cube},
        // _btm{btm}, _top{top} {this->set_corners();}
        // bh_cube(cube_t *parent_cube, vec &&btm, vec &&top) : parent{parent_cube},
        // _btm{std::move(btm)}, _top{std::move(_top)} {this->set_corners();}
        bh_cube(cube_t *parent_cube, const vec &btm, const vec &top, const bod_t *_body) : parent{parent_cube},
        _btm{btm}, _top{top}, _bod{_body} {
            this->set_corners();
            this->check_first_body();
        }
        bh_cube(cube_t *parent_cube, vec &&btm, vec &&top, const bod_t *_body) : parent{parent_cube},
        _btm{std::move(btm)}, _top{std::move(top)}, _bod{_body} {
            this->set_corners();
            this->check_first_body();
        }
        const bod_t *get_body() const noexcept {
            return this->_bod;
        }
        cube_t *root() noexcept {
            cube_t *_cpy = this;
            while (_cpy->parent != nullptr) {
                _cpy = _cpy->parent;
            }
            return _cpy;
        }
        cube_t *find_leaf(const vec &pos) { // finds the leaf node which contains the given position (if present)
            cube_t **sptr = this->_sub;
            unsigned char i = 0;
            while (i < 8) {
                if (*sptr != nullptr) {
                    if (pos.between((*sptr)->_btm, (*sptr)->_top)) {
                        if ((*sptr)->_bod != nullptr) {
                            return *sptr;
                        }
                        sptr = (*sptr)->_sub;
                        i = 0;
                        continue;
                    }
                }
                ++i;
                ++sptr;
            }
            return nullptr; // case for pos not being within any leaf node
        }
        bool add_body(const bod_t *_body) {
            if (_body == nullptr) [[unlikely]]
                return false;
            if (!_body->curr_pos.between(this->_btm, this->_top)) [[unlikely]]
                throw invalid_body_placement{_body->id};
            // if (this->_bod != nullptr) { // case for first split of cube into sub-cubes
            //     // if (!this->dispatch(this->_bod)) {
            //     //     throw invalid_body_placement{this->_bod->id};
            //     // }
            //     this->dispatch(this->_bod);
            //     this->_bod = nullptr;
            // }
            // this->_mass += _body->mass_;
            // if (!this->dispatch(_body))
            //     throw invalid_body_placement{_body->id};
#ifdef RUNNING_COM
            M prev_mass = this->_mass;
#endif
            this->_mass += _body->mass_;
#ifndef RUNNING_COM
#ifdef RUNNING_MR_SUM
            this->mr_sum += _body->mass_*_body->curr_pos;
#endif
#else
            this->_com = (this->_com*prev_mass + _body->mass_*_body->curr_pos)/this->_mass;
#endif
#ifdef BH_REC_NUM
            ++this->_num;
#endif
            this->dispatch(_body);
            return true;
        }
        M cube_mass() const {
            return this->_mass;
        }
        vec cube_com() const noexcept {
            // this->calc_com();
            return this->_com;
        }
#ifdef BH_REC_NUM
        uint64_t num_bods() const noexcept {
            return this->_num;
        }
#endif
        ~bh_cube() {
            // std::cout << "Which dtor? this = " << this << std::endl;
            if (this->parent == nullptr) [[unlikely]] { // only root node should do the destructing
                cube_t *cptr = this;
                cube_t **sptr;
                unsigned char i = 0;
                // unsigned char count = 0;
                while (true) {
                    sptr = cptr->_sub;
                    loop_start:
                    for (; i < 8; ++i, ++sptr) {
                        if (*sptr != nullptr) {
                            if ((*sptr)->_bod != nullptr) { // case for leaf node
                                delete *sptr;
                                *sptr = nullptr;
                            }
                            else {
                                cptr = *sptr;
                                goto loop_end;
                            }
                        }
                    }
                    if (cptr->parent == nullptr) // got back up to the root node/cube
                        break;
                    if (i == 8) {
                        i = 0;
                        sptr = cptr->parent->_sub;
                        while (*sptr++ != cptr) ++i; // get *sptr to point to the next cube that needs processing
                        // i = cptr - *(cptr->parent->_sub);
                        cptr = cptr->parent;
                        delete cptr->_sub[i];
                        cptr->_sub[i++] = nullptr;
                        goto loop_start;
                    }
                    loop_end:
                    i = 0;
                }
            }
        }
        const cube_t &operator[](unsigned char index) const {
            if (index > 7)
                throw std::invalid_argument{"gtd::bh_cube<M, R, T, rF>::operator[] error: index cannot be greater than "
                                            "7.\n"};
            if (this->_sub[index] == nullptr)
                throw std::invalid_argument{"Error: no gtd::bh_cube<> at index specified.\n"};
            return *(this->_sub[index]);
        }
        // template <isNumWrapper m, isNumWrapper r, isNumWrapper t, uint64_t recFreq>
        friend std::ostream &operator<< <M, R, T, rF>(std::ostream &os, const bh_cube &cube);
        friend class bh_tree<M, R, T, rF>;
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, uint64_t, uint64_t, bool>
        friend class system; // friend class declaration cannot be a partial specialisation
    };
    template <isNumWrapper m, isNumWrapper r, isNumWrapper t, uint64_t recFreq>
    std::ostream &operator<<(std::ostream &os, const bh_cube<m, r, t, recFreq> &cube) {
        /* The ONLY function in which I have allowed myself to use recursion (since it is not an important function).
         * Beware of stack overflows (though unlikely to ever occur due to shallow recursion depth)! */
        const bh_cube<m, r, t, recFreq> *const *ptr = cube._sub;
        os << "[gtd::bh_cube@" << &cube <<
#ifdef BH_REC_NUM
        ":num_bods=" << cube._num << ",mass="
#else
        ":mass="
#endif
        << cube._mass << ",com=" << cube._com << ",body=";
        if (cube._bod == nullptr)
            os << "NULL";
        else
            os << *cube._bod;
        os << ",cubes:\n";
        for (unsigned char i = 0; i < 8; ++i, ++ptr) {
            os << "cube_" << +i << '=';
            if (*ptr == nullptr)
                os << "NULL";
            else
                os << '\n' << **ptr;
            os << '\n';
        }
        return os << ']';
    }
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t rF>
    std::ostream &operator<<(std::ostream&, const bh_tree<M, R, T, rF>&);
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t rF>
    class bh_tree {
        using cube_t = typename bh_cube<M, R, T, rF>::cube_t;
        using bod_t = typename bh_cube<M, R, T, rF>::bod_t;
        using vec = typename bh_cube<M, R, T, rF>::vec;
        cube_t *_root = nullptr;
        vec _btm{};
        vec _top{};
        void find_min_max(const bod_t *_ptr, typename std::vector<bod_t>::size_type _size) {
            T min_x = _ptr->curr_pos.x;
            T min_y = _ptr->curr_pos.y;
            T min_z = _ptr->curr_pos.z;
            T max_x = _ptr->curr_pos.x;
            T max_y = _ptr->curr_pos.y;
            T max_z = _ptr++->curr_pos.z;
            while (--_size > 0) {
                if (_ptr->curr_pos.x < min_x)
                    min_x = _ptr->curr_pos.x;
                if (_ptr->curr_pos.y < min_y)
                    min_y = _ptr->curr_pos.y;
                if (_ptr->curr_pos.z < min_z)
                    min_z = _ptr->curr_pos.z;
                if (_ptr->curr_pos.x > max_x)
                    max_x = _ptr->curr_pos.x;
                if (_ptr->curr_pos.y > max_y)
                    max_y = _ptr->curr_pos.y;
                if (_ptr->curr_pos.z > max_z)
                    max_z = _ptr->curr_pos.z;
                ++_ptr;
            }
            T x_range = max_x - min_x;
            T y_range = max_y - min_y;
            T z_range = max_z - min_z;
            T longest_edge = x_range > y_range ? (x_range > z_range ? x_range : z_range) :
                                                 (y_range > z_range ? y_range : z_range);
            T pad = longest_edge/20.0l; // to leave a little wiggle room
            _btm.x = min_x - pad;
            _btm.y = min_y - pad;
            _btm.z = min_z - pad;
            longest_edge += 2*pad;
            _top.x = _btm.x + longest_edge;
            _top.y = _btm.y + longest_edge;
            _top.z = _btm.z + longest_edge;
        }
    public:
        bh_tree() = default;
        vec bottom() {
            return this->_btm;
        }
        vec top() {
            return this->_top;
        }
#ifdef BH_REC_NUM
        uint64_t num_bods() const noexcept {
            return _root == nullptr ? 0 : _root->num_bods();
        }
#endif
#ifndef RUNNING_COM
        void calc_com() {
            _root->calc_com();
        }
#endif
        void add_bodies(const std::vector<bod_t> &bods) {
            if (bods.empty())
                return;
            const bod_t *ptr = bods.data();
            typename std::vector<bod_t>::size_type _size = bods.size();
            if (_root == nullptr) {
                this->find_min_max(ptr, _size--);
                _root = new cube_t{nullptr, _btm, _top, ptr++};
            }
            while (_size --> 0) {
                _root->add_body(ptr++);
            }
        }
        void assign_bodies(const std::vector<bod_t> &bods) {
            delete _root;
            _root = nullptr;
            this->add_bodies(bods);
        }
        typename cube_t::iterator begin() {
            return {_root};
        }
        typename cube_t::iterator end() {
            return {_root, true};
        }
        typename cube_t::cell_iterator cell_begin(const typename cube_t::body_iterator &body_it,
                                                  long double opening_angle) {
            return {_root, body_it, opening_angle};
        }
        typename cube_t::nn_iterator nn_begin(const typename cube_t::body_iterator &body_it) {
            return {body_it};
        }
        void delete_tree() noexcept {
            delete _root;
            _root = nullptr;
            _btm.x = T{};
            _btm.y = T{};
            _btm.z = T{};
            _top.x = T{};
            _top.y = T{};
            _top.z = T{};
        };
        ~bh_tree() {
            delete _root; // safe to delete nullptr
        }
        friend std::ostream &operator<< <M, R, T, rF>(std::ostream&, const bh_tree<M, R, T, rF>&);
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, uint64_t, uint64_t, bool>
        friend class system;
    };
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t rF>
    std::ostream &operator<<(std::ostream &os, const bh_tree<M, R, T, rF> &tree) {
        if (tree._root == nullptr)
            return os << "NULL";
        return os << *tree._root;
    }
    template <uint64_t recFreq>
    using cube = bh_cube<long double, long double, long double, recFreq>;
    typedef bh_cube<long double, long double, long double, 0> cube_0f;
    typedef bh_cube<long double, long double, long double, 1> cube_1f;
    typedef bh_cube<long double, long double, long double, 10> cube_10f;
    typedef bh_cube<long double, long double, long double, 100> cube_100f;
    typedef bh_cube<long double, long double, long double, 1000> cube_1000f;
    template <uint64_t recFreq>
    using octree = bh_tree<long double, long double, long double, recFreq>;
    typedef bh_tree<long double, long double, long double, 0> octree_0f;
    typedef bh_tree<long double, long double, long double, 1> octree_1f;
    typedef bh_tree<long double, long double, long double, 10> octree_10f;
    typedef bh_tree<long double, long double, long double, 100> octree_100f;
    typedef bh_tree<long double, long double, long double, 1000> octree_1000f;
}
#endif
