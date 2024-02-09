/* Copyright (c) 2023 Gregor Anton Randall Hartl Watters
 * This software is protected under the MIT license. Please see the LICENSE file for more information. */

#ifndef GREGSYS_HPP
#define GREGSYS_HPP

#ifndef __cplusplus
#error "The gregsys.hpp header file is a C++ header file only.\n"
#endif

#include "greg8tree.hpp"
#include <random>

#define EMPTY

#define FUNC_TEMPL_SELECT(func, bfunc, cf, ...) \
if constexpr (collisions == overlap_coll_check || !(memFreq && fileFreq)) { \
    func<false, false>(__VA_ARGS__);            \
    bfunc \
    if constexpr (collisions == overlap_coll_check) {                       \
        if constexpr (memFreq && fileFreq) {                           \
            if (!(steps % memFreq) && !(steps % fileFreq)) { \
                this->s_coll<true, true>(); \
            } \
            else if (!(steps % memFreq)) { \
                this->s_coll<true, false>(); \
            } \
            else if (!(steps % fileFreq)) { \
                this->s_coll<false, true>(); \
                cf \
            } \
            else { \
                this->s_coll<false, false>(); \
                cf \
            } \
        }                                \
        else if constexpr (memFreq) { \
            if (steps % memFreq) { \
                this->s_coll<false, false>(); \
                cf \
            } \
            else { \
                this->s_coll<true, false>(); \
            } \
        } \
        else if constexpr (fileFreq) { \
            if (steps % fileFreq) { \
                this->s_coll<false, false>(); \
                cf \
            } \
            this->s_coll<false, true>(); \
            cf \
        }                                \
        else {                            \
            this->s_coll<false, false>();\
            cf \
        }\
    }                                    \
    else { \
        cf                               \
    }\
} \
else { \
    if constexpr (memFreq && fileFreq) { \
        if (!(steps % memFreq) && !(steps % fileFreq)) { \
            func<true, true>(__VA_ARGS__);      \
            bfunc \
        } \
        else if (!(steps % memFreq)) { \
            func<true, false>(__VA_ARGS__);     \
            bfunc \
        } \
        else if (!(steps % fileFreq)) { \
            func<false, true>(__VA_ARGS__);     \
            bfunc \
            cf \
        } \
        else { \
            func<false, false>(__VA_ARGS__);    \
            bfunc \
            cf \
        } \
    } \
    else if constexpr (memFreq) { \
        if (steps % memFreq) { \
            func<false, false>(__VA_ARGS__);    \
            bfunc \
            cf \
        } \
        else { \
            func<true, false>(__VA_ARGS__);     \
            bfunc \
        } \
    } \
    else if constexpr (fileFreq) { \
        if (steps % fileFreq) { \
            func<false, false>(__VA_ARGS__);    \
            bfunc \
            cf \
        } \
        func<false, true>(__VA_ARGS__);         \
        bfunc \
        cf \
    } \
}

#define MEM_LOOP \
if constexpr (memFreq) { \
    if (steps % memFreq) \
        for (bod_t &bod : bods) \
            bod.acc.make_zero(); \
    else {       \
        total_pe = total_ke =  T{0}; \
        for (bod_t &bod : bods) { \
            bod.acc.make_zero(); \
            bod.pe = T{0}; \
        } \
    }\
} \
else \
    for (bod_t &bod : bods) \
        bod.acc.make_zero();

namespace gtd {
    class nsys_load_error : public nbody_error {
    public:
        nsys_load_error() : nbody_error{"Error loading N-body data from .nsys file.\n"} {}
        explicit nsys_load_error(const char *msg) : nbody_error{msg} {}
    };
    class nsys_format_error : public nsys_load_error {
    public:
        nsys_format_error() : nsys_load_error{"Error loading N-body data from .nsys file.\n"} {}
        explicit nsys_format_error(const char *msg) : nsys_load_error{msg} {}
    };
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double,
              bool prog = false, bool mergeOverlappingBodies = false, int collisions = 0,
              uint64_t memFreq = 0, uint64_t fileFreq = 1, bool binaryFile = true>
    class system {
        /* The system class represents a system in 3D space that is composed of bodies which interact with each other
         * gravitationally (and through contact force in case of collisions - if this option is specified using the
         * collisions template parameter (3 or 7)). The bodies in the system are represented by gtd::body objects. */
        /* A system object is capable of evolving, in time, the positions, velocities & energies of the bodies via
         * various different integration methods, each specified by an integral constant (just below). */
        /* The 'M', 'R' & 'T' template parameters are the data types used to represent the bodies' masses, radii, &
         * positions and velocities, respectively. The other (non-type) template parameters are present in order to
         * avoid a small runtime overhead in checking which have been set and/or to avoid duplication of code. The
         * 'prog' parameter determines whether progress is printed to standard output during the main evolution of the
         * bodies. 'mergeOverlappingBodies' determines whether bodies that are added to a system object are merged
         * together if they overlap, or whether an exception is thrown instead. 'collisions' determines how collisions
         * should be predicted and how they should be dealt with (0 = not at all, 3 = simple - based on overlap, 7 =
         * advanced (predictive) collision checking). 'memFreq' represents the frequency with which position, velocity
         * and energy data values should be stored in memory (on the heap) and, similarly, fileFreq represents the
         * frequency with which they should be written to a file. Finally, 'binaryFile' represents whether the data
         * being written to a file is in binary format or text (.csv) format (if fileFreq is non-zero). */
    public:
        static constexpr int two_body = 1; // static constants to determine which integration method to perform
        static constexpr int euler = 2;
        static constexpr int modified_euler = 4;
        static constexpr int midpoint = 8;
        static constexpr int leapfrog_kdk = 16;
        static constexpr int leapfrog_dkd = 32;
        static constexpr int rk4 = 64;
        static constexpr int rk3_8 = 128;
        static constexpr int barnes_hut = 65'536; // static constants to determine the force approx. method to use
        static constexpr int fast_multipole = 131'072;
        static constexpr int no_coll_check = 0;
        static constexpr int overlap_coll_check = 3;
        static constexpr int pred_coll_check = 7;
        static constexpr long double G_SI = 0.000'000'000'066'743; // G in SI units (m^3 kg^-1 s^-2)
    private: // G, below, has not been made a static constant as it varies between objects, depending on units passed
        long double G = 66'743; // Newtonian constant of Gravitation (10^(-15) m^3 kg^-1 s^-2)
        using vec_t = vector3D<T>;
        using bod_t = body<M, R, T, memFreq>;
        using vec_size_t = typename std::vector<bod_t>::size_type;
        using sys_t = system<M, R, T, prog, mergeOverlappingBodies, collisions, memFreq, fileFreq, binaryFile>;
        using mom_t = decltype(M{}*T{});
        using ke_t = decltype(std::declval<M>()*std::declval<T>()*std::declval<T>());
        using pe_t = decltype((std::declval<long double>()*std::declval<M>()*std::declval<M>())/std::declval<T>());
        using cube_t = bh_cube<M, R, T, memFreq>;
        using tree_t = bh_tree<M, R, T, memFreq>;
        std::vector<bod_t> bods; // not using set, map or list as need fast random access to bodies
        tree_t btree{}; // Barnes-Hut body tree - only used in evolve() if Barnes-Hut or fast multiple option selected
        // std::set<bod_t, std::less<>> del_bods; // set to store bodies removed from system after collision mergers
        /* Here I declare a std::set to store bodies that are removed from the evolution in the case of a merger. I have
         * declared the type of the comparator used by the std::set object to be that of an anonymous lambda type, which
         * will cause the std::pair objects to be sorted first by the iteration at which they were "destroyed", and then
         * by their IDs (for when times of destruction are equal). */
        // std::set<std::pair<ull_t, bod_t>,
        //         decltype([](const std::pair<ull_t, bod_t> &p1, const std::pair<ull_t, bod_t> &p2){
        //             return p1.first == p2.first ? p1.second < p2.second : p1.first < p2.first;
        //         })> del_bods;
#ifdef GREGSYS_MERGERS
        static inline auto pair_bods_func = [](const std::pair<uint64_t, bod_t> &p1,
                                               const std::pair<uint64_t, bod_t> &p2) {
            return p1.first == p2.first ? p1.second < p2.second : p1.first < p2.first;
        };
        union {
            std::set<bod_t, std::less<>> *del_bods_n;
            std::set<std::pair<uint64_t, bod_t>, decltype(pair_bods_func)> *del_bods_m;
        };
#endif
#ifdef GREGSYS_SOFTENING
        long double eps = 0; // softening
        long double eps_sq = eps*eps;
#endif
        long double dt = 1;
        long double half_dt = dt/2; // I gave it its own variable, since it is a commonly used quantity
        uint64_t iterations = 1'000;
        long double prev_dt{}; // used in methods called after evolve(), since could be changed by setter
        long double time_elapsed{};
        uint64_t prev_iterations{}; // same here
#ifdef GREGSYS_MERGERS
        mom_t min_tot_com_mom{BILLION*10}; // minimum sum of magnitudes of COM momenta for two bodies to merge
#endif
        // using P = decltype((G*std::declval<M>()*std::declval<M>())/std::declval<T>()); // will be long double
        long double bh_theta = PI/3; // opening angle parameter for Barnes-Hut force approximation
        uint64_t steps{}; // defined as an instance variable since it's required in numerous functions
        T total_pe{}; // these 5 variables are defined as instance variables to avoid their redefinition in many funcs
        T total_ke{};
        uint64_t inner{};
        uint64_t outer{};
        mutable vec_size_t num_bods{};
        std::map<uint64_t, std::vector<T>> pe; // to store potential energies of bodies
        std::map<uint64_t, std::vector<T>> energy; // total energy for each body at each iteration (KE + PE)
        //std::map<unsigned long long, std::vector<T>> del_ke; // map to store KE vectors of deleted bodies from mergers
        /* it would have been nice to use std::maps to store potential and total energies for all bodies, as it would
         * have allowed lookup by body_id, but its lookup complexity is O(log(N)), whereas accessing an element in a
         * std::vector by index is O(1) */
        std::vector<T> tot_pe; // total potential energy for the entire system at each iteration
        std::vector<T> tot_ke; // total kinetic energy for the entire system at each iteration
        std::vector<T> tot_e; // total energy for the entire system at each iteration (KE + PE)
        bool evolved = false; // indicates whether a gtd::system object is in an evolved state
        uint64_t *G_vals = new uint64_t[6];
        union {
            mutable std::ofstream *ostream = nullptr; // to write the system's data to a file
            mutable std::ifstream *istream; // to read data from .nsys files
        };
        static inline String def_path = get_home_path<String>() + FILE_SEP + "System_Trajectories_";
        static inline uint64_t to_ull(String &&str) {
            if (!str.isnumeric())
                return -1;
            uint64_t total = 0;
            for (const char &c : str) {
                total *= 10;
                total += c - 48;
            }
            return total;
        }
        void parse_units_format(String &&str) {
            if (!std::regex_match(str.c_str(), std::regex(R"(^\s*(m|M)\s*(?=\d{1,18}\s*:)\d*[1-9]+\d*\s*:\s*(?=\d{1,18}\s*,)\d*[1-9]+\d*\s*,\s*(d|D)\s*(?=\d{1,18}\s*:)\d*[1-9]+\d*\s*:\s*(?=\d{1,18}\s*,)\d*[1-9]+\d*\s*,\s*(t|T)\s*(?=\d{1,18}\s*:)\d*[1-9]+\d*\s*:\s*(?=\d{1,18}\s*$)\d*[1-9]+\d*\s*$)"))) {
                str.append_front("The units_format string passed, \"").append_back("\", does not match the format "
                                                                                   "required:\n\"Ma:b,Dc:x,Ty:z\", "
                                                                                   "where 'a', 'b', 'c', 'x', 'y', and "
                                                                                   "'z' represent any positive integer "
                                                                                   "number of 1-18 characters starting "
                                                                                   "with at least one non-zero "
                                                                                   "character.\n");
                throw std::invalid_argument(str.c_str());
            }
            str.strip("MmDdTt ");
            size_t colon_index = str.find(':');
            size_t comma_index = str.find(',');
            uint64_t m_denom = to_ull(str.substr(0, colon_index)); // denominator
            uint64_t m_num = to_ull(str.substr(colon_index + 1, comma_index)); // numerator
            str.erase_chars(0, comma_index + 1);
            colon_index = str.find(':');
            comma_index = str.find(',');
            uint64_t d_num = to_ull(str.substr(0, colon_index));
            uint64_t d_denom = to_ull(str.substr(colon_index + 1, comma_index));
            str.erase_chars(0, comma_index + 1);
            colon_index = str.find(':');
            comma_index = str.find(',');
            uint64_t t_denom = to_ull(str.substr(0, colon_index));
            uint64_t t_num = to_ull(str.substr(colon_index + 1, comma_index));
            G = 66'743; // 10^(-15) m^3 kg^-1 s^-2
            G *= m_num*d_num*d_num*d_num*t_num*t_num;
            G /= m_denom*d_denom*d_denom*d_denom*t_denom*t_denom;
            G /= 1'000'000'000'000'000; // correcting for the original G being 10^(15) times larger than it should be
            *G_vals++ = m_num;
            *G_vals++ = m_denom;
            *G_vals++ = d_num;
            *G_vals++ = d_denom;
            *G_vals++ = t_num;
            *G_vals = t_denom;
            G_vals -= 5;
        }
#ifdef GREGSYS_MERGERS
        void set_set() requires (collisions != 0) {
            if constexpr (memFreq)
                del_bods_m = new std::set<std::pair<ull_t, bod_t>, decltype(pair_bods_func)>{pair_bods_func};
            else
                del_bods_n = new std::set<bod_t, std::less<>>;
        }
#endif
        void clear_bodies() requires (memFreq != 0) { // makes sure the trajectories of all bodies are deleted
            // for (bod_t &bod : bods)
            //     bod.clear();
            std::for_each(bods.begin(), bods.end(), [this](bod_t &bod){bod.clear();});
        }
        void clear_bodies(vec_size_t &as_of) requires (memFreq != 0) {
            // if (as_of >= bods.size())
            //     return;
            // vec_size_t size = bods.size();
            // while (as_of < size)
            //     bods[as_of++].clear();
            std::for_each(bods.begin() + as_of, bods.end(), [this](bod_t &bod){bod.clear();});
        }
        void clear_bodies(vec_size_t &&as_of) requires (memFreq != 0) {
            clear_bodies(as_of);
        }
        void check_overlap() {
            bod_t *outer_b;
            bod_t *inner_b;
            num_bods = bods.size();
            for (outer = 0; outer < num_bods; ++outer) {
                outer_b = bods.data() + outer;
                inner_b = outer_b + 1;
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    if constexpr (mergeOverlappingBodies)
                        if (inner == outer) { // will be true at some point if there has been a merger
                            ++inner_b;
                            continue;
                        }
                    // if ((outer_b->curr_pos - inner_b->curr_pos).magnitude() < outer_b->radius + inner_b->radius)
                    if (vec_ops::distance(outer_b->curr_pos, inner_b->curr_pos) < outer_b->radius + inner_b->radius) {
                        if constexpr (!mergeOverlappingBodies) {
                            String str = "The bodies with id=";
                            str.append_back(outer_b->id).append_back(" and id=").append_back(inner_b->id);
                            str.append_back(" that were added to this system object overlap.\n");
                            throw overlapping_bodies_error(str.c_str());
                        }
                        *outer_b += *inner_b; // merges the two overlapping bodies
                        bods.erase(bods.begin() + inner); // thus the number of total bodies is reduced by 1
                        if (inner < outer) {
                            --outer;
                            --outer_b;
                        }
                        inner = 0; // have to recalculate possible mergers for newly created body
                        inner_b = bods.data();
                        --num_bods;
                        continue;
                    }
                    ++inner_b;
                }
            }
        }
        void cumulative_acc_and_pe(bod_t &b1, bod_t &b2) requires (memFreq != 0) {
            vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
#ifdef GREGSYS_SOFTENING
            // vector3D<T> &&r12_uvec = r12.unit_vector();
            auto &&r12_sq = r12*r12;
            auto &&r12_sq_eps = r12_sq + this->eps_sq; // square of Euclidean distance in "4D" space
            r12 /= sqrtl(r12_sq);
            b1.acc += ((G*b2.mass_)/(r12_sq_eps))*r12;//_uvec;
            b2.acc -= ((G*b1.mass_)/(r12_sq_eps))*r12;//_uvec;
            auto &&pot_energy = -(G*b1.mass_*b2.mass_)/sqrtl(r12_sq_eps);
#else
            long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
            b1.acc += ((G*b2.mass_)/(r12_cubed_mag))*r12;
            b2.acc -= ((G*b1.mass_)/(r12_cubed_mag))*r12;
            auto &&pot_energy = -(G*b1.mass_*b2.mass_)/r12.magnitude();
#endif
            b1.pe += pot_energy;
            b2.pe += pot_energy;
        }
        void cumulative_acc(bod_t &b1, bod_t &b2) {
            vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
#ifdef GREGSYS_SOFTENING
            // vector3D<T> &&r12_uvec = r12.unit_vector();
            auto &&r12_sq = r12*r12;
            r12 /= r12.magnitude();
            auto &&r12_sq_eps = r12_sq + this->eps_sq; // square of Euclidean distance in "4D" space
            b1.acc += ((G*b2.mass_)/(r12_sq_eps))*r12;//_uvec;
            b2.acc -= ((G*b1.mass_)/(r12_sq_eps))*r12;//_uvec;
#else
            long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
            b1.acc += ((G*b2.mass_)/(r12_cubed_mag))*r12;
            b2.acc -= ((G*b1.mass_)/(r12_cubed_mag))*r12;
#endif
        }
        void cumulative_pe(bod_t &b1, bod_t &b2) requires (memFreq != 0) {
            vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
#ifdef GREGSYS_SOFTENING
            auto &&pot_energy = -(G*b1.mass_*b2.mass_)/sqrtl(r12*r12 + this->eps_sq);
#else
            auto &&pot_energy = -(G*b1.mass_*b2.mass_)/r12.magnitude();
#endif
            b1.pe += pot_energy;
            b2.pe += pot_energy;
        }
        void bh_acc_and_pe(bod_t *_bod, const cube_t *_src) requires (memFreq != 0) {
            vec_t &&r12 = _src->_com - _bod->curr_pos;
#ifdef GREGSYS_SOFTENING
            auto &&r12_sq = r12*r12 + eps_sq;
#else
            auto &&r12_sq = r12*r12;
#endif
            auto &&r12_mag = sqrtl(r12_sq);
            _bod->acc = this->G*(_src->_mass/(r12_sq))*(r12/r12_mag);
            _bod->pe -= -(G*_bod->mass_*_src->_mass)/r12_mag;
        }
        void bh_acc(bod_t *_bod, const cube_t *_src) {
            vec_t &&r12 = _src->_com - _bod->curr_pos;
#ifdef GREGSYS_SOFTENING
            auto &&r12_sq = r12*r12 + eps_sq;
#else
            auto &&r12_sq = r12*r12;
#endif
            _bod->acc = this->G*(_src->_mass/(r12_sq))*(r12/(sqrtl(r12_sq)));
        }
        void bh_pe(bod_t *_bod, const cube_t *_src) requires (memFreq != 0) {
            vector3D<T> &&r12 = _src->_com - _bod->curr_pos;
#ifdef GREGSYS_SOFTENING
            _bod->pe -= (G*_bod->mass_*_src->_mass)/(sqrtl(r12*r12 + eps_sq));
#else
            _bod->pe -= (G*_bod->mass_*_src->_mass)/r12.magnitude();//sqrtl(r12*r12);
#endif
        }
        template <bool mem, bool file>
        void take_euler_step() {
            for (bod_t &bod : bods) {
                bod.curr_pos += bod.curr_vel*dt;
                bod.curr_vel += bod.acc*dt;
                if constexpr (mem)
                    bod.add_pos_vel_ke();
                if constexpr (file) {
                    // WRITE TO FILE
                }
            }
        }
        template <bool mem, bool file>
        void take_modified_euler_step(const std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>>
                                      &predicted_vals) {
            outer = 0;
            for (bod_t &bod : bods) {
                bod.curr_pos += half_dt*(bod.curr_vel + std::get<1>(predicted_vals[outer]));
                bod.curr_vel += half_dt*(bod.acc + std::get<2>(predicted_vals[outer++]));
                if constexpr (mem)
                    bod.add_pos_vel_ke();
                if constexpr (file) {
                    // WRITE TO FILE
                }
            }
        }
        template <bool mem, bool file>
        void take_midpoint_step(const std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>>
                                &predicted_vals) {
            outer = 0;
            for (bod_t &bod : bods) {
                bod.curr_pos += dt*std::get<1>(predicted_vals[outer]);
                bod.curr_vel += dt*std::get<2>(predicted_vals[outer++]);
                if constexpr (mem)
                    bod.add_pos_vel_ke();
                if constexpr (file) {
                    // WRITE TO FILE
                }
            }
        }
        void create_energy_vectors() requires (memFreq != 0) {
            uint64_t iters_p1 = (iterations + 1)/memFreq;
            // pe.resize(bods_size);
            // energy.resize(bods_size);
            pe.clear();
            energy.clear();
            for (bod_t &bod : bods) {
                // pe[i].reserve(iters_p1);
                // energy[i].reserve(iters_p1);
                pe[bod.id].reserve(iters_p1);
                energy[bod.id].reserve(iters_p1);
                bod.positions.reserve(iters_p1);
                bod.velocities.reserve(iters_p1);
                bod.energies.reserve(iters_p1);
            }
            tot_pe.reserve(iters_p1);
            tot_ke.reserve(iters_p1);
            tot_e.reserve(iters_p1);
        }
        void calc_acc_and_e() {
            MEM_LOOP // zeros out total_pe, total_ke and pe of each body if memFreq, and always zeros all accelerations
            num_bods = bods.size();
            if constexpr (!memFreq) {
                for (outer = 0; outer < num_bods; ++outer) {
                    bod_t &ref = bods[outer];
                    for (inner = outer + 1; inner < num_bods; ++inner)
                        this->cumulative_acc(ref, bods[inner]);
                }
            }
            else {
                if (steps % memFreq) {
                    for (outer = 0; outer < num_bods; ++outer) {
                        bod_t &ref = bods[outer];
                        for (inner = outer + 1; inner < num_bods; ++inner)
                            this->cumulative_acc(ref, bods[inner]);
                    }
                    return;
                }
                else {
                    for (outer = 0; outer < num_bods; ++outer) {
                        bod_t &ref = bods[outer];
                        for (inner = outer + 1; inner < num_bods; ++inner)
                            this->cumulative_acc_and_pe(ref, bods[inner]);
                        ref.pe /= T{2};
                        pe[ref.id].push_back(ref.pe);
                        energy[ref.id].push_back(ref.curr_ke + ref.pe);
                        total_pe += ref.pe;
                        total_ke += ref.curr_ke;
                    }
                }
            }
            if constexpr (memFreq) {
                tot_pe.push_back(total_pe);
                tot_ke.push_back(total_ke);
                tot_e.push_back(total_pe + total_ke);
            }
        }
        template <bool mem, bool file, bool only_e = false>
        void loop_kdk() {
            num_bods = bods.size();
            for (outer = 0; outer < num_bods; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    if constexpr (mem && memFreq)
                        this->cumulative_acc_and_pe(ref, bods[inner]);
                    else if constexpr (only_e && memFreq)
                        this->cumulative_pe(ref, bods[inner]);
                    else
                        this->cumulative_acc(ref, bods[inner]);
                }
                if constexpr (only_e && memFreq) {
                    ref.pe /= T{2};
                    pe[ref.id].push_back(ref.pe);
                    energy[ref.id].push_back(ref.curr_ke + ref.pe);
                    total_pe += ref.pe;
                    total_ke += ref.curr_ke;
                }
                else {
                    /* KICK for half a step */
                    ref.curr_vel += half_dt*ref.acc;
                    if constexpr (mem && memFreq) {
                        ref.pe /= T{2};
                        ref.add_pos_vel_ke(); // store new particle position, velocity and kinetic energy
                        pe[ref.id].push_back(ref.pe);
                        energy[ref.id].push_back(ref.curr_ke + ref.pe);
                        total_pe += ref.pe;
                        total_ke += ref.curr_ke;
                    }
                    if constexpr (file && fileFreq) {
                        // WRITE TO FILE
                    }
                }
            }
        }
        void leapfrog_kdk_acc_e_and_step() {
            MEM_LOOP
            FUNC_TEMPL_SELECT(loop_kdk, EMPTY, return;)
            /* repeatedly evaluating the same if constexpr conditions DOES NOT MATTER (apart from slightly increasing
             * compilation time) given that they are evaluated at compile time, so there is never any runtime overhead*/
            if constexpr (collisions == overlap_coll_check) {
                loop_kdk<false, false, true>(); // only energy calculations
            }
            if constexpr (memFreq) { // in the FUNC_TEMPL_SELECT macro it is checked whether not steps % memFreq
                tot_pe.push_back(total_pe);
                tot_ke.push_back(total_ke);
                tot_e.push_back(total_pe + total_ke);
            }
        }
        template <typename bodsFuncT> // no need to re-check requirements here, already checked in evolve()
        void leapfrog_kdk_acc_e_and_step(const bodsFuncT &func) {
            MEM_LOOP
            FUNC_TEMPL_SELECT(loop_kdk, func(this->bods);, return;)
            if constexpr (collisions == overlap_coll_check)
                loop_kdk<false, false, true>();
            if constexpr (memFreq) {
                tot_pe.push_back(total_pe);
                tot_ke.push_back(total_ke);
                tot_e.push_back(total_pe + total_ke);
            }
        }
        template <bool mem, bool file>
        void loop_dkd() {
            for (outer = 0; outer < num_bods; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner)
                    cumulative_acc(ref, bods[inner]);
                /* KICK for a full step */
                ref.curr_vel += dt*ref.acc;
                /* DRIFT for half a step */
                ref.curr_pos += half_dt*ref.curr_vel;
                /* the reason it is possible to update the outer body's position within the outer loop (without
                 * affecting the synchronisation of the particles) is because, by here, its effect (at its
                 * now-previous position) on all the other particles in the system has been calculated (all
                 * subsequent updates of the accelerations of the other particles no longer depend on the position
                 * of the outer body) */
                if constexpr (mem && memFreq) {
                    ref.add_pos_vel_ke(); // store new particle position, velocity and kinetic energy
                    total_ke += ref.curr_ke;
                }
                if (file && fileFreq) {
                    // WRITE TO FILE
                }
            }
        }
        void leapfrog_dkd_acc_e_and_step() {
            MEM_LOOP
            num_bods = bods.size();
            FUNC_TEMPL_SELECT(loop_dkd, EMPTY, return;)
            /* a second loop is required to compute the potential energies based on the updated positions */
            if constexpr (memFreq) { // again, only reached if not steps % memFreq (checked within FUNC_TEMPL_SELECT)
                for (outer = 0; outer < num_bods; ++outer) {
                    bod_t &ref = bods[outer];
                    for (inner = outer + 1; inner < num_bods; ++inner)
                        cumulative_pe(ref, bods[inner]);
                    ref.pe /= T{2};
                    pe[ref.id].push_back(ref.pe);
                    energy[ref.id].push_back(ref.curr_ke + ref.pe);
                    total_pe += ref.pe;
                    if constexpr (collisions == overlap_coll_check)
                        total_ke += ref.curr_ke;
                }
                tot_pe.push_back(total_pe);
                tot_ke.push_back(total_ke);
                tot_e.push_back(total_pe + total_ke);
            }
        }
        template <typename bodsFuncT>
        void leapfrog_dkd_acc_e_and_step(const bodsFuncT &func) {
            MEM_LOOP
            num_bods = bods.size();
            FUNC_TEMPL_SELECT(loop_dkd, func(this->bods);, return;)
            if constexpr (memFreq) {
                for (outer = 0; outer < num_bods; ++outer) {
                    bod_t &ref = bods[outer];
                    for (inner = outer + 1; inner < num_bods; ++inner)
                        cumulative_pe(ref, bods[inner]);
                    ref.pe /= T{2};
                    pe[ref.id].push_back(ref.pe);
                    energy[ref.id].push_back(ref.curr_ke + ref.pe);
                    total_pe += ref.pe;
                    if constexpr (collisions == overlap_coll_check)
                        total_ke += ref.curr_ke;
                }
                tot_pe.push_back(total_pe);
                tot_ke.push_back(total_ke);
                tot_e.push_back(total_pe + total_ke);
            }
        }
        void calc_energy() requires (memFreq != 0) {
            num_bods = bods.size();
            for (bod_t &bod : bods)
                bod.pe = T{};
            total_pe = total_ke = T{};
            for (outer = 0; outer < num_bods; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    // auto &&pot_energy = (G*ref.mass_*bods[inner].mass_)/(ref.curr_pos-bods[inner].curr_pos).magnitude();
                    // ref.pe -= pot_energy;
                    // bods[inner].pe -= pot_energy;
                    cumulative_pe(ref, bods[inner]);
                }
                ref.pe /= T{2};
                pe[ref.id].push_back(ref.pe);
                energy[ref.id].push_back(ref.curr_ke + ref.pe);
                total_pe += ref.pe;
                total_ke += ref.curr_ke;
            }
            tot_pe.push_back(total_pe);
            tot_ke.push_back(total_ke);
            tot_e.push_back(total_pe + total_ke);
        }
        void bh_s_coll() { // removed static from map
            std::map<long double, std::tuple<bod_t*, decltype(M{}*T{}), decltype(M{}*T{}), vector3D<T>>> overlapping;
            bod_t *bod_o;
            bod_t *bod_i;
            R rad_dist{};
            long double dist;
            mom_t axis_mom_o;
            mom_t axis_mom_i;
            auto it = btree.begin();
            while (it) {
                overlapping.clear();
                typename cube_t::nn_iterator nit{it};
                bod_o = const_cast<bod_t*>(it._cptr->_bod);
                while (nit) {
                    bod_i = const_cast<bod_t*>(nit++._cptr->_bod);
                    rad_dist = bod_o->radius + bod_i->radius;
                    vector3D<T> &&r12 = bod_i->curr_pos - bod_o->curr_pos;
                    dist = r12.magnitude(); // DEAL WITH ZERO DISTANCE CASE
                    r12.x /= dist; r12.y /= dist; r12.z /= dist; // more efficient than calling normalise()
                    if (rad_dist > dist) {
                        axis_mom_o = bod_o->momentum()*r12;
                        axis_mom_i = bod_i->momentum()*r12;
                        overlapping.emplace(std::piecewise_construct, std::forward_as_tuple(dist),
                                            std::forward_as_tuple(bod_i, axis_mom_o, axis_mom_i, r12));
                    }
                }
                if (!overlapping.empty()) {
                    T o_vel;
                    T i_vel;
                    T o_minus_i;
                    T new_o_vel;
                    T new_i_vel;
                    long double avg_rest;
                    for (const auto &[_, tup] : overlapping) {
                        o_vel = std::get<1>(tup)/bod_o->mass_; // recalculating is cheaper than adding them to the tuple
                        i_vel = std::get<2>(tup)/std::get<0>(tup)->mass_; // up above
                        o_minus_i = o_vel - i_vel; // v1 - v2
                        avg_rest = (bod_o->rest_c + std::get<0>(tup)->rest_c)/2;
                        if (o_vel > 0 || i_vel < 0) {
                            if (o_vel <= 0 && i_vel < 0)
                                if (o_vel < i_vel) // bodies are already separating
                                    continue; // case for bodies having passed through each other
                            if (o_vel > 0 && i_vel >= 0)
                                if (o_vel < i_vel)
                                    continue; // case for bodies having passed through each other
                            new_o_vel = (std::get<1>(tup) + std::get<2>(tup) -
                                         std::get<0>(tup)->mass_*avg_rest*o_minus_i)/(bod_o->mass_ +
                                                                                      std::get<0>(tup)->mass_);
                            new_i_vel = (std::get<1>(tup) + std::get<2>(tup) +
                                         bod_o->mass_*avg_rest*o_minus_i)/(bod_o->mass_ + std::get<0>(tup)->mass_);
                            bod_o->curr_vel += (new_o_vel - o_vel)*std::get<3>(tup);
                            std::get<0>(tup)->curr_vel += (new_i_vel - i_vel)*std::get<3>(tup);
                        }
                    }
                }
                ++it;
            }
        }
        // bool have_nan = false;
        template <bool mem, bool file>
        void s_coll() { // "simple" (ahem, ahem) collision detection and evolution | removed static from map
            std::map<long double, std::tuple<bod_t*, decltype(M{}*T{}), decltype(M{}*T{}), vector3D<T>>> overlapping;
            outer = 0;
            num_bods = bods.size();
            bod_t *bod_o;
            bod_t *bod_i;
            bod_t *merging_bod;
            R rad_dist{};
            long double dist;
#ifdef GREGSYS_MERGERS
            mom_t axis_mom_o;
            mom_t axis_mom_i;
            mom_t max_mom{min_tot_com_mom};
            mom_t curr_mom{};
            // long double min_dist = HUGE_VALL; // usually expands to infinity
            vec_size_t merging_index;
            bool first_merged;
            bod_t *merged;
            alignas(bod_t) char fake_body[sizeof(bod_t)];
            std::set<uint64_t> merged_ids;
            merged_ids.clear();
#else
            bod_o = bods.data();
#endif
            while (outer < num_bods) {
#ifdef GREGSYS_MERGERS
                first_merged = false;
                merging_bod = nullptr;
                bod_o = bods.data() + outer;
#endif
                overlapping.clear();
                bod_i = bod_o + 1; // to point to body just after bod_o
                for (inner = outer + 1; inner < num_bods; ++inner, ++bod_i) {
                    // bod_i = bods.data() + inner;
#ifdef GREGSYS_MERGERS
                    inner_loop:
                    if (inner == outer)
                        continue;
#endif
                    rad_dist = bod_o->radius + bod_i->radius;
                    vector3D<T> &&r12 = bod_i->curr_pos - bod_o->curr_pos;
                    dist = r12.magnitude(); // DEAL WITH ZERO DISTANCE CASE
                    /* if (std::isnan(dist) && !have_nan) {
                        std::cout << "THERE IS A NAN HERE!!!!!" << std::endl;
                        std::cout << "bod_o: " << *bod_o << std::endl;
                        std::cout << "bod_i: " << *bod_i << std::endl;
                        for (const auto &bodd : *this)
                            std::cout << bodd << std::endl;
                        abort();
                        have_nan = true;
                    } */
                    r12.x /= dist; r12.y /= dist; r12.z /= dist; // more efficient than calling normalise()
                    if (rad_dist > dist) {
                        // vector3D<T> &&com_vel = vel_com(bod_o, bod_i);
#ifdef GREGSYS_MERGERS
                        axis_mom_o = bod_o->momentum()*r12;
                        axis_mom_i = bod_i->momentum()*r12;
                        if((curr_mom = axis_mom_o - axis_mom_i - (bod_o->mass_ - bod_i->mass_)*
                                (vel_com(bod_o, bod_i)*r12)) >= max_mom) {
                            max_mom = curr_mom;
                            merging_bod = bod_i;
                            merging_index = inner;
                            continue;
                        }
                        if (merging_bod == nullptr) {
                            overlapping.emplace(std::piecewise_construct, std::forward_as_tuple(dist),
                                                std::forward_as_tuple(bod_i, axis_mom_o, axis_mom_i, r12));
                        }
#else
                        overlapping.emplace(std::piecewise_construct, std::forward_as_tuple(dist),
                                    std::forward_as_tuple(bod_i, bod_o->momentum()*r12, bod_i->momentum()*r12, r12));
#endif
                    }
                }
#ifdef GREGSYS_MERGERS
                if (merging_bod != nullptr) {
                    /* I have opted not to use my += overload for adding a body onto another and updating it (this
                     * would avoid performing 2 deletions from the std::vector) because the position, velocity and KE
                     * data for the body before the merger would be mixed with its new "self" after the merger. In
                     * addition, there would then be inconsistency between which bodies end up in the "deleted" bodies
                     * std::set (retaining all their history) and those that remain in the main std::vector. */
                    // bod_t &&merged = *bod_o + *merging_bod; // create body that is the result of the merger
                    if (!first_merged) {
                        merged = new(fake_body) bod_t{bod_o, merging_bod};//create body that is the result of the merger
                        if constexpr (memFreq) {
                            del_bods_m->emplace(steps, std::move(*bod_o));
                        }
                        else {
                            del_bods_n->emplace(std::move(*bod_o));
                        }
                        // bods.erase(bods.begin() + outer); // erase the merged bodies from the std::vector
                        // bod_o->~body();
                        bod_o = merged;
                        first_merged = true;
                    } else {
                        *bod_o += *merging_bod;
                    }
                    // must check that the merging body is not one from a previous merger at this time-step:
                    if (!merged_ids.contains(merging_bod->id)) {
                        if constexpr (memFreq) {
                            del_bods_m->emplace(steps, std::move(*merging_bod));
                        }
                        else {
                            del_bods_n->emplace(std::move(*merging_bod));
                        }
                    }
                    bods.erase(bods.begin() + merging_index);// - 1); // -1 because outer body was just deleted
                    // new position of primary merged body will be one less if merging body comes before
                    if (merging_index < outer) {
                        --outer;
                    }
                    --num_bods; // there is a net loss of 1 body
                    inner = 0;
                    bod_i = bods.data();
                    overlapping.clear();
                    merging_bod = nullptr;
                    goto inner_loop;
                }
#endif
                if (!overlapping.empty()) {
                    T o_vel;
                    T i_vel;
                    T o_minus_i;
                    T new_o_vel;
                    T new_i_vel;
                    long double avg_rest;
                    for (const auto &[_, tup] : overlapping) {
                        o_vel = std::get<1>(tup)/bod_o->mass_; // recalculating is cheaper than adding them to the tuple
                        i_vel = std::get<2>(tup)/std::get<0>(tup)->mass_; // up above
                        o_minus_i = o_vel - i_vel; // v1 - v2
                        avg_rest = (bod_o->rest_c + std::get<0>(tup)->rest_c)/2;
                        if (o_vel > 0 || i_vel < 0) {
                            if (o_vel <= 0 && i_vel < 0)
                                if (o_vel < i_vel) // bodies are already separating
                                    continue; // case for bodies having passed through each other
                            if (o_vel > 0 && i_vel >= 0)
                                if (o_vel < i_vel)
                                    continue; // case for bodies having passed through each other
                            new_o_vel = (std::get<1>(tup) + std::get<2>(tup) -
                                         std::get<0>(tup)->mass_*avg_rest*o_minus_i)/(bod_o->mass_ +
                                                                                      std::get<0>(tup)->mass_);
                            new_i_vel = (std::get<1>(tup) + std::get<2>(tup) +
                                         bod_o->mass_*avg_rest*o_minus_i)/(bod_o->mass_ + std::get<0>(tup)->mass_);
                            bod_o->curr_vel += (new_o_vel - o_vel)*std::get<3>(tup);
                            std::get<0>(tup)->curr_vel += (new_i_vel - i_vel)*std::get<3>(tup);
                        }
                    }
                }
#ifdef GREGSYS_MERGERS
                if (first_merged) {
                    // bod_o->rec_counter = steps;
                    // move the merged bodies into the deleted bodies set:
                    // bods.emplace_back(std::move(merged)); // move the new body into the std::vector storing all bodies
                    merged_ids.insert(bod_o->id);
                    bods[outer].~body(); // have to manually call dtor of body being moved into
                    copy(bods.data() + outer, fake_body, sizeof(bod_t)); // bypass constructors
                    // bod_o->~body_counter();
                    // if constexpr (file)
                    //     bod_o = bods.data() + outer;
                }
#endif
                // if constexpr (file) {
                //     // WRITE TO FILE
                // }
                ++outer;
#ifndef GREGSYS_MERGERS
                ++bod_o;
#endif
            }
            if constexpr (memFreq || file) { // repeated constexpr conditionals don't make it to runtime!!
                // the repeated loop below is to avoid re-checking the conditions repeatedly inside the loop
#ifdef GREGSYS_MERGERS
                if (!merged_ids.empty()) { // hmmmm... I think I must remove the negation, it was so long I don't recall
                    // printf("Merged IDs:\n");
                    // for (const auto &id : merged_ids)
                    //     std::cout << id << std::endl;
#endif
                    if constexpr (mem || file) {
                        for (bod_t &b : bods) {
                            if constexpr (mem) {
                                b.add_pos_vel_ke();
                            }
                            if constexpr (file) {
                                // WRITE TO FILE
                            }
                        }
                    }
                    return;
#ifdef GREGSYS_MERGERS
                }
                if (merged_ids.size() < pe.size()) {
                    for (bod_t &b : bods) {
                        if constexpr (mem) {
                            b.add_pos_vel_ke();
                        }
                        // bod_o->rec_counter = steps;
                        if constexpr (memFreq) { // can remove this if constexpr - will do it later
                            if (!merged_ids.contains(b.id)) {
                                pe.emplace(b.id, std::vector<T>{}); // add a std::vector to store its potential energies
                                energy.emplace(b.id, std::vector<T>{}); // same for its total energy (stores its own KE)
                            }
                        }
                        if constexpr (file) {
                            // WRITE TO FILE
                        }
                    }
                    return;
                }
                for (bod_t &b : bods) {
                    if constexpr (mem) {
                        b.add_pos_vel_ke();
                    }
                    // bod_o->rec_counter = steps;
                    if constexpr (memFreq) { // also must remove this one
                        pe.emplace(b.id, std::vector<T>{}); // add a std::vector to store its potential energies
                        energy.emplace(b.id, std::vector<T>{}); // same for its total energy (it stores its own KE)
                    }
                    if constexpr (file) {
                        // WRITE TO FILE
                    }
                }
#endif
            }
        }
#ifdef GREGSYS_MERGERS
    public: // the 3 following methods are for debugging - mark for removal
        void analysis() requires (memFreq != 0) {
            printf("---------------ANALYSIS--------------\nBods:\n");
            for (const auto& bod : bods) {
                std::cout << "\nBody: " << bod <<
                          "\nPositions size: " << bod.positions.size() <<
                          ", Velocity size: " << bod.velocities.size() <<
                          "\nKE size: " << bod.energies.size() << ", PE size: " << pe[bod.id].size() <<
                          ", Energy size: " << energy[bod.id].size() << '\n' << std::endl;
            }
            printf("-----------\nDel_Bods:\n");
            for (const auto& [it, bod] : *del_bods_m) {
                std::cout << "Iteration: " << it << "\nBody: " << bod <<
                "\nPositions size: " << bod.positions.size() <<
                ", Velocity size: " << bod.velocities.size() <<
                "\nKE size: " << bod.energies.size() << ", PE size: " << pe[bod.id].size() <<
                ", Energy size: " << energy[bod.id].size() << '\n' << std::endl;
            }
        }
        std::vector<bod_t>& get_bods() {
            return this->bods;
        }
        std::set<std::pair<uint64_t, bod_t>, decltype(pair_bods_func)>& get_del_bods() requires (memFreq != 0) {
            return *this->del_bods_m;
        }
#endif
    private:
        static inline uint64_t nsys_file_size(uint64_t n_bods) noexcept {
            char rem = (char) ((n_bods*sizeof(nsys_chunk)) % 4);
            return sizeof(nsys_header) + n_bods*sizeof(nsys_chunk) + (rem ? 4 - rem : 0);
        }
        // path is guaranteed not to be nullptr here:
        void load_from_nsys(const char *path, bool alloc_chunks, bool check_nan_inf)
        requires (std::is_fundamental_v<M> && std::is_fundamental_v<R> && std::is_fundamental_v<T> && CHAR_BIT == 8) {
            auto final_handler = [this](uint64_t index, const char *str) { // removed static
                this->istream->close();
                delete this->istream;
                char msg[128];
                snprintf(msg, 128, "Error: NaN or INF encountered in chunk with index %" PRIu64 // uint64_t format spec.
                ". Invalid field: \"%s\"\n", index, str);
                throw nsys_load_error{msg};
            };
#ifndef _WIN32
            /* Unfortunately, there is absolutely no POSIX-conforming standard way of ensuring that files larger than
             * 2GB can be worked with. The size of off_t can be checked to ensure it is 64 bits wide, but nothing can be
             * done if it is not. There are many non-standard extensions, such as defining the _FILE_OFFSET_BITS macro
             * as 64 (as a GNU extension), but as these are not part of the standard, I do not include them. */
            struct stat buffer{};
            if (stat(path, &buffer) == -1)
                throw nsys_load_error{"Error obtaining .nsys file information. No data loaded.\n"};
            if (!S_ISREG(buffer.st_mode))
                throw nsys_load_error{".nsys file provided is not a regular file. Could not load data.\n"};
            if (buffer.st_size < (off_t) sizeof(nsys_header))
                throw nsys_load_error{"Insufficient .nsys file size. No data loaded.\n"};
#else
            WIN32_FILE_ATTRIBUTE_DATA buffer{};
            if (!GetFileAttributesExA(path, GetFileExInfoStandard, &buffer))
                throw nsys_load_error{"Error obtaining .nsys file attributes. No n-body data loaded.\n"};
            if (buffer.dwFileAttributes != FILE_ATTRIBUTE_NORMAL)
                throw nsys_load_error{".nsys file provided is not a regular file. No data loaded.\n"};
            uint64_t fileSize = buffer.dwFileSizeLow + (buffer.dwFileSizeHigh << 32);
            if (fileSize < sizeof(nsys_header))
                throw nsys_load_error{".nsys file size insufficient. Could not load data."};
#endif
            istream = new std::ifstream{path, std::ios_base::in | std::ios_base::binary};
            if (!*istream) {
                delete istream;
                throw nsys_load_error{"Error opening .nsys file. Could not load data.\n"};
            }
            nsys_header header;
            istream->read((char *) &header, sizeof(nsys_header));
            if (header.signature[0] != 'N' && header.signature[1] != 'S') {
                istream->close();
                delete istream;
                throw nsys_format_error{"Invalid .nsys format: invalid file signature.\n"};
            }
            if ((header.d_types & 0b1000000) != 0) { // 0b10000000 == 128, but I use the binary repr. for clarity
                istream->close();
                delete istream;
                throw nsys_format_error{"Invalid .nsys format: msb occupied in d_types field (should be zero).\n"};
            } // the following macro saves about 60 lines of code
#define NSYS_TYPE_CHECK(type, bin_num1, bin_num2, str) \
            if constexpr (std::floating_point<type>) { \
                if (!(header.d_types & bin_num1)) {    \
                    istream->close(); \
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          " floating-point format, but template parameter "#type" is integral.\n"}; \
                } \
            } \
            else if constexpr (std::signed_integral<type>) { \
                if (header.d_types & bin_num1) {       \
                    istream->close(); \
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          " floating-point format, but template parameter "#type" is integral.\n"}; \
                }/* I have a separate conditional for the case below so that a different error message can be logged */\
                if (!(header.d_types & bin_num2)) {    \
                    istream->close(); \
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          " unsigned integral format, but template parameter "#type" is signed " \
                                          "integral.\n"}; \
                } \
            } \
            else { \
                if (header.d_types & bin_num1) {       \
                    istream->close();\
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          " floating-point format, but template parameter "#type" is integral.\n"}; \
                } \
                if (header.d_types & bin_num2) {       \
                    istream->close(); \
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          "signed integral format, but template parameter "#type" is unsigned " \
                                          "integral.\n"}; \
                } \
            } \
            if (header.type##_size != sizeof(type)) {  \
                istream->close(); \
                delete istream; \
                char error_msg[128]; \
                snprintf(error_msg, 128, "Size mismatch error: size of "#type" template parameter (%zu bytes) does " \
                                         "not match its reported size (%d bytes) in the .nsys file.\n", \
                        sizeof(type), header.type##_size); /* includes null terminator */ \
                throw nsys_load_error{error_msg}; \
            }
            NSYS_TYPE_CHECK(M, 0b00000010, 0b00010000, mass)
            NSYS_TYPE_CHECK(R, 0b00000100, 0b00100000, radius)
            NSYS_TYPE_CHECK(T, 0b00001000, 0b01000000, position and velocity)
#undef NSYS_TYPE_CHECK
            if (header.rest_coeff_size > sizeof(long double)) {
                istream->close();
                delete istream;
                char error_msg[156];
                snprintf(error_msg, 156, "Size mismatch error: coefficient of restitution data size (%d bytes) larger "
                                         "than widest floating point data type available (long double = %zu bytes).\n",
                         header.rest_coeff_size, sizeof(long double)); // includes null terminator
                throw nsys_format_error{error_msg};
            }
            void (*rest_c_func)(long double *, const char *);
            if (header.rest_coeff_size == sizeof(long double))
                rest_c_func = +[](long double *coeff, const char *ptr) {
                char *coeff_ptr = (char *) coeff;
                for (unsigned char counter = 0; counter < sizeof(long double); ++counter)
                    *coeff_ptr++ = *ptr++;
            };
            else {
                if (header.rest_coeff_size == sizeof(double)) {
                    rest_c_func = +[](long double *coeff, const char *ptr) {
                        double coeff_d; // removed static
                        char *coeff_ptr = (char *) &coeff_d;
                        for (unsigned char counter = 0; counter < sizeof(double); ++counter)
                            *coeff_ptr++ = *ptr++;
                        *coeff = coeff_d;
                    };
                } else if (header.rest_coeff_size == sizeof(float)) {
                    rest_c_func = +[](long double *coeff, const char *ptr) {
                        float coeff_f; // removed static
                        char *coeff_ptr = (char *) &coeff_f;
                        for (unsigned char counter = 0; counter < sizeof(float); ++counter)
                            *coeff_ptr++ = *ptr++;
                        *coeff = coeff_f;
                    };
                } else {
                    istream->close();
                    delete istream;
                    char error_msg[192];
                    snprintf(error_msg, 192, "Size mismatch error: coefficient of restitution data size (%d bytes) "
                                             "does not match any floating point data type available:\n"
                                             "float = %zu bytes,\n"
                                             "double = %zu bytes,\n"
                                             "long double = %zu bytes.\n",
                             header.rest_coeff_size, sizeof(float), sizeof(double), sizeof(long double));
                    throw nsys_format_error{error_msg};
                }
            }
#ifndef _WIN32
            if (buffer.st_size != (off_t) nsys_file_size(header.num_bodies)) {
#else
            if (fileSize != nsys_file_size(header.num_bodies)) {
#endif
                istream->close();
                delete istream;
                char error_msg[156];
                snprintf(error_msg, 156, "Error: file size does not match expected file size based on number of bodies."
                                         "\nExpected: %llu bytes. Instead got: %zu bytes. No data loaded.\n",
#ifndef _WIN32
                         (unsigned long long) nsys_file_size(header.num_bodies), (size_t) buffer.st_size);
#else
                         (unsigned long long) nsys_file_size(header.num_bodies), (size_t) fileSize);
#endif
                throw nsys_format_error{error_msg};
            }
            if (header.d_types & 0b00000001) { // case for floating-point time value
                if constexpr (sizeof(long double) > sizeof(double) && sizeof(double) > sizeof(float)) {
                    switch (header.time_size) {
                        case sizeof(long double):
#define NSYS_LDBL_TIME_CASE \
                            copy(&this->time_elapsed, header.time_point, sizeof(long double)); \
                            copy(&this->dt, header.delta_t, sizeof(long double));
                            NSYS_LDBL_TIME_CASE
                            break;
                        case sizeof(double): {
#define NSYS_DBL_TIME_CASE \
                            double time_val; \
                            copy(&time_val, header.time_point, sizeof(double)); \
                            this->time_elapsed = time_val; \
                            copy(&time_val, header.delta_t, sizeof(double)); \
                            this->dt = time_val;
                            NSYS_DBL_TIME_CASE
                            break;
                        }
                        case sizeof(float): { // almost certainly == 4
#define NSYS_FLT_TIME_CASE \
                            float time_val; \
                            copy(&time_val, header.time_point, sizeof(float)); \
                            this->time_elapsed = time_val; \
                            copy(&time_val, header.delta_t, sizeof(float)); \
                            this->dt = time_val;
                            NSYS_FLT_TIME_CASE
                            break;
                        }
                        default:
#define NSYS_TIME_SIZE_ERROR \
                            istream->close(); \
                            delete istream; \
                            char error_msg[212]; \
                            snprintf(error_msg, 212, "Size mismatch error: reported size of data type for time values" \
                                                     " in .nsys file (%d bytes) does not match any floating " \
                                                     "point data type available:\n" \
                                                     "float = %zu bytes,\n" \
                                                     "double = %zu bytes,\n" \
                                                     "long double = %zu bytes.\n", \
                                     header.time_size, sizeof(float), sizeof(double), sizeof(long double)); \
                            throw nsys_format_error{error_msg};
                            NSYS_TIME_SIZE_ERROR
                            // this->time_elapsed = std::numeric_limits<long double>::has_quiet_NaN ?
                            //                      std::numeric_limits<long double>::quiet_NaN() :
                            //                      std::numeric_limits<long double>::max(); // error case
                    }
                }
                else if constexpr (sizeof(long double) == sizeof(double) && sizeof(double) > sizeof(float)) {
                    if (header.time_size == sizeof(long double)) {
                        NSYS_LDBL_TIME_CASE
                    }
                    else if (header.time_size == sizeof(float)) {
                        NSYS_FLT_TIME_CASE
                    }
                    else {
                        NSYS_TIME_SIZE_ERROR
                    }
                }
                else if constexpr (sizeof(long double) > sizeof(double) && sizeof(double) == sizeof(float)) {
                    if (header.time_size == sizeof(long double)) {
                        NSYS_LDBL_TIME_CASE
                    }
                    else if (header.time_size == sizeof(double)) {
                        NSYS_DBL_TIME_CASE
                    }
                    else {
                        NSYS_TIME_SIZE_ERROR
                    }
                }
                else { // case for sizeof(long double) == sizeof(double) == sizeof(float)
                    if (header.time_size == sizeof(long double)) {
                        NSYS_LDBL_TIME_CASE
                    }
                    else {
                        NSYS_TIME_SIZE_ERROR
                    }
                }
#undef NSYS_TIME_SIZE_ERROR
#undef NSYS_FLT_TIME_CASE
#undef NSYS_DBL_TIME_CASE
#undef NSYS_LDBL_TIME_CASE
                /*
                switch (header.time_size) {
                    case sizeof(long double):
                        copy(&this->time_elapsed, header.time_point, sizeof(long double));
                        copy(&this->dt, header.delta_t, sizeof(long double));
                        break;
                    case sizeof(double): {
                        double time_val;
                        copy(&time_val, header.time_point, sizeof(double));
                        this->time_elapsed = time_val;
                        copy(&time_val, header.delta_t, sizeof(double));
                        this->dt = time_val;
                        break;
                    }
                    case sizeof(float): { // almost certainly == 4
                        float time_val;
                        copy(&time_val, header.time_point, sizeof(float));
                        this->time_elapsed = time_val;
                        copy(&time_val, header.delta_t, sizeof(float));
                        this->dt = time_val;
                        break;
                    }
                    default:
                        istream->close();
                        delete istream;
                        char error_msg[212];
                        snprintf(error_msg, 212, "Size mismatch error: reported size of data type for time values in "
                                                 ".nsys file (%d bytes) does not match any floating point data type "
                                                 "available:\n"
                                                 "float = %zu bytes,\n"
                                                 "double = %zu bytes,\n"
                                                 "long double = %zu bytes.\n",
                                 header.time_size, sizeof(float), sizeof(double), sizeof(long double));
                        throw nsys_format_error{error_msg};
                        // this->time_elapsed = std::numeric_limits<long double>::has_quiet_NaN ?
                        //                      std::numeric_limits<long double>::quiet_NaN() :
                        //                      std::numeric_limits<long double>::max(); // error case
                } */
            }
            else { // case for integral time value
                if constexpr (sizeof(unsigned long long) >= sizeof(uint64_t)) {
                    if (header.time_size <= sizeof(unsigned long long)) {
                        // this->time_elapsed = std::numeric_limits<long double>::has_quiet_NaN ?
                        //                      std::numeric_limits<long double>::quiet_NaN() :
                        //                      std::numeric_limits<long double>::max(); // error case
                        unsigned long long time_val = 0;
                        copy(&time_val, header.time_point, header.time_size);
                        this->time_elapsed = time_val;
                        copy(&time_val, header.delta_t, header.time_size);
                        this->dt = time_val;
                        goto correct;
                    }
                }
                else {
                    if (header.time_size <= sizeof(uint64_t)) {
                        // this->time_elapsed = std::numeric_limits<long double>::has_quiet_NaN ?
                        //                      std::numeric_limits<long double>::quiet_NaN() :
                        //                      std::numeric_limits<long double>::max(); // error case
                        uint64_t time_val = 0;
                        copy(&time_val, header.time_point, header.time_size);
                        this->time_elapsed = time_val;
                        copy(&time_val, header.delta_t, header.time_size);
                        this->dt = time_val;
                        goto correct;
                    }
                }
                istream->close();
                delete istream;
                char error_msg[204];
                snprintf(error_msg, 204, "Size mismatch error: reported size of data type for time values in "
                                         ".nsys file (%d bytes) larger than largest unsigned integral data types "
                                         "available:\n"
                                         "uint64_t = %zu bytes,\n"
                                         "unsigned long long = %zu bytes\n",
                         header.time_size, sizeof(uint64_t), sizeof(unsigned long long));
                throw nsys_format_error{error_msg};
            }
            correct:
            if (check_nan_inf && (std::isnan(this->dt) || std::isinf(this->dt))) {
                istream->close();
                delete istream;
                throw nsys_format_error{"Error: time-step is NaN or INF."};
            }
            /* The above check is not performed on the time_elapsed variable as it is not necessary for the correct
             * functioning of a gtd::system<> instance. */
            this->check_dt();
            // this->half_dt = this->dt/2;
            this->G = 66'743;
            this->G *= header.m_num*header.d_num*header.d_num*header.d_num*header.t_num*header.t_num;
            this->G /= header.m_denom*header.d_denom*header.d_denom*header.d_denom*header.t_denom*header.t_denom;
            this->G /= 1'000'000'000'000'000;
            *G_vals++ = header.m_num;
            *G_vals++ = header.m_denom;
            *G_vals++ = header.d_num;
            *G_vals++ = header.d_denom;
            *G_vals++ = header.t_num;
            *G_vals = header.t_denom;
            G_vals -= 5;
            this->bods.clear();
            if (!header.num_bodies) {
                istream->close();
                delete istream;
                return;
            }
            this->bods.reserve(header.num_bodies); // avoids memory reallocation during the loop
            long double restitution{};
            uint64_t counter = 0;
            if (alloc_chunks) {
                /* Option 1: read in all data from chunks into allocated memory, and construct bodies in the "bods"
                 * std::vector using the data from the chunks. Much faster than option 2, but uses much more memory. */
                nsys_chunk *chunks = new nsys_chunk[header.num_bodies];
                nsys_chunk *chunk = chunks;
                istream->read((char *) chunks, header.num_bodies*sizeof(nsys_chunk));
                istream->close();
                delete istream;
                if (check_nan_inf) {
                    while (counter < header.num_bodies) {
                        rest_c_func(&restitution, chunk->r_coeff);
                        if (std::isnan(restitution) || std::isinf(restitution))
                            final_handler(counter, "r_coeff");
                        if constexpr (std::is_floating_point_v<M>)
                            if (std::isnan(chunk->mass) || std::isinf(chunk->mass))
                                final_handler(counter, "mass");
                        if constexpr (std::is_floating_point_v<R>)
                            if (std::isnan(chunk->radius) || std::isinf(chunk->radius))
                                final_handler(counter, "radius");
                        if constexpr (std::is_floating_point_v<T>) {
                            if (std::isnan(chunk->xpos) || std::isinf(chunk->xpos))
                                final_handler(counter, "xpos");
                            if (std::isnan(chunk->ypos) || std::isinf(chunk->ypos))
                                final_handler(counter, "ypos");
                            if (std::isnan(chunk->zpos) || std::isinf(chunk->zpos))
                                final_handler(counter, "zpos");
                            if (std::isnan(chunk->xvel) || std::isinf(chunk->xvel))
                                final_handler(counter, "xvel");
                            if (std::isnan(chunk->yvel) || std::isinf(chunk->yvel))
                                final_handler(counter, "yvel");
                            if (std::isnan(chunk->zvel) || std::isinf(chunk->zvel))
                                final_handler(counter, "zvel");
                        }
                        /* I created a new, nasty constructor in the gtd::body<> class which accepts individual
                         * components of the position and velocity (rather than two vectors) to avoid making any
                         * unnecessary copies. */
                        this->bods.emplace_back(chunk->mass, chunk->radius, chunk->xpos, chunk->ypos, chunk->zpos,
                                                chunk->xvel, chunk->yvel, chunk->zvel, restitution);
                        ++chunk;
                        ++counter;
                    }
                    delete [] chunks;
                    return;
                }
                while (counter < header.num_bodies) {
                    rest_c_func(&restitution, chunk->r_coeff);
                    this->bods.emplace_back(chunk->mass, chunk->radius, chunk->xpos, chunk->ypos, chunk->zpos,
                                            chunk->xvel, chunk->yvel, chunk->zvel, restitution);
                    ++chunk;
                    ++counter;
                }
                delete [] chunks;
                return;
            }
            /* Option 2: read data in from the file chunk by chunk, constructing each gtd::body<> within the "bods"
             * std::vector after each read. Slower than option 1 as requires N times more I/O function calls (where N is
             * number of bodies), but reduces memory footprint. */
            nsys_chunk chunk;
            if (check_nan_inf) {
                while (counter < header.num_bodies) {
                    istream->read((char *) &chunk, sizeof(nsys_chunk));
                    rest_c_func(&restitution, chunk.r_coeff);
                    if (std::isnan(restitution) || std::isinf(restitution))
                        final_handler(counter, "r_coeff");
                    if constexpr (std::is_floating_point_v<M>)
                        if (std::isnan(chunk.mass) || std::isinf(chunk.mass))
                            final_handler(counter, "mass");
                    if constexpr (std::is_floating_point_v<R>)
                        if (std::isnan(chunk.radius) || std::isinf(chunk.radius))
                            final_handler(counter, "radius");
                    if constexpr (std::is_floating_point_v<T>) {
                        if (std::isnan(chunk.xpos) || std::isinf(chunk.xpos))
                            final_handler(counter, "xpos");
                        if (std::isnan(chunk.ypos) || std::isinf(chunk.ypos))
                            final_handler(counter, "ypos");
                        if (std::isnan(chunk.zpos) || std::isinf(chunk.zpos))
                            final_handler(counter, "zpos");
                        if (std::isnan(chunk.xvel) || std::isinf(chunk.xvel))
                            final_handler(counter, "xvel");
                        if (std::isnan(chunk.yvel) || std::isinf(chunk.yvel))
                            final_handler(counter, "yvel");
                        if (std::isnan(chunk.zvel) || std::isinf(chunk.zvel))
                            final_handler(counter, "zvel");
                    }
                    this->bods.emplace_back(chunk.mass, chunk.radius, chunk.xpos, chunk.ypos, chunk.zpos,
                                            chunk.xvel, chunk.yvel, chunk.zvel, restitution);
                    ++counter;
                }
            } else {
                while (counter < header.num_bodies) {
                    istream->read((char *) &chunk, sizeof(nsys_chunk));
                    rest_c_func(&restitution, chunk.r_coeff);
                    this->bods.emplace_back(chunk.mass, chunk.radius, chunk.xpos, chunk.ypos, chunk.zpos,
                                            chunk.xvel, chunk.yvel, chunk.zvel, restitution);
                    ++counter;
                }
            }
            istream->close();
            delete istream;
        }
        void check_dt() {
            if (this->dt <= 0)
                throw std::invalid_argument{"Time-step cannot be zero or negative.\n"};
            this->half_dt = this->dt/2;
        }
        static inline bool check_option(int option) noexcept {
            uint16_t loword = option & 0x0000ffff;
            uint16_t hiword = option >> 16;
            return !(loword & (loword - 1)) || !loword || loword > rk4 ||
                   !(hiword & (hiword - 1)) || hiword > 2;
        }
        // constexpr void check_coll() {
        //     static_assert(collisions <= pred_coll_check && !(collisions & (collisions + 1)),
        //                   "Invalid collision-checking option.");
        // }
        template <bool iters = true>
        void print_progress() const noexcept requires (prog) {
            if constexpr (iters)
#ifndef _WIN32
                printf(CYAN_TXT_START "Iteration " BLUE_TXT_START "%" PRIu64 RED_TXT_START "/" MAGENTA_TXT_START
                   "%" PRIu64"\r", this->steps, this->iterations);
#else
                printf("Iteration %" PRIu64"/%" PRIu64"\r", steps, iterations);
#endif
            else
#ifndef _WIN32
                printf(CYAN_TXT_START "Iteration " BLUE_TXT_START "%" PRIu64"\r", this->steps);
#else
                printf("Iteration %" PRIu64"\r", steps);
#endif
        }
        void print_conclusion(const std::chrono::time_point<std::chrono::high_resolution_clock> &start,
                              std::chrono::nanoseconds &total) requires (prog) {
            total = std::chrono::high_resolution_clock::now() - start;
            std::cout << RESET_TXT_FLAGS << BLACK_TXT("\n--------------------Done--------------------\n") <<
            UNDERLINED_TXT_START GREEN_TXT_START << this->method_str() <<
            RESET_TXT_FLAGS WHITE_TXT(" - time elapsed: ") BOLD_TXT_START YELLOW_TXT_START << total.count()/BILLION <<
            RESET_TXT_FLAGS WHITE_TXT_START " second" << "s"[total == std::chrono::seconds{1}] <<
            RESET_TXT_FLAGS << std::endl;
        }
    public:
        static_assert(collisions <= pred_coll_check && !(collisions & (collisions + 1)),
                      "Invalid collision-checking option.");
        typedef struct nsys_file_header {
            const char signature[2] = {'N', 'S'};
            const char d_types = 1 + (std::floating_point<M> << 1) +
                                     (std::floating_point<R> << 2) +
                                     (std::floating_point<T> << 3) +
                                     (std::signed_integral<M> << 4) +
                                     (std::signed_integral<R> << 5) +
                                     (std::signed_integral<T> << 6);
            const char time_size = sizeof(long double);
            char time_point[16]{};
            char delta_t[16]{};
            const unsigned char M_size = sizeof(M);
            const unsigned char R_size = sizeof(R);
            const unsigned char T_size = sizeof(T);
            const unsigned char rest_coeff_size = sizeof(long double);
            uint64_t m_num{};
            uint64_t m_denom{};
            uint64_t d_num{};
            uint64_t d_denom{};
            uint64_t t_num{};
            uint64_t t_denom{};
            uint64_t num_bodies{};
        } nsys_header; // size == 96 bytes
        typedef struct nsys_file_chunk {
            char r_coeff[16];
            M mass{};
            R radius{};
            T xpos{};
            T ypos{};
            T zpos{};
            T xvel{};
            T yvel{};
            T zvel{};
        } nsys_chunk;
        system() : G{G_SI} {
            this->parse_units_format("M1:1,D1:1,T1:1");
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        explicit system(const char *nsys_path, bool alloc_chunks = true, bool check_nan_inf = true,
                        bool allow_overlapping_bodies = true) {
            if (nsys_path == nullptr)
                throw std::invalid_argument{"nullptr passed as .nsys file path.\n"};
            this->load_from_nsys(nsys_path, alloc_chunks, check_nan_inf);
            if (!allow_overlapping_bodies)
                this->check_overlap();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        system(const std::initializer_list<bod_t> &list) : bods{list}, G{G_SI} {
            check_overlap();
            // check_coll();
            this->parse_units_format("M1:1,D1:1,T1:1");
            if constexpr (memFreq)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        system(std::initializer_list<bod_t> &&list) : bods{std::move(list)}, G{G_SI} {
            check_overlap();
            // check_coll();
            this->parse_units_format("M1:1,D1:1,T1:1");
            if constexpr (memFreq)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        explicit system(long double timestep, uint64_t num_iterations, const char *units_format) :
                dt{timestep}, iterations{num_iterations} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            // check_coll();
            this->check_dt();
            parse_units_format(units_format);
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        template <uint64_t mF>
        system(const std::vector<body<M, R, T, mF>> &bodies, long double timestep = 1,
               uint64_t num_iterations = 1'000, const char *units_format = "M1:1,D1:1,T1:1",
               bool allow_overlapping_bodies = false) :
               bods{bodies.begin(), bodies.end()}, dt{timestep}, iterations{num_iterations} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            this->check_dt();
            if (!allow_overlapping_bodies)
                check_overlap();
            // check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq && mF) // no need to clear bodies if the ones passed do not record history
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        system(std::vector<bod_t> &&bodies, long double timestep = 1,
               uint64_t num_iterations = 1'000, const char *units_format = "M1:1,D1:1,T1:1",
               bool allow_overlapping_bodies = false) :
               bods{std::move(bodies)}, dt{timestep}, iterations{num_iterations} {
            this->check_dt();
            if (!allow_overlapping_bodies)
                check_overlap();
            // check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        template <uint64_t mF>
        system(std::vector<body<M, R, T, mF>> &&bodies, long double timestep = 1,
               uint64_t num_iterations = 1'000, const char *units_format = "M1:1,D1:1,T1:1",
               bool allow_overlapping_bodies = false) :
               dt{timestep}, iterations{num_iterations} {
            this->check_dt();
            std::for_each(bodies.begin(), bodies.end(),
                          [this](body<M, R, T, mF> &b){bods.emplace_back(std::move(b));});
            if (!allow_overlapping_bodies)
                check_overlap();
            // check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq && mF)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        system(const sys_t &other) : bods{other.bods.begin(), other.bods.end()}, dt{other.dt}, half_dt{dt/2},
        iterations{other.iterations}, G{other.G} {
            check_overlap();
            // check_coll();
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        system(sys_t &&other) : bods{std::move(other.bods)}, dt{other.dt}, half_dt{dt/2}, iterations{other.iterations},
        G{other.G} {
            check_overlap();
            // check_coll();
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        /* copy constructors are made to only copy the variables seen below: it does not make sense for a new system
         * object (even when copy-constructed) to be "evolved" if it has not gone through the evolution itself */
        template <bool prg, bool mrg, int coll, uint64_t mF, uint64_t fF, bool bF>
        // so a system can be constructed from another without same checks
        system(const system<M, R, T, prg, mrg, coll, mF, fF, bF> &other,
               bool allow_overlapping_bodies = true) : bods{other.bods.begin(), other.bods.end()}, dt{other.dt},
               half_dt{dt/2}, iterations{other.iterations}, G{other.G} {
            if (!allow_overlapping_bodies)
                check_overlap();
            // check_coll();
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq && mF)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        template <bool prg, bool mrg, int coll, uint64_t fF, bool bF>
        system(system<M, R, T, prg, mrg, coll, memFreq, fF, bF> &&other, bool allow_overlapping_bodies = true) :
               bods{std::move(other.bods)}, dt{other.dt}, half_dt{dt/2}, iterations{other.iterations}, G{other.G} {
            if (!allow_overlapping_bodies)
                check_overlap();
            // check_coll();
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        template <bool prg, bool mrg, int coll, uint64_t mF, uint64_t fF, bool bF>
        system(system<M, R, T, prg, mrg, coll, mF, fF, bF> &&other, bool allow_overlapping_bodies = true) :
               dt{other.dt}, half_dt{dt/2}, iterations{other.iterations}, G{other.G} {
            std::for_each(other.bods.begin(), other.bods.end(),
                          [this](body<M, R, T, mF> &b){bods.emplace_back(std::move(b));});
            if (!allow_overlapping_bodies)
                check_overlap();
            // check_coll();
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq && mF)
                clear_bodies();
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                this->set_set();
#endif
        }
        vec_size_t num_bodies() const noexcept {
            return bods.size();
        }
        long double elapsed_time() const noexcept {
            return this->time_elapsed; // will either be NaN or LDBL_MAX if bad value read in from .nsys file
        }
        long double G_val() const noexcept {
            return this->G;
        }
        uint64_t iters() const noexcept {
            return this->iterations;
        }
        uint64_t iters(uint64_t number) noexcept {
            if (!number) // cannot have zero iterations
                return this->iterations;
            this->iterations = number;
            return number;
        }
        long double timestep() const noexcept {
            return this->dt;
        }
        long double timestep(const long double &delta_t) noexcept {
            if (delta_t <= 0)
                return this->dt;
            this->dt = delta_t;
            this->half_dt = dt/2;
            return delta_t;
        }
        void set_G_units(const char *units_format) {
            this->parse_units_format(units_format);
        } // the size the .nsys file would be at the time the function is called:
        uint64_t nsys_file_size() requires (std::is_fundamental_v<M> &&
                                            std::is_fundamental_v<R> &&
                                            std::is_fundamental_v<T> && CHAR_BIT == 8) {
            return nsys_file_size(this->bods.size());
        }
#ifdef GREGSYS_SOFTENING
        long double softening() const noexcept {
            return eps;
        }
        long double softening(long double _epsilon) noexcept {
            if (_epsilon < sqrtl(std::numeric_limits<long double>::min())) // to be improved
                return eps;
            eps = _epsilon;
            eps_sq = _epsilon*_epsilon; // perhaps check this too
            return _epsilon;
        }
#endif
        long double bh_opening_angle(long double _theta) const noexcept {
            return bh_theta;
        }
        long double bh_opening_angle(long double _theta) noexcept {
            if (_theta < 0 || _theta >= PI)
                return bh_theta;
            bh_theta = _theta;
            return _theta;
        }
#ifdef GREGSYS_MERGERS
        bool set_min_tot_com_mom(const mom_t& total_momentum) requires (collisions != 0) {
            if (total_momentum < mom_t{})
                return false;
            this->min_tot_com_mom = total_momentum;
            return true;
        }
#endif
        sys_t &add_body(const bod_t &bod, bool chk_overlap = true) {
            bods.emplace_back(bod);
            if constexpr (memFreq)
                bods.back().clear(); // all bodies within a system object must start out without an evolution
            if (chk_overlap)
                check_overlap();
            return *this;
        }
        sys_t &add_body(bod_t &&bod, bool chk_overlap = true) {
            if constexpr (memFreq)
                bod.clear();
            bods.emplace_back(std::move(bod));
            if (chk_overlap)
                check_overlap();
            return *this;
        }
        template <typename ...Args>
        sys_t &emplace_body(bool chk_overlap, Args&& ...args) { // to allow a body to be constructed in-place
            bods.emplace_back(std::forward<Args>(args)...);
            if (chk_overlap)
                this->check_overlap();
            return *this;
        }
        sys_t &add_bodies(const std::vector<bod_t> &bodies, bool chk_overlap = true) {
            if (bodies.empty())
                return *this;
            num_bods = bods.size();
            bods.insert(bods.end(), bodies.begin(), bodies.end());
            if constexpr (memFreq)
                clear_bodies(num_bods);
            if (chk_overlap)
                this->check_overlap();
            return *this;
        }
        sys_t &add_bodies(std::vector<bod_t> &&bodies, bool chk_overlap = true) {
            if (bodies.empty())
                return *this;
            for (auto &b : bodies) {
                if constexpr (memFreq)
                    b.clear();
                bods.emplace_back(std::move(b));
            }
            if (chk_overlap)
                this->check_overlap();
            return *this;
        }
        template <bool progress, bool merging, int coll, uint64_t mF, uint64_t fF, bool bin>
        sys_t &add_system(const system<M, R, T, progress, merging, coll, mF, fF, bin> &other, bool chk_overlap = true) {
            return this->add_bodies(other.bods, chk_overlap);
        }
        template <bool progress, bool merging, int coll, uint64_t mF, uint64_t fF, bool bin>
        sys_t &add_system(system<M, R, T, progress, merging, coll, mF, fF, bin> &&other, bool chk_overlap = true) {
            return this->add_bodies(std::move(other.bods), chk_overlap);
        }
        bod_t &back() {
            if (!this->bods.empty())
                return this->bods.back();
            throw std::logic_error{"gtd::system<>::back cannot be called on an empty gtd::system object.\n"};
        }
        const bod_t &back() const {
            if (!this->bods.empty())
                return this->bods.back();
            throw std::logic_error{"gtd::system<>::back cannot be called on an empty gtd::system object.\n"};
        }
        bod_t &front() {
            if (!this->bods.empty())
                return this->bods.front();
            throw std::logic_error{"gtd::system<>::front cannot be called on an empty gtd::system object.\n"};
        }
        const bod_t &front() const {
            if (!this->bods.empty())
                return this->bods.front();
            throw std::logic_error{"gtd::system<>::front cannot be called on an empty gtd::system object.\n"};
        }
        bod_t &get_body(uint64_t id) {
            for (auto &b : *this) // this is where it would be nice to be using a set!!
                if (b.id == id)
                    return b;
            throw std::invalid_argument("The id passed does not correspond to any body present in this system "
                                        "object.\n");
        }
        const bod_t &get_body(uint64_t id) const {
            for (const auto &b : *this)
                if (b.id == id)
                    return b;
            throw std::invalid_argument("The id passed does not correspond to any body present in this system "
                                        "object.\n");
        }
        bool remove_body(uint64_t id) {
            auto end_it = bods.cend();
            for (typename std::vector<bod_t>::const_iterator it = bods.cbegin(); it < end_it; ++it)
                if (it->id == id) {
                    bods.erase(it);
                    return true;
                }
            return false;
        }
        ke_t kinetic_energy() const noexcept {
            ke_t tot{};
            for (const bod_t &b : this->bods)
                tot += 0.5l*b.mass_*(b.curr_vel*b.curr_vel);
            return tot;
        }
        pe_t potential_energy() const noexcept {
            if (this->bods.empty())
                return {};
            pe_t tot{};
            const bod_t *obod = this->bods.data();
            const bod_t *ibod;
            uint64_t o_counter = 0;
            uint64_t i_counter;
            this->num_bods = this->bods.size();
            uint64_t num_m1 = this->num_bods - 1;
            while (o_counter < num_m1) {
                i_counter = ++o_counter;
                ibod = obod + 1;
                while (i_counter++ < this->num_bods)
                    tot -= (obod->mass_*ibod->mass_)/(obod->curr_pos - ibod->curr_pos).magnitude();
            }
            return this->G*tot;
        }
        decltype(std::declval<ke_t>() + std::declval<pe_t>()) total_energy() const noexcept {
            return this->kinetic_energy() + this->potential_energy();
        }
        vec_t com_pos() const {
            if (this->bods.empty())
                throw empty_system_error{};
            vec_t _pos;
            M mass{};
            for (const bod_t &b : this->bods) {
                _pos += b.mass_*b.curr_pos;
                mass += b.mass_;
            }
            return _pos/mass;
        }
        vec_t com_vel() const {
            if (this->bods.empty())
                throw empty_system_error{};
            vec_t _vel;
            M mass{};
            for (const bod_t &b : this->bods) {
                _vel += b.mass_*b.curr_vel;
                mass += b.mass_;
            }
            return _vel/mass;
        }
        std::pair<vec_t, vec_t> com_pos_vel() const { // convenience method, avoids double mass-summation
            if (this->bods.empty())
                throw empty_system_error{};
            vec_t _pos;
            vec_t _vel;
            M mass{};
            for (const bod_t &b : this->bods) {
                _pos += b.mass_*b.curr_pos;
                _vel += b.mass_*b.curr_vel;
                mass += b.mass_;
            }
            return {_pos/mass, _vel/mass};
        }
        static inline sys_t random_comet(const vec_t &pos, // I will likely remove all these newlines eventually, but,
                                         const vec_t &vel, // for now, I like the clarity
                                         const vec_t &_omega,
                                         uint64_t _n,
                                         const T &bounding_rad,
                                         const M &b_mass,
                                         const R &b_rad,
                                         long double r_coeff = 1.0l,
                                         int integration_method = sys_t::leapfrog_kdk,
                                         long double ev_r_coeff = 0.5l,
                                         uint64_t ev_iters = 10'000,
                                         long double ev_dt = 0.25,
                                         uint64_t iters = 1'000,
                                         long double dt = 1,
                                         const char *units_format = "M1:1,D1:1,T1:1")
                                         requires (std::convertible_to<T, long double>) {
            if (bounding_rad <= 0)
                throw negative_radius_error{"Error: bounding radius must be larger than zero.\n"};
            if (bounding_rad < b_rad)
                throw std::invalid_argument{"Error: the bounding radius must exceed the radius of a body.\n"};
            if (SPHERE_VOLUME(bounding_rad) < SPHERE_VOLUME(b_rad)*_n*4/HCP_PACKING_FRACTION)
                throw std::invalid_argument{"Error: insufficient bounding sphere volume for placement of bodies.\n"};
#define RAND_COM_CHECKS(coefficient) \
            using vec = vector3D<T>; \
            if (!_n) \
                throw std::invalid_argument{"Error: number of bodies cannot be zero.\n"}; \
            if (b_mass <= 0) \
                throw negative_radius_error{"Error: the mass each body will have must be larger than zero.\n"}; \
            if (b_rad <= 0) \
                throw negative_radius_error{"Error: the radius each body will have must be larger than zero.\n"}; \
            if (dt <= 0 || ev_dt <= 0) \
                throw std::invalid_argument{"Error: time-steps cannot be negative or zero.\n"}; \
            if (r_coeff < 0 || r_coeff > 1 || coefficient < 0 || coefficient > 1) \
                throw std::invalid_argument{"Error: coefficients of restitution must be between 0 and 1.\n"}; \
            check_option(integration_method);
            RAND_COM_CHECKS(ev_r_coeff)
            std::uniform_real_distribution<long double> _r{0.0l, 1.0l};
#define RAND_COM_BODY(com1, com2, com3, ...) \
            const uint64_t n_ = _n; \
            std::uniform_real_distribution<long double> _phi{0, 2*PI}; \
            std::uniform_real_distribution<long double> _theta{0, 1}; \
            std::mt19937_64 _mtw{seed()}; /* mersenne-twister engine */ \
            std::vector<bod_t> _bods; \
            _bods.reserve(_n); \
            vec _bpos; \
            long double r_val; \
            long double phi_val; \
            long double theta_val; \
            long double sin_theta; \
            bod_t *const cptr = _bods.data(); \
            bod_t *ptr; \
            uint64_t counter; \
            uint64_t curr_size = 0; \
            R b_diam = 2*b_rad; \
            com1 \
            while (_n --> 0) { \
                loop_start: \
                r_val = bounding_rad*cbrtl(_r(_mtw)); \
                phi_val = _phi(_mtw); \
                theta_val = acosl(1 - 2*_theta(_mtw)); \
                sin_theta = sinl(theta_val); \
                _bpos.x = r_val*sin_theta*cosl(phi_val); \
                _bpos.y = r_val*sin_theta*sinl(phi_val); \
                _bpos.z = r_val*cosl(theta_val); \
                ptr = cptr; /* called reserve() at the beginning, guarantees no new allocations (since size < cap.) */ \
                counter = 0; \
                while (counter++ < curr_size) { \
                    if (vec_ops::distance(_bpos, ptr->curr_pos) < b_diam) \
                        goto loop_start; \
                    ++ptr; \
                } \
                _bods.emplace_back(__VA_ARGS__); \
                com2 \
                ++curr_size; \
            } \
            com3 /* all bodies are same mass, so mass can be taken out of COM equation */
#define RAND_COM_END(loop_bit, ...) \
            retsys.set_timestep(dt); \
            retsys.set_iterations(iters); \
            spin_bodies(retsys.bods, _com, _omega, true);\
            _com = pos - _com; /* No longer COM, but offset between COM and comet position */ \
            for (bod_t &_b : retsys) {    \
                loop_bit \
                _b.rest_c = r_coeff; \
                _b.curr_pos += _com; \
                _b.curr_vel += vel; \
            } \
            __VA_ARGS__
            RAND_COM_BODY(vec _com{};, _com += _bpos;, _com /= n_;, b_mass, b_rad, _bpos, vec{}, ev_r_coeff)
            sys_t retsys{std::move(_bods), ev_dt, ev_iters, units_format, true};
            retsys.evolve(integration_method);
            RAND_COM_END(EMPTY, return retsys;)
        }
        // GAUSSIAN COMETS ARE DOWN FOR MAINTENANCE AT THE MOMENT
        /*
        static inline sys_t random_comet(const vec_t &pos,
                                         const vec_t &vel,
                                         const vec_t &_omega,
                                         long double sd, // sd cannot be after _n, else the method call would be ambig.
                                         uint64_t _n,
                                         const M &b_mass,
                                         const R &b_rad,
                                         long double r_coeff,
                                         int integration_method = sys_t::leapfrog_kdk,
                                         long double ev_r_coeff = 0.5l,
                                         uint64_t ev_iters = 10'000,
                                         long double ev_dt = 0.25,
                                         uint64_t iters = 1'000,
                                         long double dt = 1,
                                         const char *units_format = "M1:1,D1:1,T1:1")
        requires (std::convertible_to<T, long double>) {
            if (sd <= 0)
                throw std::invalid_argument{"Error: standard deviation must be positive.\n"};
            RAND_COM_CHECKS(ev_r_coeff)
            std::normal_distribution<long double> _r{0, sd};
            RAND_COM_BODY(vec _com{};, _com += _bpos;, _com /= n_;, b_mass, b_rad, _bpos, vec{}, ev_r_coeff)
            sys_t retsys{std::move(_bods), ev_dt, ev_iters, units_format, true};
            retsys.evolve(integration_method);
            RAND_COM_END(EMPTY, return retsys;)
        } */
        template <bool simple_rad = false>
        static inline std::tuple<sys_t, long double, T> // returns comet as system, packing fraction, effective radius
                            random_comet(const vec_t &pos,
                                         const vec_t &vel,
                                         const vec_t &_omega,
                                         uint64_t _n,
                                         const T &bounding_rad,
                                         const M &b_mass,
                                         const R &b_rad,
                                         long double r_coeff,
                                         long double n_exp,
                                         long double d_scale = 0,
                                         int integration_method = sys_t::leapfrog_kdk,
                                         long double min_cor = 0.5l,
                                         long double max_cor = 1.0l,
                                         long double ev_dt = 0.25l,
                                         uint64_t probing_iters = 1'000,
                                         long double dist_tol = 0,
                                         uint64_t iters = 1'000,
                                         long double dt = 1,
                                         const char *units_format = "M1:1,D1:1,T1:1") {
            if (bounding_rad <= 0)
                throw negative_radius_error{"Error: bounding radius must be larger than zero.\n"};
            if (bounding_rad < b_rad)
                throw std::invalid_argument{"Error: the bounding radius must exceed the radius of a body.\n"};
            R small_vol;
            if (SPHERE_VOLUME(bounding_rad) < (small_vol = SPHERE_VOLUME(b_rad))*_n*4/HCP_PACKING_FRACTION)
                throw std::invalid_argument{"Error: insufficient bounding sphere volume for placement of bodies.\n"};
            if (max_cor > 1 || max_cor <= min_cor)
                throw std::invalid_argument{"Error: maximum COR must be greater than minimum COR.\n"};
            RAND_COM_CHECKS(min_cor)
            std::uniform_real_distribution<long double> _r{0.0l, 1.0l};
            RAND_COM_BODY(vec _com{};, _com += _bpos;, _com /= n_;, b_mass, b_rad, _bpos, vec{})
            long double max_dist = cbrtl(((small_vol*n_)/SCP_PACKING_FRACTION)*(3/(4*_PI_)));
            if (d_scale <= 0) d_scale = max_dist;
            long double scale_factor = powl(d_scale, n_exp);
            max_dist += dist_tol*max_dist;
            long double max_dist_exp = powl(max_dist, n_exp);
            long double min_max_cor = max_cor - (max_cor - min_cor)*(max_dist_exp/(max_dist_exp + scale_factor));
#ifdef GREGSYS_MERGERS
            M tot_mass;
#endif
            long double dist_exp;
            for (bod_t &_b : _bods) {
                dist_exp = powl((_b.curr_pos - _com).magnitude(), n_exp);
                _b.rest_c = max_cor - (max_cor - min_cor)*(dist_exp/(dist_exp + scale_factor));
#ifdef GREGSYS_MERGERS
                tot_mass += _b.mass_;
#endif
            }
            sys_t retsys{std::move(_bods), ev_dt, iters, units_format, true};
            uint64_t _count = 0; // haven't made it static inside lambda to avoid overhead of setting to zero after ev.
            if constexpr (prog)
                std::cout << '\n';
            retsys.evolve(integration_method,
            [&retsys, &probing_iters, &min_max_cor, &_count](){
                if constexpr (prog) {
                    long double min_cor = 1.0l;
                    for (const bod_t &_b: retsys.bods) {
                        if (_b.rest_c < min_cor) {
                            min_cor = _b.rest_c;
                        }
                    }
                    std::cout <<
                              #ifndef _WIN32
                              "\033[F" // move cursor up by one line
                              #endif
                    BOLD_TXT_START BLUE_TXT_START "COR: " RESET_TXT_FLAGS
                    BOLD_TXT_START GREEN_TXT_START << min_cor << RESET_TXT_FLAGS BLACK_TXT_START "\t->\t"
                    RESET_TXT_FLAGS BOLD_TXT_START MAGENTA_TXT_START << min_max_cor << RESET_TXT_FLAGS "\n";
                    if (min_cor < min_max_cor) {
                        _count = 0;
                        return true;
                    }
                }
                else
                    for (const bod_t &_b : retsys.bods) {
                        if (_b.rest_c < min_max_cor) { // means at least 1 body is still too far out from COM
                            _count = 0;
                            return true; // so continue with iterations
                        }
                    }
                return _count++ < probing_iters;
            }, nullptr,
#ifndef GREGSYS_MERGERS
            [&min_cor, &max_cor, &scale_factor, &_com, &n_, &n_exp](auto &bodies){
#else
            [&min_cor, &max_cor, &scale_factor, &tot_mass, &_com, &n_, &n_exp](auto &bodies){
#endif
                // com is recalculated at each iteration to account for any drift caused by f.p. errors
                for (const bod_t &_b : bodies)
#ifndef GREGSYS_MERGERS
                    _com += _b.curr_pos;
                _com /= n_;
#else
                    _com += _b.mass_*_b.curr_pos;
                _com /= tot_mass;
#endif
                long double dist_exp;
                for (bod_t &_b : bodies) {
                    dist_exp = powl((_b.curr_pos - _com).magnitude(), n_exp);
                    _b.rest_c = max_cor - (max_cor - min_cor)*(dist_exp/(dist_exp + scale_factor));
                }
            });
            vec furthest;
            min_max_cor = 1;
            long double pf;
            if constexpr (simple_rad) {
                RAND_COM_END(if (_b.rest_c < min_max_cor) {min_max_cor = _b.rest_c; furthest = _b.curr_pos;},
                             T s_rad = (furthest - _com).magnitude() + b_rad;
                             return {std::move(retsys), (small_vol*n_)/SPHERE_VOLUME(s_rad), s_rad};)
            } else {
                for (const bod_t &_b : retsys.bods)
                    if (_b.rest_c < min_max_cor) {min_max_cor = _b.rest_c; furthest = _b.curr_pos;}
                pf = comet_pf(retsys, small_vol, furthest, _com);
                RAND_COM_END(EMPTY, EMPTY)
                return {std::move(retsys), pf, cbrtl((3*n_*small_vol)/(4*_PI_*pf))};
            }
        }
    private:
        static inline long double comet_pf(const sys_t &sys, const T &b_vol, const vec_t &furthest, const vec_t &_com) {
            uint64_t _inner_count = 0;
            uint64_t tot_count = 0;
            std::multiset<T> _outer;
            T lower = 0;
            T upper = vec_ops::distance(_com, furthest);
            T middle = upper/2.0l;
            T dist;
            for (const bod_t &_b : sys.bods) {
                dist = vec_ops::distance(_com, _b.curr_pos);
                if (dist <= middle) {
                    ++_inner_count;
                    continue;
                }
                _outer.insert(dist);
            }
            std::pair<typename std::multiset<T>::iterator, typename std::multiset<T>::iterator> it_pair;
            while (_inner_count < _outer.size()) {
                tot_count += _inner_count;
                _inner_count = 0;
                lower = middle;
                middle = (lower + upper)/2.0l;
                it_pair = _outer.equal_range(middle);
                if (!_outer.contains(middle))
                    ++it_pair.second;
                _inner_count = std::distance(_outer.begin(), it_pair.first);
                _outer.erase(_outer.begin(), it_pair.second);
            }
            tot_count += _inner_count;
            return (b_vol*tot_count)/SPHERE_VOLUME(middle);
        } /*
        static inline long double comet_pf(const sys_t &sys, const T &dr, const vec_t &_com) { // comet packing fraction
            std::vector<std::pair<T, uint64_t>> shells;
            T _prad = 0;
            T _rad = dr;
            uint64_t count;
            T dist;
            goto start_loop;
            while (count) {
                shells.emplace_back(_rad, count);
                start_loop: // to avoid placing zero into std::vector for final iteration
                count = 0;
                for (const bod_t &_b : sys.bods) {
                    dist = vec_ops::distance(_com, _b.curr_pos);
                    if (_prad <= dist && dist < _rad) {
                        ++count;
                    }
                }
                _prad = _rad;
                _rad += dr;
            }
            uint64_t count2;
        } */
    public:
        /*
        template <bool simple_rad = false>
        static inline std::tuple<sys_t, long double, T>
                            random_comet(const vec_t &pos,
                                         const vec_t &vel,
                                         const vec_t &_omega,
                                         long double sd,
                                         uint64_t _n,
                                         const M &b_mass,
                                         const R &b_rad,
                                         long double r_coeff,
                                         long double n_exp,
                                         long double d_scale,
                                         int integration_method = sys_t::leapfrog_kdk,
                                         long double min_cor = 0.5l,
                                         long double max_cor = 1.0l,
                                         long double ev_dt = 0.25l,
                                         uint64_t probing_iters = 1'000,
                                         long double dist_tol = 0,
                                         uint64_t iters = 1'000,
                                         long double dt = 1,
                                         const char *units_format = "M1:1,D1:1,T1:1") {
            if (sd <= 0)
                throw std::invalid_argument{"Error: standard deviation must be positive.\n"};
            RAND_COM_CHECKS(min_cor)
            std::normal_distribution<long double> _r{0, sd};
            RAND_COM_BODY(vec _com{};, _com += _bpos;, _com /= n_;, b_mass, b_rad, _bpos, vec{})
            R small_vol = SPHERE_VOLUME(b_rad); // is only a redundant variable if simple_rad == true
            long double max_dist = cbrtl(((small_vol*n_)/SCP_PACKING_FRACTION)*(3/(4*_PI_)));
            // MUST EVENTUALLY MOVE THE DUPLICATED CODE INTO A FUNCTION!!!
            if (d_scale <= 0) d_scale = max_dist;
            long double scale_factor = powl(d_scale, n_exp);
            max_dist += dist_tol*max_dist;
            long double max_dist_exp = powl(max_dist, n_exp);
            long double min_max_cor = max_cor - (max_cor - min_cor)*(max_dist_exp/(max_dist_exp + scale_factor));
#ifdef GREGSYS_MERGERS
            M tot_mass;
#endif
            long double dist_exp;
            for (bod_t &_b : _bods) {
                dist_exp = powl((_b.curr_pos - _com).magnitude(), n_exp);
                _b.rest_c = max_cor - (max_cor - min_cor)*(dist_exp/(dist_exp + scale_factor));
#ifdef GREGSYS_MERGERS
                tot_mass += _b.mass_;
#endif
            }
            sys_t retsys{std::move(_bods), ev_dt, iters, units_format, true};
            uint64_t _count = 0;
            if constexpr (prog)
                std::cout << '\n';
            retsys.evolve(integration_method,
                          [&retsys, &probing_iters, &min_max_cor, &_count](){
                              if constexpr (prog) {
                                  long double min_cor = 1.0l;
                                  for (const bod_t &_b: retsys.bods) {
                                      if (_b.rest_c < min_cor) {
                                          min_cor = _b.rest_c;
                                      }
                                  }
                                  std::cout <<
                                            #ifndef _WIN32
                                            "\033[F"
                                            #endif
                                            BOLD_TXT_START BLUE_TXT_START "COR: " RESET_TXT_FLAGS
                                            BOLD_TXT_START GREEN_TXT_START << min_cor << RESET_TXT_FLAGS BLACK_TXT_START
                                            "\t->\t" RESET_TXT_FLAGS BOLD_TXT_START MAGENTA_TXT_START << min_max_cor <<
                                            RESET_TXT_FLAGS "\n";
                                  if (min_cor < min_max_cor) {
                                      _count = 0;
                                      return true;
                                  }
                              }
                              else
                                  for (const bod_t &_b : retsys.bods) {
                                      if (_b.rest_c < min_max_cor) {
                                          _count = 0;
                                          return true;
                                      }
                                  }
                              return _count++ < probing_iters;
                          }, nullptr,
#ifndef GREGSYS_MERGERS
                          [&min_cor, &max_cor, &scale_factor, &_com, &n_, &n_exp](auto &bodies){
#else
                              [&min_cor, &max_cor, &scale_factor, &tot_mass, &_com, &n_, &n_exp](auto &bodies){
#endif
                              for (const bod_t &_b : bodies)
#ifndef GREGSYS_MERGERS
                                  _com += _b.curr_pos;
                              _com /= n_;
#else
                                _com += _b.mass_*_b.curr_pos;
                              _com /= tot_mass;
#endif
                              long double dist_exp;
                              for (bod_t &_b : bodies) {
                                  dist_exp = powl((_b.curr_pos - _com).magnitude(), n_exp);
                                  _b.rest_c = max_cor - (max_cor - min_cor)*(dist_exp/(dist_exp + scale_factor));
                              }
                          });
            vec furthest;
            min_max_cor = 1;
            long double pf;
            if constexpr (simple_rad) {
                RAND_COM_END(if (_b.rest_c < min_max_cor) {min_max_cor = _b.rest_c; furthest = _b.curr_pos;},
                             T s_rad = (furthest - _com).magnitude() + b_rad;
                             return {std::move(retsys), (small_vol*n_)/SPHERE_VOLUME(s_rad), s_rad};)
            } else {
                for (const bod_t &_b : retsys.bods)
                    if (_b.rest_c < min_max_cor) {min_max_cor = _b.rest_c; furthest = _b.curr_pos;}
                pf = comet_pf(retsys, small_vol, furthest, _com);
                RAND_COM_END(EMPTY, EMPTY)
                return {std::move(retsys), pf, cbrtl((3*n_*small_vol)/(4*_PI_*pf))};
            }
        } */
#undef RAND_COM_BODY
#undef RAND_COM_CHECKS
        static inline sys_t hcp_comet(const R &radius, long double bulk_density) {

        }
        template <bool strict_rad = false>
        static inline std::pair<sys_t, T>
                            hcp_comet(const T &rad,
                                      const vector3D<T> &pos,
                                      const vector3D<T> &vel,
                                      const R &spacing,
                                      const M &b_mass,
                                      const R &b_rad,
                                      long double r_coeff,
                                      const vec_t &orientation = vec_t::zero,
                                      const vec_t &_omega = vec_t::zero,
                                      bool force_central_body = true,
                                      bool adjust_to_com = true,
                                      long double timestep = 1,
                                      uint64_t num_iterations = 1'000,
                                      const char *units_format = "M1:1,D1:1,T1:1") {
            using vec = vector3D<T>;
            if (rad <= 0)
                throw negative_radius_error{"Error: comet radius must be larger than zero.\n"};
            if (b_mass <= 0)
                throw negative_mass_error{"Error: the mass each body will have must be larger than zero.\n"};
            if (b_rad <= 0)
                throw negative_radius_error{"Error: the radius each body will have must be larger than zero.\n"};
            if (spacing < 0)
                throw std::invalid_argument{"Error: spacing between bodies cannot be negative.\n"};
            if (rad < b_rad)
                throw std::invalid_argument{"Error: the radius of the comet must exceed the radius of a body.\n"};
            if (timestep <= 0)
                throw std::invalid_argument{"Error: time-step cannot be negative or zero.\n"};
            if (r_coeff < 0 || r_coeff > 1)
                throw std::invalid_argument{"Error: coefficient of restitution must be between 0 and 1.\n"};
            const R b_sep = b_rad*2 + spacing; // separation between the centres of adjacent bodies
            const R b_halfsep = b_rad + spacing/2.0l; // haven't done b_sep/2 to minimise f.p. error
            const T diaml = 2*(rad + rad*0.0625l); // slightly larger diameter to create hcp cube
            uint64_t lbods = llroundl(diaml / (b_sep)) + 1; // number of bodies along length of cube
            uint64_t whbods = llroundl((diaml*sqrtl(3)) / b_sep) + 1; // along width and along height of cube
            // std::cout << "b_sep: " << b_sep << std::endl;
            // std::cout << "b_halfsep: " << b_halfsep << std::endl;
            // std::cout << "radl: " << diaml << std::endl;
            // std::cout << "lbods: " << lbods << std::endl;
            // std::cout << "whbods: " << whbods << std::endl;
            std::vector<vec> first_row; // first row of bodies
            first_row.reserve(lbods);
            uint64_t counter = 0;
            while (counter < lbods) // create position vectors for all bodies in the first row
                first_row.emplace_back(counter++*b_sep, 0); // z-coordinate automatically zero
            std::vector<vec> second_row;
            second_row.reserve(lbods);
            T y_sep = sqrtl(3)*b_halfsep; // separation between adjacent rows, also the second row's y-coordinate
            counter = 0;
            while (counter < lbods)
                second_row.emplace_back(b_halfsep + counter++*b_sep, y_sep);
            std::vector<std::vector<vec>> plane_A;
            plane_A.reserve(whbods);
            plane_A.push_back(first_row);
            plane_A.push_back(second_row);
            std::vector<vec> row;
            T y_coord{};
            counter = 1;
            while (++counter < whbods) {
                row.reserve(lbods);
                y_coord = counter*y_sep;
                if (!(counter % 2))
                    for (const vec &v : first_row)
                        row.emplace_back(v.x, y_coord);
                else
                    for (const vec &v : second_row)
                        row.emplace_back(v.x, y_coord);
                plane_A.emplace_back(std::move(row));
                row.clear();
            }
            y_coord = b_halfsep / sqrtl(3);
            T y_coord2 = y_coord + y_sep;
            T z_sep = b_sep*sqrtl(2.0l/3.0l);
            vec *ptr1 = first_row.data();
            vec *ptr2 = second_row.data();
            counter = 0;
            while (counter++ < lbods) {
                ptr1->x += b_halfsep;
                ptr1->y = y_coord;
                ptr1++->z = z_sep;
                ptr2->x += b_halfsep;
                ptr2->y = y_coord2;
                ptr2++->z = z_sep;
            }
            std::vector<std::vector<vec>> plane_B;
            plane_B.emplace_back(first_row);
            plane_B.emplace_back(second_row);
            counter = 1;
            while (++counter < whbods) {
                row.reserve(lbods);
                y_coord2 = y_coord + counter*y_sep;
                if (!(counter % 2))
                    for (const vec &v : first_row)
                        row.emplace_back(v.x, y_coord2, z_sep);
                else
                    for (const vec &v : second_row)
                        row.emplace_back(v.x, y_coord2, z_sep);
                plane_B.emplace_back(std::move(row));
                row.clear();
            }
            std::vector<std::vector<std::vector<vec>>> cube = {plane_A, plane_B};
            std::vector<std::vector<vec>> plane;
            counter = 1;
            while (++counter < whbods) {
                plane.reserve(whbods);
                y_coord2 = counter*z_sep; // reusing y_coord2 - actually represents z-coordinate in this loop
                if (!(counter % 2))
                    for (const std::vector<vec> &_row : plane_A) {
                        row.reserve(lbods);
                        for (const vec &v : _row) {
                            row.emplace_back(v.x, v.y, y_coord2);
                        }
                        plane.emplace_back(std::move(row));
                        row.clear();
                    }
                else
                    for (const std::vector<vec> &_row : plane_B) {
                        row.reserve(lbods);
                        for (const vec &v : _row) {
                            row.emplace_back(v.x, v.y, y_coord2);
                        }
                        plane.emplace_back(std::move(row));
                        row.clear();
                    }
                cube.emplace_back(std::move(plane));
                plane.clear();
            }
            /* By this point, the positions for the entire HCP-cube will have been generated. Now the central position
             * will be selected, and the points that fall within the comet's radius of the central point will be used to
             * construct gtd::body<> objects. */
            vec centre;
            if (force_central_body) { // comet will be centred on a particular body
                centre = cube[whbods/2][whbods/2][lbods/2];
            }
            else { // comet will be centred on the geometrical centre of the HCP cube generated
                centre.x = cube[0][1][lbods - 1].x/2.0l;// - cube[0][0][0].x;
                centre.y = cube[0][whbods - 1][0].y/2.0l;// - cube[0][0][0].y;
                centre.z = cube[whbods - 1][0][0].z/2.0l;// - cube[0][0][0].z;
            }
            std::vector<bod_t> bodies;
            uint64_t half_lbods = lbods/2;
            vec *_beg;
            vec *_end;
            vec *_h1;
            vec *_h2;
            /* Now comes the juicy part. In this loop, every row of every plane in the cube is iterated over. The
             * pointers _h1 and _h2 ( == half 1 and half 2) are made to point to the central two vectors within a single
             * row. Whilst the vectors fall within the radius of the comet (by calculating their distances to its
             * centre) _h1 and _h2 are decremented and incremented, respectively, whilst constructing bodies at those
             * points. Once _h1 or _h2 point to a vector that falls outside the radius, the iteration stops, so as not
             * to compute distances for the following vectors pointlessly, and the next row is moved on to. */
            /* If lbods is even, lbods/2 iterations will be performed for _h1 and _h2. If lbods is odd, lbods/2.0 - 0.5
             * iterations will be performed for _h1 and lbods/2.0 + 0.5 iterations for _h2 (1 more for _h2). */
            counter = 0;
            typename std::vector<bod_t>::size_type n_bods;
            // discard bodies outside sphere and displace comet:
            if ((!orientation || orientation == vec::up) && !_omega && !adjust_to_com) {
                const vec offset = pos - centre;
                for (std::vector<std::vector<vec>> &_plane: cube) {
                    for (std::vector<vec> &_row: _plane) {
                        _beg = _row.data(); // first vector3D in row
                        _end = _beg + lbods; // past the last vector3D in row
                        _h1 = _beg + half_lbods;
                        _h2 = _h1--;
                        if constexpr (strict_rad) {
                            while (_h1 >= _beg && (vec_ops::distance(*_h1, centre) + b_rad) <= rad)
                                bodies.emplace_back(b_mass, b_rad, *_h1-- + offset, vel, r_coeff);
                            while (_h2 < _end && (vec_ops::distance(*_h2, centre) + b_rad) <= rad)
                                bodies.emplace_back(b_mass, b_rad, *_h2++ + offset, vel, r_coeff);
                        } else {
                            while (_h1 >= _beg && vec_ops::distance(*_h1, centre) <= rad) {
                                bodies.emplace_back(b_mass, b_rad, *_h1-- + offset, vel, r_coeff);
                            }
                            while (_h2 < _end && vec_ops::distance(*_h2, centre) <= rad) {
                                bodies.emplace_back(b_mass, b_rad, *_h2++ + offset, vel, r_coeff);
                            }
                        }
                    }
                }
                n_bods = bodies.size();
                return {std::piecewise_construct,
                        std::forward_as_tuple(std::move(bodies), timestep, num_iterations, units_format),
                        std::forward_as_tuple(cbrtl((SPHERE_VOLUME(b_rad)*n_bods*9*_ROOT_2_)/(4*_PI_*_PI_))*
                                              (1 + spacing/(2*b_rad)))};
            }
            for (std::vector<std::vector<vec>> &_plane: cube) { // discard bodies but do not displace comet
                for (std::vector<vec> &_row: _plane) {
                    _beg = _row.data();
                    _end = _beg + lbods;
                    _h1 = _beg + half_lbods;
                    _h2 = _h1--;
                    if constexpr (strict_rad) {
                        while (_h1 >= _beg && (vec_ops::distance(*_h1, centre) + b_rad) <= rad)
                            bodies.emplace_back(b_mass, b_rad, *_h1--, vel, r_coeff);
                        while (_h2 < _end && (vec_ops::distance(*_h2, centre) + b_rad) <= rad)
                            bodies.emplace_back(b_mass, b_rad, *_h2++, vel, r_coeff);
                    } else {
                        while (_h1 >= _beg && vec_ops::distance(*_h1, centre) <= rad)
                            bodies.emplace_back(b_mass, b_rad, *_h1--, vel, r_coeff);
                        while (_h2 < _end && vec_ops::distance(*_h2, centre) <= rad)
                            bodies.emplace_back(b_mass, b_rad, *_h2++, vel, r_coeff);
                    }
                }
            }
            n_bods = bodies.size();
            if (adjust_to_com) { // centre will now be made the centre-of-mass (COM)
                centre.make_zero();
                for (const bod_t &_b : bodies)
                    centre += _b.curr_pos;
                centre /= n_bods; // for bodies of equal mass, COM is simply the sum of positions div. by num. bodies
            }
            rotate_bodies(bodies, centre, vec::up, orientation); // these functions will simply return if nothing is to
            spin_bodies(bodies, centre, _omega); // be performed
            /*
            if (_omega) {
                const vec omega_hat = _omega.unit_vector();
                const auto omega_mag = _omega.magnitude();
                vec r_r;
                vec tan_vel; // tangential velocity of body w.r.t centre of comet
                if (orientation && orientation != vec::up) {
                    const vec rot_vec = vec_ops::cross(vec::up, orientation);
                    const long double _theta = vec_ops::angle_between(vec::up, orientation);
                    for (std::vector<std::vector<vec>> &_plane: cube) {
                        for (std::vector<vec> &_row: _plane) {
                            _beg = _row.data(); // first vector3D in row
                            _end = _beg + lbods; // past the last vector3D in row
                            _h1 = _beg + half_lbods;
                            _h2 = _h1--;
                            if constexpr (strict_rad) {
                                while (_h1 >= _beg && ((r_r = *_h1-- - centre).magnitude() + b_rad) <= rad) {
                                    tan_vel = (r_r - (r_r*omega_hat)*omega_hat).magnitude();
                                    bodies.emplace_back(b_mass, b_rad, r_r.rodrigues_rotate(rot_vec, _theta) + pos, vel,
                                                        r_coeff);
                                }
                                while (_h2 < _end && ((r_r = *_h2++ - centre).magnitude() + b_rad) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad, r_r.rodrigues_rotate(rot_vec, _theta) + pos, vel,
                                                        r_coeff);
                                }
                            } else {
                                while (_h1 >= _beg && vec_ops::distance(*_h1, centre) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad,
                                                        (*_h1-- - centre).rodrigues_rotate(rot_vec, _theta) + pos, vel,
                                                        r_coeff);
                                }
                                while (_h2 < _end && vec_ops::distance(*_h2, centre) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad,
                                                        (*_h2++ - centre).rodrigues_rotate(rot_vec, _theta) + pos, vel,
                                                        r_coeff);
                                }
                            }
                        }
                    }
                } else {
                    for (std::vector<std::vector<vec>> &_plane: cube) {
                        for (std::vector<vec> &_row: _plane) {
                            _beg = _row.data(); // first vector3D in row
                            _end = _beg + lbods; // past the last vector3D in row
                            _h1 = _beg + half_lbods;
                            _h2 = _h1--;
                            if constexpr (strict_rad) {
                                while (_h1 >= _beg && (vec_ops::distance(*_h1, centre) + b_rad) <= rad)
                                    bodies.emplace_back(b_mass, b_rad, *_h1-- - centre + pos, vel, r_coeff);
                                while (_h2 < _end && (vec_ops::distance(*_h2, centre) + b_rad) <= rad)
                                    bodies.emplace_back(b_mass, b_rad, *_h2++ - centre + pos, vel, r_coeff);
                            } else {
                                while (_h1 >= _beg && vec_ops::distance(*_h1, centre) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad, *_h1-- - centre + pos, vel, r_coeff);
                                }
                                while (_h2 < _end && vec_ops::distance(*_h2, centre) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad, *_h2++ - centre + pos, vel, r_coeff);
                                }
                            }
                        }
                    }
                }
            } else {
                if (orientation && orientation != vec::up) {
                    const vec rot_vec = vec_ops::cross(vec::up, orientation);
                    const long double _theta = vec_ops::angle_between(vec::up, orientation);
                    for (std::vector<std::vector<vec>> &_plane: cube) {
                        for (std::vector<vec> &_row: _plane) {
                            _beg = _row.data(); // first vector3D in row
                            _end = _beg + lbods; // past the last vector3D in row
                            _h1 = _beg + half_lbods;
                            _h2 = _h1--;
                            if constexpr (strict_rad) {
                                while (_h1 >= _beg && (vec_ops::distance(*_h1, centre) + b_rad) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad,
                                                        (*_h1-- - centre).rodrigues_rotate(rot_vec, _theta) + pos, vel,
                                                        r_coeff);
                                }
                                while (_h2 < _end && (vec_ops::distance(*_h2, centre) + b_rad) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad,
                                                        (*_h2++ - centre).rodrigues_rotate(rot_vec, _theta) + pos, vel,
                                                        r_coeff);
                                }
                            } else {
                                while (_h1 >= _beg && vec_ops::distance(*_h1, centre) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad,
                                                        (*_h1-- - centre).rodrigues_rotate(rot_vec, _theta) + pos, vel,
                                                        r_coeff);
                                }
                                while (_h2 < _end && vec_ops::distance(*_h2, centre) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad,
                                                        (*_h2++ - centre).rodrigues_rotate(rot_vec, _theta) + pos, vel,
                                                        r_coeff);
                                }
                            }
                        }
                    }
                } else {
                    for (std::vector<std::vector<vec>> &_plane: cube) {
                        for (std::vector<vec> &_row: _plane) {
                            _beg = _row.data(); // first vector3D in row
                            _end = _beg + lbods; // past the last vector3D in row
                            _h1 = _beg + half_lbods;
                            _h2 = _h1--;
                            if constexpr (strict_rad) {
                                while (_h1 >= _beg && (vec_ops::distance(*_h1, centre) + b_rad) <= rad)
                                    bodies.emplace_back(b_mass, b_rad, *_h1-- - centre + pos, vel, r_coeff);
                                while (_h2 < _end && (vec_ops::distance(*_h2, centre) + b_rad) <= rad)
                                    bodies.emplace_back(b_mass, b_rad, *_h2++ - centre + pos, vel, r_coeff);
                            } else {
                                while (_h1 >= _beg && vec_ops::distance(*_h1, centre) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad, *_h1-- - centre + pos, vel, r_coeff);
                                }
                                while (_h2 < _end && vec_ops::distance(*_h2, centre) <= rad) {
                                    bodies.emplace_back(b_mass, b_rad, *_h2++ - centre + pos, vel, r_coeff);
                                }
                            }
                        }
                    }
                }
            }
            */
            centre = pos - centre; // displacement offset (re-using centre variable)
            for (bod_t &_b : bodies)
                _b.curr_pos += centre;
            /* As of C++17, copy elision GUARANTEES that a returned prvalue is constructed directly in the memory
             * location of the return value, thus avoiding unnecessary copies. */
            return {std::piecewise_construct,
                    std::forward_as_tuple(std::move(bodies), timestep, num_iterations, units_format),
                    std::forward_as_tuple(cbrtl((SPHERE_VOLUME(b_rad)*n_bods*9*_ROOT_2_)/(4*_PI_*_PI_))*
                    (1 + spacing/(2*b_rad)))};
        }
        static void rotate_bodies(std::vector<bod_t> &bodies, const vec_t &about, const vec_t &_i_ornt,
                                  const vec_t &_f_ornt) {
            if (bodies.empty() || !_i_ornt || !_f_ornt || _i_ornt == _f_ornt)
                return;
            const vec_t rot_vec = vec_ops::cross(_i_ornt, _f_ornt);
            const long double _theta = vec_ops::angle_between(_i_ornt, _f_ornt);
            for (bod_t &_b : bodies)
                _b.curr_pos = about + (_b.curr_pos - about).rodrigues_rotate(rot_vec, _theta);
        }
        static void spin_bodies(std::vector<bod_t> &bodies, const vec_t &about, const vec_t &_omega,bool assign=false) {
            if (bodies.empty() || !_omega)
                return;
            const auto omega_mag = _omega.magnitude(); // magnitude of the angular velocity
            const vec_t omega_hat = _omega / omega_mag; // ang.vel unit vec - more efficient than calling .unit_vector()
            vec_t rel_pos; // position of a given body relative to about
            if (assign) // discard previous velocity
                for (bod_t &_b : bodies) {
                    rel_pos = _b.curr_pos - about;
                    _b.curr_vel = vec_ops::cross(_omega, rel_pos - (rel_pos*omega_hat)*omega_hat);
                }
            else // add to previous velocity
                for (bod_t &_b : bodies) {
                    rel_pos = _b.curr_pos - about;
                    _b.curr_vel += vec_ops::cross(_omega, rel_pos - (rel_pos*omega_hat)*omega_hat);
                }
        }
        void reset() requires (memFreq != 0) {
            // for (bod_t &bod : bods)
            //     bod.reset(true);
#ifdef GREGSYS_MERGERS
            if constexpr (collisions) {
                bods.insert(bods.end(), std::make_move_iterator(del_bods_m->begin()),
                            std::make_move_iterator(del_bods_m->end()));
                del_bods_m->clear();
            }
#endif
            std::for_each(bods.begin(), bods.end(), [this](bod_t &bod){bod.reset(true);});
            pe.clear();
            energy.clear();
            tot_pe.clear();
            tot_ke.clear();
            tot_e.clear();
            evolved = false;
            time_elapsed = 0;
        }
        /* reset() deletes the evolution of the bodies AND returns them to their original positions, velocities, and
         * kinetic energies, whilst clear_evolution() deletes the evolution BUT retains the current positions,
         * velocities and kinetic energies of the bodies. */
        void clear_evolution() requires (memFreq != 0) {
            // for (bod_t &bod : bods)
            //     bod.clear();
            // bods.insert(bods.end(), std::make_move_iterator(del_bods_m->begin()),
            //             std::make_move_iterator(del_bods_m->end()));
#ifdef GREGSYS_MERGERS
            if constexpr (collisions)
                del_bods_m->clear();
#endif
            std::for_each(bods.begin(), bods.end(), [](bod_t &bod){bod.clear();});
            pe.clear();
            energy.clear();
            tot_pe.clear();
            tot_ke.clear();
            tot_e.clear();
            evolved = false;
            time_elapsed = 0;
        }
    private:
        const char *&method_str() requires (prog) {
            static const char *ptr{};
            return ptr;
        }
    public:
        template <typename evFuncT = std::nullptr_t,
                  typename bodFuncT = std::nullptr_t,
                  typename bodsFuncT = std::nullptr_t>
        requires (std::same_as<evFuncT, std::nullptr_t> ||
        requires (evFuncT f1) {
            {f1()} -> std::same_as<bool>;
        }) && (std::same_as<bodFuncT, std::nullptr_t> || requires (bodFuncT f2, bod_t &bod) {
            {f2(bod)} -> std::same_as<void>;
        }) && (std::same_as<bodsFuncT, std::nullptr_t> || requires (bodsFuncT f2, std::vector<bod_t> &_bodies) {
            {f2(_bodies)} -> std::same_as<void>;
        })
        const std::chrono::nanoseconds evolve(int integration_method = leapfrog_kdk,
                                               const evFuncT &continue_ev_if = nullptr,
                                               const bodFuncT &for_each_bod = nullptr,
                                               const bodsFuncT &for_all_bods = nullptr) {
            std::chrono::time_point<std::chrono::high_resolution_clock> start; // removed static
            std::chrono::nanoseconds total; // removed static
            num_bods = bods.size();
            if (!num_bods || !check_option(integration_method))
                return total = std::chrono::nanoseconds::zero();
            // time_t start = time(nullptr);
            if constexpr (std::same_as<evFuncT, std::nullptr_t> || memFreq || fileFreq || prog)
                steps = 0;
            if constexpr (memFreq) {
                if (evolved)
                    this->clear_evolution();
                this->create_energy_vectors();
            }
#ifndef _WIN32
            if constexpr (prog)
                std::cout << BOLD_TXT_START;
#endif
            start = std::chrono::high_resolution_clock::now();
            if ((integration_method & two_body) == two_body) {
                if (num_bods != 2)
                    throw two_body_error();
                /* no force approximation techniques performed here, as the integration is just for two bodies */
            }
            else if ((integration_method & euler) == euler) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                    typename cube_t::iterator it;
                    typename cube_t::cell_iterator cit{bh_theta};
                    while (steps < iterations) {
                        btree.add_bodies(this->bods);
                        it = btree.begin();
                        cit.set_root(btree._root);
                        if constexpr (memFreq && fileFreq) {
                            if (!(steps % memFreq) || !(steps % fileFreq)) {
                                while (it) {
                                    cit = it;
                                    while (cit) {
                                        this->bh_acc_and_pe(it->_cptr->_bod, cit++->_cptr);
                                    }
                                    ++it;
                                }
                            }
                            else {

                            }
                        }
                        else if constexpr (memFreq) {

                        }
                        else if constexpr (fileFreq) {

                        }
                        else {

                        }
                        btree.delete_tree(); // has to be reconstructed at every step since particles positions change
                    }
                }
                else {
                    if constexpr (std::same_as<evFuncT, std::nullptr_t>) {
                        while (steps < iterations) {
                            if constexpr (collisions == pred_coll_check) {

                            }
                            this->calc_acc_and_e();
                            ++steps;
                            // this->take_euler_step();
                            if constexpr (!std::same_as<bodsFuncT, std::nullptr_t>) {
                                FUNC_TEMPL_SELECT(this->take_euler_step, for_all_bods(this->bods);, EMPTY)
                            } else {
                                FUNC_TEMPL_SELECT(this->take_euler_step, EMPTY, EMPTY)
                            }
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            /* by having "prog" as a template parameter, it avoids having to re-evaluate whether it is true
                             * or not within the loop, since the "if constexpr ()" branches get rejected at compile-time */
                            // if constexpr (collisions == overlap_coll_check)
                            //     this->s_coll();
                            if constexpr (prog)
                                this->print_progress();
                        }
                    } else {
                        while (continue_ev_if()) {
                            if constexpr (collisions == pred_coll_check) {

                            }
                            this->calc_acc_and_e();
                            if constexpr (memFreq || fileFreq || prog)
                                ++steps;
                            if constexpr (!std::same_as<bodsFuncT, std::nullptr_t>) {
                                FUNC_TEMPL_SELECT(this->take_euler_step, for_all_bods(this->bods);, EMPTY)
                            } else {
                                FUNC_TEMPL_SELECT(this->take_euler_step, EMPTY, EMPTY)
                            }
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress<false>();
                        }
                    }
                    if constexpr (prog)
                        this->method_str() = "Euler";
                }
            }
            else if ((integration_method & modified_euler) == modified_euler) {
                /* The below std::vector will hold the euler-predicted values for position, velocity and acceleration,
                 * respectively. These are subsequently used to determine the corrected values. */
                std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>> predicted{num_bods};
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    if constexpr (std::same_as<evFuncT, std::nullptr_t>) {
                        while (steps < iterations) {
                            this->calc_acc_and_e();
                            for (outer = 0; outer < num_bods; ++outer) {
                                bod_t &ref = bods[outer];
                                std::get<0>(predicted[outer]) = ref.curr_pos + dt*ref.curr_vel;
                                std::get<1>(predicted[outer]) = ref.curr_vel + dt*ref.acc;
                            }
                            for (outer = 0; outer < num_bods; ++outer) {
                                bod_t &ref = bods[outer];
                                std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &pred_outer = predicted[outer];
                                for (inner = outer + 1; inner < num_bods; ++inner) {
                                    bod_t &ref_i = bods[inner];
                                    vector3D<T> &&r12 = std::get<0>(predicted[inner]) - std::get<0>(pred_outer);
                                    long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
                                    std::get<2>(pred_outer) += ((G*ref_i.mass_)/(r12_cubed_mag))*r12;
                                    std::get<2>(predicted[inner]) -= ((G*ref.mass_)/(r12_cubed_mag))*r12;
                                }
                            }
                            // this->take_modified_euler_step(predicted);
                            ++steps;
                            if constexpr (!std::same_as<bodsFuncT, std::nullptr_t>) {
                                FUNC_TEMPL_SELECT(this->take_modified_euler_step, for_all_bods(this->bods);, EMPTY,
                                                  predicted)
                            } else {
                                FUNC_TEMPL_SELECT(this->take_modified_euler_step, EMPTY, EMPTY, predicted)
                            }
                            for (std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &tup: predicted)
                                std::get<2>(tup).make_zero();
                            if constexpr (collisions == overlap_coll_check) {
                                // this->s_coll();
                                num_bods = bods.size(); // number of bodies can change through mergers
                            }
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress();
                        }
                    } else {
                        while (continue_ev_if()) {
                            this->calc_acc_and_e();
                            for (outer = 0; outer < num_bods; ++outer) {
                                bod_t &ref = bods[outer];
                                std::get<0>(predicted[outer]) = ref.curr_pos + dt*ref.curr_vel;
                                std::get<1>(predicted[outer]) = ref.curr_vel + dt*ref.acc;
                            }
                            for (outer = 0; outer < num_bods; ++outer) {
                                bod_t &ref = bods[outer];
                                std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &pred_outer = predicted[outer];
                                for (inner = outer + 1; inner < num_bods; ++inner) {
                                    bod_t &ref_i = bods[inner];
                                    vector3D<T> &&r12 = std::get<0>(predicted[inner]) - std::get<0>(pred_outer);
                                    long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
                                    std::get<2>(pred_outer) += ((G*ref_i.mass_)/(r12_cubed_mag))*r12;
                                    std::get<2>(predicted[inner]) -= ((G*ref.mass_)/(r12_cubed_mag))*r12;
                                }
                            }
                            if constexpr (memFreq || fileFreq || prog)
                                ++steps;
                            if constexpr (!std::same_as<bodsFuncT, std::nullptr_t>) {
                                FUNC_TEMPL_SELECT(this->take_modified_euler_step, for_all_bods(this->bods);, EMPTY,
                                                  predicted)
                            } else {
                                FUNC_TEMPL_SELECT(this->take_modified_euler_step, EMPTY, EMPTY, predicted)
                            }
                            for (std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &tup: predicted)
                                std::get<2>(tup).make_zero();
                            if constexpr (collisions == overlap_coll_check) {
                                num_bods = bods.size();
                            }
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress<false>();
                        }
                    }
                    if constexpr (prog)
                        this->method_str() = "Modified Euler";
                }
            }
            else if ((integration_method & midpoint) == midpoint) {
                std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>> predicted{num_bods};
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    if constexpr (std::same_as<evFuncT, std::nullptr_t>) {
                        while (steps < iterations) {
                            this->calc_acc_and_e();
                            for (outer = 0; outer < num_bods; ++outer) {
                                bod_t &ref = bods[outer];
                                std::get<0>(predicted[outer]) = ref.curr_pos + half_dt*ref.curr_vel;
                                std::get<1>(predicted[outer]) = ref.curr_vel + half_dt*ref.acc;
                            }
                            for (outer = 0; outer < num_bods; ++outer) {
                                bod_t &ref = bods[outer];
                                std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &pred_outer = predicted[outer];
                                for (inner = outer + 1; inner < num_bods; ++inner) {
                                    bod_t &ref_i = bods[inner];
                                    vector3D<T> &&r12 = std::get<0>(predicted[inner]) - std::get<0>(pred_outer);
                                    long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
                                    std::get<2>(pred_outer) += ((G*ref_i.mass_)/(r12_cubed_mag))*r12;
                                    std::get<2>(predicted[inner]) -= ((G*ref.mass_)/(r12_cubed_mag))*r12;
                                }
                            }
                            // this->take_midpoint_step(predicted);
                            ++steps;
                            if constexpr (!std::same_as<bodsFuncT, std::nullptr_t>) {
                                FUNC_TEMPL_SELECT(this->take_midpoint_step, for_all_bods(this->bods);, EMPTY, predicted)
                            } else {
                                FUNC_TEMPL_SELECT(this->take_midpoint_step, EMPTY, EMPTY, predicted)
                            }
                            for (std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &tup: predicted)
                                std::get<2>(tup).make_zero();
                            if constexpr (collisions == overlap_coll_check) {
                                // this->s_coll();
                                num_bods = bods.size();
                            }
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress();
                        }
                    } else {
                        while (continue_ev_if()) {
                            this->calc_acc_and_e();
                            for (outer = 0; outer < num_bods; ++outer) {
                                bod_t &ref = bods[outer];
                                std::get<0>(predicted[outer]) = ref.curr_pos + half_dt*ref.curr_vel;
                                std::get<1>(predicted[outer]) = ref.curr_vel + half_dt*ref.acc;
                            }
                            for (outer = 0; outer < num_bods; ++outer) {
                                bod_t &ref = bods[outer];
                                std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &pred_outer = predicted[outer];
                                for (inner = outer + 1; inner < num_bods; ++inner) {
                                    bod_t &ref_i = bods[inner];
                                    vector3D<T> &&r12 = std::get<0>(predicted[inner]) - std::get<0>(pred_outer);
                                    long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
                                    std::get<2>(pred_outer) += ((G*ref_i.mass_)/(r12_cubed_mag))*r12;
                                    std::get<2>(predicted[inner]) -= ((G*ref.mass_)/(r12_cubed_mag))*r12;
                                }
                            }
                            // this->take_midpoint_step(predicted);
                            if constexpr (memFreq || fileFreq || prog)
                                ++steps;
                            if constexpr (!std::same_as<bodsFuncT, std::nullptr_t>) {
                                FUNC_TEMPL_SELECT(this->take_midpoint_step, for_all_bods(this->bods);, EMPTY, predicted)
                            } else {
                                FUNC_TEMPL_SELECT(this->take_midpoint_step, EMPTY, EMPTY, predicted)
                            }
                            for (std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &tup: predicted)
                                std::get<2>(tup).make_zero();
                            if constexpr (collisions == overlap_coll_check) {
                                // this->s_coll();
                                num_bods = bods.size();
                            }
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress<false>();
                        }
                    }
                    if constexpr (prog)
                        this->method_str() = "Midpoint method";
                }
            }
            else if ((integration_method & leapfrog_kdk) == leapfrog_kdk) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                    btree.delete_tree();
                    btree.add_bodies(this->bods);
                    typename cube_t::iterator it = btree.begin();
                    // typename cube_t::cell_iterator cit{btree._root, this->bh_theta};
                    // cit.set_root(btree._root);
                    while (it) { // need an initial acc.
                        // cit = it;
                        typename cube_t::cell_iterator cit{btree._root, it, this->bh_theta};
                        while (cit) {
                            // std::cout << "COM: " << cit._cptr->_com << std::endl;
                            this->bh_acc(const_cast<bod_t*>(it._cptr->_bod), cit++._cptr);
                        }
                        // std::cout << *it << std::endl;
                        ++it;
                    }
                    for (bod_t &bod : bods) {
                        /* KICK for half a step */
                        bod.curr_vel += half_dt*bod.acc;
                    }
                    while (steps++ < iterations) {
                        btree.delete_tree();
                        for (bod_t &bod : bods) {
                            /* DRIFT for a full step */
                            bod.curr_pos += dt*bod.curr_vel;
                        }
                        btree.add_bodies(this->bods);
                        it = btree.begin();
                        // cit.set_root(btree._root);
                        while (it) {
                            // cit = it;
                            typename cube_t::cell_iterator cit{btree._root, it, this->bh_theta};
                            while (cit) {
                                // std::cout << "COM: " << cit._cptr->_com << std::endl;
                                this->bh_acc(const_cast<bod_t*>(it._cptr->_bod), cit++._cptr);
                            }
                            // std::cout << *it << std::endl;
                            ++it;
                        }
                        for (bod_t &bod : bods) {
                            /* KICK for full step */
                            bod.curr_vel += dt*bod.acc;
                        }
                        if constexpr (collisions) {
                            this->bh_s_coll();
                        }
                        if constexpr (prog)
                            this->print_progress();
                    }
                    if constexpr (prog)
                        this->method_str() = "Leapfrog KDK and Barnes-Hut";
                }
                else {
                    this->calc_acc_and_e(); // need an initial acceleration
                    if constexpr (std::same_as<evFuncT, std::nullptr_t>) {
                        while (steps++ < iterations) {
                            for (bod_t &bod: bods) {
                                /* KICK for half a step */
                                bod.curr_vel += half_dt*bod.acc;
                                /* DRIFT for a full step */
                                bod.curr_pos += dt*bod.curr_vel;
                            }
                            /* update accelerations and energies based on new positions and KICK for half a step */
                            if constexpr (std::same_as<bodsFuncT, std::nullptr_t>)
                                this->leapfrog_kdk_acc_e_and_step();
                            else
                                this->leapfrog_kdk_acc_e_and_step(for_all_bods);
                            // if constexpr (collisions == overlap_coll_check)
                            //     this->s_coll();
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress();
                        }
                    } else {
                        while (continue_ev_if()) {
                            if constexpr (memFreq || fileFreq || prog)
                                ++steps;
                            for (bod_t &bod: bods) {
                                bod.curr_vel += half_dt*bod.acc;
                                bod.curr_pos += dt*bod.curr_vel;
                            }
                            if constexpr (std::same_as<bodsFuncT, std::nullptr_t>)
                                this->leapfrog_kdk_acc_e_and_step();
                            else
                                this->leapfrog_kdk_acc_e_and_step(for_all_bods);
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress<false>();
                        }
                    }
                    if constexpr (prog)
                        this->method_str() = "Leapfrog KDK";
                }
                goto bye;
            }
            else if ((integration_method & leapfrog_dkd) == leapfrog_dkd) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    if constexpr (memFreq)
                        this->calc_energy();
                    if constexpr (std::same_as<evFuncT, std::nullptr_t>) {
                        while (steps++ < iterations) {
                            for (bod_t &bod: bods)
                                /* DRIFT for half a step */
                                bod.curr_pos += half_dt*bod.curr_vel;
                            /* update acc., then KICK for a full step and DRIFT for half a step, then update E */
                            if constexpr (std::same_as<bodsFuncT, std::nullptr_t>)
                                this->leapfrog_dkd_acc_e_and_step();
                            else
                                this->leapfrog_dkd_acc_e_and_step(for_all_bods);
                            // if constexpr (collisions == overlap_coll_check)
                            //     this->s_coll();
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress();
                        }
                    } else {
                        while (continue_ev_if()) {
                            if constexpr (memFreq || fileFreq || prog)
                                ++steps;
                            for (bod_t &bod: bods)
                                bod.curr_pos += half_dt*bod.curr_vel;
                            if constexpr (std::same_as<bodsFuncT, std::nullptr_t>)
                                this->leapfrog_dkd_acc_e_and_step();
                            else
                                this->leapfrog_dkd_acc_e_and_step(for_all_bods);
                            if constexpr (!std::same_as<bodFuncT, std::nullptr_t>)
                                for (bod_t &_b : this->bods)
                                    for_each_bod(_b);
                            if constexpr (prog)
                                this->print_progress<false>();
                        }
                    }
                    if constexpr (prog)
                        this->method_str() = "Leapfrog DKD";
                }
                goto bye;
            }
            else if ((integration_method & rk4) == rk4) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {

                }
            }
            else if ((integration_method & rk3_8) == rk3_8) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {

                }
            }
            if constexpr (memFreq) {
                if (!(steps % memFreq)) {
                    this->calc_energy();
                }
            }
            bye:
            if constexpr (prog)
                print_conclusion(start, total);
            else
                total = std::chrono::high_resolution_clock::now() - start;
            evolved = true;
            prev_dt = dt;
            prev_iterations = iterations;
            time_elapsed = iterations*dt;
            return total;
        }
        std::pair<T*, long double> *pe_extrema() const requires (memFreq != 0) {
            /* This method returns a pointer to a 2-element array of std::pairs. The first pair concerns the minima and
             * the second pair concerns the maxima. The first element of the pairs is a pointer to a 2-element array of
             * type T. The first element of the array is the actual extreme value, and the second is the difference
             * between the extreme value and starting PE. The second element of the pairs is the time at which this
             * extreme occurs. */
            if (bods.empty())
                throw empty_system_error("pe_diff() cannot be called on empty system objects (no bodies present).\n");
            if (!evolved)
                throw no_evolution_error("pe_diff() can only be called after a system object has been evolved.\n");
            T minima[2]; // first the min. PE, then the difference between it and the starting PE | removed static
            T maxima[2]; // same but for max. | removed static
            std::pair<T*, long double> min_max[2] = {{minima, 0}, {maxima, 0}}; // removed static
            unsigned long long min_index{};
            unsigned long long max_index{};
            unsigned long long count{};
            const T &first_element = tot_pe[0]; // avoids 3 function calls
            T min = first_element;
            T max = first_element;
            for (const T &pot_e : tot_pe) {
                if (pot_e < min) {
                    min = pot_e;
                    min_index = count;
                }
                if (pot_e > max) {
                    max = pot_e;
                    max_index = count;
                }
                ++count;
            }
            *minima = min; // assign extrema
            *maxima = max;
            *(minima + 1) = first_element - min; // assign differences from initial PE
            *(maxima + 1) = max - first_element;
            min_max->second = min_index*prev_dt*memFreq; // assign times at which extrema occurred
            (min_max + 1)->second = max_index*prev_dt*memFreq;
            return min_max;
        }
        std::pair<T*, long double> *ke_extrema() const requires (memFreq != 0) {
            /* Same as the above method, but for KE. */
            if (bods.empty())
                throw empty_system_error("pe_diff() cannot be called on empty system objects (no bodies present).\n");
            if (!evolved)
                throw no_evolution_error("pe_diff() can only be called after a system object has been evolved.\n");
            T minima[2]; // first the min. PE, then the difference between it and the starting PE | removed static
            T maxima[2]; // same but for max. | removed static
            std::pair<T*, long double> min_max[2] = {{minima, 0}, {maxima, 0}}; // removed static
            unsigned long long min_index{};
            unsigned long long max_index{};
            unsigned long long count{};
            const T &first_element = tot_ke[0];
            T min = first_element;
            T max = first_element;
            for (const T &kin_e : tot_ke) {
                if (kin_e < min) {
                    min = kin_e;
                    min_index = count;
                }
                if (kin_e > max) {
                    max = kin_e;
                    max_index = count;
                }
                ++count;
            }
            *minima = min;
            *maxima = max;
            *(minima + 1) = first_element - min;
            *(maxima + 1) = max - first_element;
            min_max->second = min_index*prev_dt*memFreq;
            (min_max + 1)->second = max_index*prev_dt*memFreq;
            return min_max;
        }
        std::pair<T*, long double> *tot_e_extrema() const requires (memFreq != 0) {
            /* Same as the above method, but for KE. */
            if (bods.empty())
                throw empty_system_error("pe_diff() cannot be called on empty system objects (no bodies present).\n");
            if (!evolved)
                throw no_evolution_error("pe_diff() can only be called after a system object has been evolved.\n");
            T minima[2]; // first the min. PE, then the difference between it and the starting PE | removed static
            T maxima[2]; // same but for max. | removed static
            std::pair<T*, long double> min_max[2] = {{minima, 0}, {maxima, 0}}; // removed static
            unsigned long long min_index{};
            unsigned long long max_index{};
            unsigned long long count{};
            const T &first_element = tot_e[0];
            T min = first_element;
            T max = first_element;
            for (const T &_e : tot_e) {
                if (_e < min) {
                    min = _e;
                    min_index = count;
                }
                if (_e > max) {
                    max = _e;
                    max_index = count;
                }
                ++count;
            }
            *minima = min;
            *maxima = max;
            *(minima + 1) = first_element - min;
            *(maxima + 1) = max - first_element;
            min_max->second = min_index*prev_dt*memFreq;
            (min_max + 1)->second = max_index*prev_dt*memFreq;
            return min_max;
        }
        auto begin() {
            return bods.begin();
        }
        auto end() {
            return bods.end();
        }
        auto begin() const {
            return bods.cbegin();
        }
        auto end() const {
            return bods.cend();
        }
        auto cbegin() const {
            return bods.cbegin();
        }
        auto cend() const {
            return bods.cend();
        }
        auto rbegin() {
            return bods.rbegin();
        }
        auto rend() {
            return bods.rend();
        }
        auto rbegin() const {
            return bods.crbegin();
        }
        auto rend() const {
            return bods.crend();
        }
        auto crbegin() const {
            return bods.crbegin();
        }
        auto crend() const {
            return bods.crend();
        }
        std::streamoff to_nsys(const char *path) const requires (std::is_fundamental_v<M> &&
                                                                 std::is_fundamental_v<R> &&
                                                                 std::is_fundamental_v<T> && CHAR_BIT == 8) {
            nsys_header header; // removed static
            nsys_chunk chunk; // removed static
            std::ofstream::pos_type tell; // removed static
            if (path == nullptr)
                return 0;
            ostream = new std::ofstream{path, std::ios_base::trunc | std::ios_base::binary};
            if (!*ostream) {
                delete ostream;
                return 0;
            }
            copy(header.time_point, &this->time_elapsed, sizeof(long double));
            copy(header.delta_t, &this->dt, sizeof(long double));
            copy(&header.m_num, this->G_vals, 6*sizeof(uint64_t)); // I include sizeof just so intent is clear
            header.num_bodies = this->bods.size();
            ostream->write((char *) &header, sizeof(nsys_header));
            if (!header.num_bodies) { // case for no bodies - no reason to prohibit this
                ostream->close();
                delete ostream;
                return sizeof(nsys_header);
            }
            for (const bod_t &bod : bods) {
                // chunk.r_coeff = bod.rest_c;
                copy(chunk.r_coeff, &bod.rest_c, sizeof(long double));
                chunk.mass = bod.mass_;
                chunk.radius = bod.radius;
                chunk.xpos = bod.curr_pos.x;
                chunk.ypos = bod.curr_pos.y;
                chunk.zpos = bod.curr_pos.z;
                chunk.xvel = bod.curr_vel.x;
                chunk.yvel = bod.curr_vel.y;
                chunk.zvel = bod.curr_vel.z;
                ostream->write((char *) &chunk, sizeof(nsys_chunk));
            }
            tell = ostream->tellp();
            if (char rem = (char) (tell % 4); rem) {
                ostream->write("\0\0\0", 4 - rem); // align EOF with dword boundary - at most 3 bytes would be written
                tell += (std::streamoff) (4 - rem);
            }
            ostream->close();
            delete ostream;
            return tell;
        }
        // static std::pair<nsys_header, std::vector<bod_t>> load_nsys(const char *path, bool only_header = false) {
        //     if (path == nullptr)
        //         return {};
        // }
        bool write_trajectories(const String &path = def_path, bool binary = false) requires (memFreq != 0) {
            /* This method can only be called if memFreq is non-zero, or else there would be no history to output. */
            /* Writes the trajectories (historically) of all the bodies to a file, including all the energy values. For
             * text files, units are not written as these will depend entirely on the units_format string passed to the
             * system<M, R, T> object in its constructor. */
            if (!evolved || !bods.size())
                return false;
            if (binary) {
                if (&path == &def_path)
                    def_path.append_back(get_date_and_time()).append_back(".nbod");
                if (&path == &def_path)
                    def_path.erase_chars(def_path.get_length() - 30);
                return true;
            } // needs work - does NOT output energies
            if (&path == &def_path)
                def_path.append_back(get_date_and_time()).append_back(".csv");
            std::ofstream out(path.c_str(), std::ios_base::trunc);
            if (!out)
                return false;
            // unsigned long long count = 0;
            vec_size_t i;
            const ull_t num_max_els = (prev_iterations + (memFreq != 1))/memFreq; // CHECK THIS IS CORRECT
            // std::cout << "size of bods: " << bods.size() << ", size of del_bods: " << del_bods_m->size() << std::endl;
#ifdef GREGSYS_MERGERS
            if constexpr (collisions) {
                ull_t d_i;
                vec_size_t max_index;
                for (const auto &[destruction_it, bod] : *del_bods_m) {
                    max_index = bod.positions.size();
                    // d_i = (destruction_it + memFreq)/memFreq - max_index - 1;
                    if (destruction_it % memFreq) {
                        d_i = destruction_it/memFreq + 1 - max_index;
                    } else {
                        d_i = destruction_it/memFreq - max_index;
                    }
                    // d_i = destruction_it/memFreq;
                    // d_i += (max_index > d_i) - max_index;
                    // std::cout << "Positions size: " << bod.positions.size() << ", pe size: " << pe[bod.id].size()
                    // << std::endl;
                    std::cout << "body id: " << bod.id << std::endl << std::endl;
                    out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                           "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,"
                           "kinetic_energy,potential_energy,total_energy\r\n";
                    for (i = 0; i < max_index; ++i)
                        out << d_i++*prev_dt*memFreq << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                            << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                            << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                            << energy[bod.id][i] << "\r\n";
                    out << ",,,,,,,,,\r\n";
                }
                for (const bod_t &bod : bods) {
                    max_index = bod.positions.size();
                    d_i = prev_iterations/memFreq - max_index + 1;
                    // std::cout << "max_index: " << max_index << ", d_i: " << d_i << std::endl;
                    // std::cout << "body id: " << bod.id << std::endl;
                    out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                           "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,"
                           "kinetic_energy,potential_energy,total_energy\r\n";
                    for (i = 0; i < max_index; ++i)
                        out << d_i++*prev_dt*memFreq << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                            << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                            << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                            << energy[bod.id][i] << "\r\n";
                    out << ",,,,,,,,,\r\n";
                }
            }
            else {
#endif
                for (const bod_t &bod : bods) { // simple case of all bodies lasting from beginning to end
                    out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                           "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,"
                           "kinetic_energy,potential_energy,total_energy\r\n";
                    for (i = 0; i <= num_max_els; ++i)
                        out << i*prev_dt*memFreq << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                            << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                            << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                            << energy[bod.id][i] << "\r\n";
                    out << ",,,,,,,,,\r\n";
                }
#ifdef GREGSYS_MERGERS
            }
#endif
            out << "time_elapsed,system_PE,system_KE,system_E\r\n";
            for (i = 0; i <= num_max_els; ++i)
                out << i*prev_dt*memFreq << ',' << tot_pe[i] << ',' << tot_ke[i] << ',' << tot_e[i] << "\r\n";
            out << ",,,,,,,,,\r\n";
            out << "energy_type,min,max,min_abs_diff_from_start_val,"
                   "max_abs_diff_from_start_val,time_of_min,time_of_max\r\n";
            std::pair<T*, long double> *ptr = this->pe_extrema();
            out << "system_PE," << *ptr->first << ',' << *((ptr + 1)->first) << ',' << *(ptr->first + 1) << ','
                << *((ptr + 1)->first + 1) << ',' << ptr->second << ',' << (ptr + 1)->second << "\r\n";
            ptr = this->ke_extrema();
            out << "system_KE," << *ptr->first << ',' << *((ptr + 1)->first) << ',' << *(ptr->first + 1) << ','
                << *((ptr + 1)->first + 1) << ',' << ptr->second << ',' << (ptr + 1)->second << "\r\n";
            ptr = this->tot_e_extrema();
            out << "system_E," << *ptr->first << ',' << *((ptr + 1)->first) << ',' << *(ptr->first + 1) << ','
                << *((ptr + 1)->first + 1) << ',' << ptr->second << ',' << (ptr + 1)->second << "\r\n";
            out.close();
            if (&path == &def_path)
                def_path.erase_chars(def_path.get_length() - 29);
            // std::cout << "potential energy size: " << pe.size() << std::endl;
            // std::cout << "energy size: " << energy.size() << std::endl;
            // std::cout << "potential energy v size: " << pe[0].size() << std::endl;
            // std::cout << "energy v size: " << energy[0].size() << std::endl;
            // std::cout << "Total potential energy size: " << tot_pe.size() << std::endl;
            // std::cout << "Total kinetic energy size: " << tot_ke.size() << std::endl;
            // std::cout << "Total energy size: " << tot_e.size() << std::endl;
            return true;
        }
        ~system() {
#ifdef GREGSYS_MERGERS
            if constexpr (collisions) {
                if constexpr (memFreq)
                    delete del_bods_m;
                else
                    delete del_bods_n;
            }
#endif
            delete [] G_vals;
        }
        bod_t &operator[](vec_size_t index) {
            if (index >= this->bods.size())
                throw std::out_of_range("The requested body does not exist (index out of range).\n");
            return this->bods[index];
        }
        const bod_t &operator[](vec_size_t index) const {
            if (index >= this->bods.size())
                throw std::out_of_range("The requested body does not exist (index out of range).\n");
            return this->bods[index];
        }
        // implicitly-declared copy assignment operator is deleted, so must define my own:
        sys_t &operator=(const sys_t &other) {
            if (this == &other)
                return *this;
            this->dt = other.dt;
            this->half_dt = other.half_dt;
            this->iterations = other.iterations;
            this->G = other.G;
            this->bods = other.bods;
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq)
                this->clear_evolution();
            return *this;
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t,
                  bool prg, bool mrg, int coll, uint64_t mF, uint64_t fF, bool bF>
        sys_t &operator=(const system<m, r, t, prg, mrg, coll, mF, fF, bF> &other) {
            if constexpr (std::same_as<sys_t, decltype(other)>) {
                if (this == &other)
                    return *this;
            }
            this->dt = other.dt;
            this->half_dt = other.half_dt;
            this->iterations = other.iterations;
            this->G = other.G;
            this->bods.assign(other.bods.begin(), other.bods.end());
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq) {
                this->clear_evolution();
                // if constexpr (mF) {
                //     this->del_bods_m->insert(other.del_bods_m.begin(), other.del_bods_m.end());
                // }
            }
            return *this;
        }
        sys_t &operator=(sys_t &&other) noexcept {
            if (this == &other)
                return *this;
            this->dt = std::move(other.dt);
            this->half_dt = other.half_dt;
            this->iterations = std::move(other.iterations);
            this->G = std::move(other.G);
            this->bods.assign(std::make_move_iterator(other.bods.begin()), std::make_move_iterator(other.bods.end()));
            move(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq)
                this->clear_evolution();
            return *this;
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t,
                  bool prg, bool mrg, int coll, uint64_t mF, uint64_t fF, bool bF>
        friend std::ostream &operator<<(std::ostream&, const system<m, r, t, prg, mrg, coll, mF, fF, bF>&);
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prg1, bool prg2, bool mrg1, bool mrg2,
                  int c1, int c2, uint64_t mF1, uint64_t mF2, uint64_t fF1, uint64_t fF2, bool bF1, bool bF2>
        friend system<m, r, t, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2), bF1 & bF2>
        operator+(const system<m, r, t, prg1, mrg1, c1, mF1, fF1, bF1>&,
                  const system<m, r, t, prg2, mrg2, c2, mF2, fF2, bF2>&);
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, uint64_t, uint64_t, bool>
        friend class system;
        friend class body_tracker<M, R, T, memFreq>;
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                isNumWrapper, bool, uint64_t>
        friend class astro_scene;
    }; // all these template parameters are driving me coocoo
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T,bool prg,bool mrg,int coll,uint64_t mF,uint64_t fF,bool bF>
    std::ostream &operator<<(std::ostream &os, const system<M, R, T, prg, mrg, coll, mF, fF, bF> &sys) {
        typename std::vector<body<M, R, T, mF>>::size_type count = sys.bods.size();
        os << "[gtd::system@" << &sys << ",num_bodies=" << count;
        if (!count)
            return os << ']';
        os << ",bodies:";
        count = 0;
        for (const body<M, R, T, mF> &b : sys)
            os << "\n body_" << count++ << '=' << b;
        return os << ']';
    }
    template<isNumWrapper M, isNumWrapper R, isNumWrapper T, bool prg1, bool prg2, bool mrg1, bool mrg2, int c1, int c2,
             uint64_t mF1, uint64_t mF2, uint64_t fF1, uint64_t fF2, bool bF1, bool bF2>
    system<M, R, T, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2), bF1 & bF2>
    operator+(const system<M, R, T, prg1, mrg1, c1, mF1, fF1, bF1> &sys1,
              const system<M, R, T, prg2, mrg2, c2, mF2, fF2, bF2> &sys2) {
        typename std::vector<body<M, R, T, mF1>>::size_type size1 = sys1.bods.size();
        typename std::vector<body<M, R, T, mF2>>::size_type size2 = sys2.bods.size();
        if ((!size1 && !size2) || sys1.G != sys2.G)
            return {{}};
        if (!size1)
            return {sys2};
        if (!size2)
            return {sys1};
        system<M, R, T, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2), bF1 & bF2>
                ret_sys(sys1.bods, (sys1.dt + sys2.dt)/2.0l, (sys1.iterations + sys2.iterations)/2);
        // for (typename std::vector<body<M, R, T>>::size_type i = 0; i < size1; ++i) {
        //     for (typename std::vector<body<M, R, T>>::size_type j = 0; j < size2; ++j) {
        //         if ((sys1.bods[i].curr_pos-sys2.bods[j].curr_pos).magnitude()<sys1.bods[i].radius+sys2.bods[j].radius) {
        //             if (!system<M, R, T>::merge_if_overlapped)
        //                 throw overlapping_bodies_error();
        //             ret_sys.bods[i];
        //         }
        //     }
        // }
        ret_sys.G = sys1.G;
        copy(ret_sys.G_vals, sys1.G_vals, 6*sizeof(uint64_t));
        ret_sys.add_bodies(sys2.bods);
        return ret_sys;
    }
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t rF>
    class body_tracker {
        /* Rudimentary convenience class - exists solely for the purpose of calculating the centre-of-mass of a subset
         * of bodies within a gtd::system<> object (determined via the unary predicate passed to the ctor below). Note:
         * usage of this class assumes that no mergers will take place within the gtd::system<> object, and that no
         * subsequent bodies will be added to the system (which could result in a reallocation, invalidating all
         * pointers to bodies within this class). */
        using bod_t = body<M, R, T, rF>;
        template <bool prog, bool merge, int coll, uint64_t fF, bool bin>
        using sys_t = system<M, R, T, prog, merge, coll, rF, fF, bin>;
        using vec = vector3D<decltype(std::declval<M>()*std::declval<T>()/std::declval<M>())>;
        using vec_size_t = typename std::vector<const bod_t*>::size_type;
        std::vector<const bod_t*> _bods; // pointers to all the bodies
        // vec _com; // centre of mass of all bodies
        // M _tmass; // total mass of all bodies
        // void calc_mass() {
        //
        // }
        // void calc_com() {
        //     this->calc_mass();
        // }
    public:
        template <bool prog, bool merge, int coll, uint64_t fF, bool bin, unaryPredicate<bod_t> func>
        body_tracker(const sys_t<prog, merge, coll, fF, bin> &sys, const func &_f, vec_size_t reserve_bods = 0) {
            if (sys.bods.empty())
                return;
            if (reserve_bods)
                _bods.reserve(reserve_bods);
            for (const bod_t &_b : sys)
                if (_f(_b))
                    _bods.push_back(&_b);
        }
        uint64_t num_bods() const noexcept {
            return this->_bods.size();
        }
        M mass() const {
            M _tmass{}; // total mass of bodies
            for (const bod_t* const &_ptr : _bods)
                _tmass += _ptr->mass_;
            return _tmass;
        }
        vec com() const {
            if (_bods.empty())
                return {};
            M _tmass{}; // total mass of bodies
            vec _mr_sum;
            for (const bod_t* const &_ptr : _bods)
                _tmass += _ptr->mass_;
            for (const bod_t* const &_ptr : _bods)
                _mr_sum += _ptr->mass_*_ptr->curr_pos;
            return _tmass == M{} ? vec{} : _mr_sum / _tmass; // guaranteed copy elision
        }
        vec com_vel() const {
            if (_bods.empty())
                return {};
            M _tmass{};
            vec _mv_sum;
            for (const bod_t* const &_ptr : _bods)
                _tmass += _ptr->mass_;
            for (const bod_t* const &_ptr : _bods)
                _mv_sum += _ptr->mass_*_ptr->curr_vel;
            return _tmass == M{} ? vec{} : _mv_sum / _tmass;
        }
        T mean_dist_to(const vec &point) const {
            if (_bods.empty())
                return {};
            T sum{};
            for (const bod_t *const &_b : _bods)
                sum += (_b->curr_pos - point).magnitude();
            return sum/_bods.size();
        }
        T mean_sep() const {
            if (_bods.empty())
                throw nbody_error{"Error: no bodies being tracked.\n"};
            typename std::vector<bod_t*>::size_type _size = _bods.size();
            if (_size == 1)
                return {};
            T _tot_dist{};
            const bod_t *const *outer = _bods.data();
            const bod_t *const *inner{};
            const bod_t *const *const _end = outer + _size;
            while (outer < _end) {
                inner = outer + 1;
                while (inner < _end) {
                    _tot_dist += vec_ops::distance((*inner)->curr_pos, (*outer)->curr_pos);
                    ++inner;
                }
                ++outer;
            }
            return _tot_dist/((_size*(_size - 1))/2);
        }
        const bod_t *const *front() const noexcept {
            if (_bods.empty())
                return nullptr;
            return _bods.data();
        }
        const bod_t *const *back() const noexcept {
            if (_bods.empty())
                return nullptr;
            return _bods.data() + _bods.size() - 1;
        }
        auto begin() {
            return this->_bods.begin();
        }
        auto end() {
            return this->_bods.end();
        }
        auto begin() const {
            return this->_bods.cbegin();
        }
        auto end() const {
            return this->_bods.cend();
        }
        bod_t *&operator[](typename std::vector<const bod_t*>::size_type index) { // no bounds-checking performed
            return this->_bods[index];
        }
        const bod_t *const &operator[](typename std::vector<const bod_t*>::size_type index) const { // here neither
            return this->_bods[index];
        }
    };
    // std::set<unsigned long long> body_counter::ids;
    // unsigned long long body_counter::count = 0;
    typedef system<long double, long double, long double> sys;
    typedef system<long double, long double, long double, true> sys_p;
    typedef system<long double, long double, long double, false, true> sys_m;
    typedef system<long double, long double, long double, false, false, 3> sys_c;
    typedef system<long double, long double, long double, false, false, 7> sys_C;
    typedef system<long double, long double, long double, true, true> sys_pm;
    typedef system<long double, long double, long double, true, false, 3> sys_pc;
    typedef system<long double, long double, long double, true, false, 7> sys_pC;
    typedef system<long double, long double, long double, false, true, 3> sys_mc;
    typedef system<long double, long double, long double, false, true, 7> sys_mC;
    typedef system<long double, long double, long double, true, true, 3> sys_pmc;
    typedef system<long double, long double, long double, true, true, 7> sys_pmC;
    typedef system<long double, long double, long double, false, false, 0, 0, 0, false> sys_n;
    typedef system<long double, long double, long double, true, false, 0, 0, 0, false> sys_p_n;
    typedef system<long double, long double, long double, false, true, 0, 0, 0, false> sys_m_n;
    typedef system<long double, long double, long double, false, false, 3, 0, 0, false> sys_c_n;
    typedef system<long double, long double, long double, false, false, 7, 0, 0, false> sys_C_n;
    typedef system<long double, long double, long double, true, true, 0, 0, 0, false> sys_pm_n;
    typedef system<long double, long double, long double, true, false, 3, 0, 0, false> sys_pc_n;
    typedef system<long double, long double, long double, true, false, 7, 0, 0, false> sys_pC_n;
    typedef system<long double, long double, long double, false, true, 3, 0, 0, false> sys_mc_n;
    typedef system<long double, long double, long double, false, true, 7, 0, 0, false> sys_mC_n;
    typedef system<long double, long double, long double, true, true, 3, 0, 0, false> sys_pmc_n;
    typedef system<long double, long double, long double, true, true, 7, 0, 0, false> sys_pmC_n;
    typedef system<long double, long double, long double, false, false, 0, 0, 1, false> sys_txt;
    typedef system<long double, long double, long double, true, false, 0, 0, 1, false> sys_p_txt;
    typedef system<long double, long double, long double, false, true, 0, 0, 1, false> sys_m_txt;
    typedef system<long double, long double, long double, false, false, 3, 0, 1, false> sys_c_txt;
    typedef system<long double, long double, long double, false, false, 7, 0, 1, false> sys_C_txt;
    typedef system<long double, long double, long double, true, true, 0, 0, 1, false> sys_pm_txt;
    typedef system<long double, long double, long double, true, false, 3, 0, 1, false> sys_pc_txt;
    typedef system<long double, long double, long double, true, false, 7, 0, 1, false> sys_pC_txt;
    typedef system<long double, long double, long double, false, true, 3, 0, 1, false> sys_mc_txt;
    typedef system<long double, long double, long double, false, true, 7, 0, 1, false> sys_mC_txt;
    typedef system<long double, long double, long double, true, true, 3, 0, 1, false> sys_pmc_txt;
    typedef system<long double, long double, long double, true, true, 7, 0, 1, false> sys_pmC_txt;
    typedef system<long double, long double, long double, false, false, 0, 1, 0> sys_M;
    typedef system<long double, long double, long double, true, false, 0, 1, 0> sys_p_M;
    typedef system<long double, long double, long double, false, true, 0, 1, 0> sys_m_M;
    typedef system<long double, long double, long double, false, false, 3, 1, 0> sys_c_M;
    typedef system<long double, long double, long double, false, false, 7, 1, 0> sys_C_M;
    typedef system<long double, long double, long double, true, true, 0, 1, 0> sys_pm_M;
    typedef system<long double, long double, long double, true, false, 3, 1, 0> sys_pc_M;
    typedef system<long double, long double, long double, true, false, 7, 1, 0> sys_pC_M;
    typedef system<long double, long double, long double, false, true, 3, 1, 0> sys_mc_M;
    typedef system<long double, long double, long double, false, true, 7, 1, 0> sys_mC_M;
    typedef system<long double, long double, long double, true, true, 3, 1, 0> sys_pmc_M;
    typedef system<long double, long double, long double, true, true, 7, 1, 0> sys_pmC_M;
    typedef system<long double, long double, long double, false, false, 0, 1, 1> sys_MF;
    typedef system<long double, long double, long double, true, false, 0, 1, 1> sys_p_MF;
    typedef system<long double, long double, long double, false, true, 0, 1, 1> sys_m_MF;
    typedef system<long double, long double, long double, false, false, 3, 1, 1> sys_c_MF;
    typedef system<long double, long double, long double, false, false, 7, 1, 1> sys_C_MF;
    typedef system<long double, long double, long double, true, true, 0, 1, 1> sys_pm_MF;
    typedef system<long double, long double, long double, true, false, 3, 1, 1> sys_pc_MF;
    typedef system<long double, long double, long double, true, false, 7, 1, 1> sys_pC_MF;
    typedef system<long double, long double, long double, false, true, 3, 1, 1> sys_mc_MF;
    typedef system<long double, long double, long double, false, true, 7, 1, 1> sys_mC_MF;
    typedef system<long double, long double, long double, true, true, 3, 1, 1> sys_pmc_MF;
    typedef system<long double, long double, long double, true, true, 7, 1, 1> sys_pmC_MF;
    typedef system<long double, long double, long double, false, false, 0, 1, 1, false> sys_MF_txt;
    typedef system<long double, long double, long double, true, false, 0, 1, 1, false> sys_p_MF_txt;
    typedef system<long double, long double, long double, false, true, 0, 1, 1, false> sys_m_MF_txt;
    typedef system<long double, long double, long double, false, false, 3, 1, 1, false> sys_c_MF_txt;
    typedef system<long double, long double, long double, false, false, 7, 1, 1, false> sys_C_MF_txt;
    typedef system<long double, long double, long double, true, true, 0, 1, 1, false> sys_pm_MF_txt;
    typedef system<long double, long double, long double, true, false, 3, 1, 1, false> sys_pc_MF_txt;
    typedef system<long double, long double, long double, true, false, 7, 1, 1, false> sys_pC_MF_txt;
    typedef system<long double, long double, long double, false, true, 3, 1, 1, false> sys_mc_MF_txt;
    typedef system<long double, long double, long double, false, true, 7, 1, 1, false> sys_mC_MF_txt;
    typedef system<long double, long double, long double, true, true, 3, 1, 1, false> sys_pmc_MF_txt;
    typedef system<long double, long double, long double, true, true, 7, 1, 1, false> sys_pmC_MF_txt;
    typedef body_tracker<long double, long double, long double, 0> btrk_0f;
    typedef body_tracker<long double, long double, long double, 10> btrk_10f;
    typedef body_tracker<long double, long double, long double, 10> btrk_10f;
    typedef body_tracker<long double, long double, long double, 100> btrk_100f;
    typedef body_tracker<long double, long double, long double, 1000> btrk_1000f;
    template <uint64_t rF>
    using btrk = body_tracker<long double, long double, long double, rF>;
}
#undef MEM_LOOP
#undef FUNC_TEMPL_SELECT
#undef EMPTY
#endif
