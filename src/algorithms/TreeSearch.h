/* Copyright 2024 Henning Woydt. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/

#ifndef SUBSETOPTIMIZATION_TREESEARCH_H
#define SUBSETOPTIMIZATION_TREESEARCH_H

#define DEPTH_STAY 0
#define DEPTH_UP 1
#define DEPTH_DOWN 2

#include <chrono>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <cstring>
#include <numeric>
#include <utility>
#include <vector>
#include <cmath>
#include <typeinfo>

#include "../utility/Util.h"
#include "../utility/AlgorithmConfiguration.h"
#include "../utility/VectorOfVectors.h"
#include "UpperBoundManager.h"
#include "CandidateManager.h"
#include "SICache.h"
#include "PBFAlgorithm.h"
#include "UB2DAlgorithm.h"

namespace SMSM {

/**
 * Algorithm to find the set S with k elements that maximizes a score function.
 *
 * @tparam T The structure holding the n elements.
 */
    template<class T, typename TypeSF>
    class TreeSearch {
    public:
        T &t;
        size_t n;
        size_t k;
        AlgorithmConfiguration ac;

        // vectors to store s
        std::vector<uint32_t> s;
        std::vector<uint32_t> initial_s;
        std::vector<uint32_t> best_s;
        std::vector<uint32_t> best_greedy_s;
        size_t s_size;
        size_t tree_max_depth;

        // scores
        TypeSF initial_score = -std::numeric_limits<TypeSF>::max();
        TypeSF best_score = -std::numeric_limits<TypeSF>::max();
        TypeSF best_greedy_score = -std::numeric_limits<TypeSF>::max();

        std::vector<TypeSF> scores;

        // vectors to store 1d managers
        std::vector<CandidateManager<TypeSF>> c_managers;

        // vectors for the UB2D managers
        std::vector<UpperBoundManager<TypeSF>> ub2d_managers;

        // vectors to store the UB2D Interfaces
        std::vector<UB2DAlgorithm<TypeSF>> ub2d_algorithms;

        // vectors for the PBF managers
        std::vector<UpperBoundManager<TypeSF>> pbf_managers;

        // vectors for the PBF algorithms
        std::vector<PBFAlgorithm<TypeSF>> pbf_algorithms;

        // holds the object for recursive call
        AlgorithmConfiguration rec_ac;
        std::vector<UpperBoundManager<TypeSF>> rec_managers;

        // vectors to hold SICache
        std::vector<SICache<TypeSF>> si_caches;

        // Temporary vectors
        VectorOfVectors<size_t> p1_sets_even;
        std::vector<TypeSF> p1_si_even;
        VectorOfVectors<size_t> p1_sets_a;
        std::vector<TypeSF> p1_si_a;
        VectorOfVectors<size_t> p1_sets_b;
        std::vector<TypeSF> p1_si_b;
        std::vector<size_t> set;
        std::vector<size_t> comb;

        // Variables to measure oracles
        // UB2D
        std::chrono::time_point<std::chrono::steady_clock> sp_ub2d;
        double UB2D_success_time = 0.0;
        double UB2D_failure_time = 0.0;
        size_t UB2D_n_success = 0;
        size_t UB2D_n_failure = 0;

        // Partial Brute Force
        std::chrono::time_point<std::chrono::steady_clock> sp_pbf;
        double PBF_success_time = 0.0;
        double PBF_failure_time = 0.0;
        size_t PBF_n_success = 0;
        size_t PBF_n_failure = 0;

        // Recursive
        std::chrono::time_point<std::chrono::steady_clock> sp_rec;
        double REC_success_time = 0.0;
        double REC_failure_time = 0.0;
        size_t REC_n_success = 0;
        size_t REC_n_failure = 0;

        // variables needed for timing
        double time_needed = 0.0;
        std::chrono::steady_clock::time_point time_sp;
        size_t n_tried_calls = 0;
        bool time_exceeded = false;

        /**
         * Constructor.
         *
         * @param t The structure (it will be used with the score function).
         * @param ac Algorithm Configuration.
         */
        TreeSearch(T &t,
                   AlgorithmConfiguration ac_temp) : t(t), ac(std::move(ac_temp)) {
            n = t.get_n();
            k = ac.k;

            // allocate space for s
            s.resize(k);
            initial_s.resize(k);
            best_s.resize(k);
            best_greedy_s.resize(k);
            size_t t_original_depth = t.depth;
            for (size_t i = 0; i < ac.partial_s.size(); ++i) {
                s[i] = ac.partial_s[i];
                initial_s[i] = ac.partial_s[i];
                best_s[i] = ac.partial_s[i];
                best_greedy_s[i] = ac.partial_s[i];
                if (t_original_depth == 0) {
                    t.visit_new_depth(s, i + 1);
                }
            }
            s_size = ac.partial_s.size();
            tree_max_depth = k - s_size;

            scores.resize(tree_max_depth, -std::numeric_limits<TypeSF>::max());

            // allocate space for Candidate Managers
            for (size_t i = 0; i < tree_max_depth; ++i) { c_managers.emplace_back(n); }
            c_managers[0].fill_candidates(n, ac.partial_s, ac.init_candidates, ac.init_si, ac.init_acc);

            if (ac.UB2D_enabled) {
                // allocate space for the UB2D Manager
                for (size_t i = 0; i < tree_max_depth; ++i) { ub2d_managers.emplace_back(ac.UB2D_lazy_skip_start_value, ac.UB2D_lazy_skip_add_value, ac.UB2D_safe_skip_enabled); }

                // allocate space for UB2D Interfaces
                size_t max_l = ac.determine_UB2D_max_l(n, k);
                for (size_t i = 0; i < tree_max_depth; ++i) { ub2d_algorithms.emplace_back(max_l, ac.UB2D_alg_type); }
            }

            if (ac.PBF_enabled) {
                // allocate space for pbf managers
                for (size_t i = 0; i < tree_max_depth; ++i) { pbf_managers.emplace_back(ac.PBF_lazy_skip_start_value, ac.PBF_lazy_skip_add_value, ac.PBF_safe_skip_enabled); }

                // allocate space for pbf algorithms
                size_t max_n = ac.determine_PBF_max_n(n, k);
                size_t max_l = ac.determine_PBF_max_l(n, k);
                for (size_t i = 0; i < tree_max_depth; ++i) { pbf_algorithms.emplace_back(max_n, max_l, k, ac.PBF_alg_type, ac.PBF_safe_skip_enabled); }
            }

            if (ac.REC_enabled) {
                // allocate space for pbf managers
                for (size_t i = 0; i < tree_max_depth; ++i) { rec_managers.emplace_back(ac.REC_lazy_skip_start_value, ac.REC_lazy_skip_add_value, ac.REC_safe_skip_enabled); }

            }

            if (ac.UB2D_enabled || ac.PBF_enabled) {
                // allocate space for caches
                for (size_t i = 0; i < k; ++i) { si_caches.emplace_back(); }
            }

            if (ac.REC_enabled) {
                rec_ac = ac;
                rec_ac.REC_manual_depth_enabled.clear();
                rec_ac.write_output = false;
                rec_ac.REC_enabled = false;
                rec_ac.measure_oracle_time = false;
            }

            for (size_t i = s_size; i < k; ++i) {
                initial_s[i] = c_managers[0].get_c(i - s_size);
            }
        };

        /**
         * Will start the search algorithm and return the best found set. The set
         * will be sorted ascending.
         *
         * @return The best found set.
         */
        std::vector<uint32_t> search() {
            time_sp = get_time_point();

            // initial solution
            initial_score = t.evaluate_general(initial_s, k);
            std::copy(initial_s.begin(), initial_s.end(), best_s.begin());
            best_score = initial_score;

            // initialize the score of the empty set
            scores[0] = t.evaluate_general(s, s_size);

            // special vars
            size_t r = k - s_size;
            size_t n_remaining = c_managers[0].size - c_managers[0].offset;

            // Dynamic Candidate Ordering on Depth 0
            c_managers[0].size = DCO_depth_0();

            // special cases
            if (r == 1) {
                best_s[s_size] = c_managers[0].get_c(0);
                best_score = t.evaluate_general(best_s, k);

            } else if (n_remaining == r) {
                r_remaining_candidates(0);
            } else if (n_remaining == r + 1) {
                one_candidate_too_much(1);
            } else if (n_remaining <= ac.bf_threshold_n && r <= ac.bf_threshold_r) {
                brute_force(1);
            } else {
                // initial greedy solution
                if (ac.greedy_local_search) {
                    greedy_local_search();
                } else {
                    simple_greedy();
                }

                // initialize structures
                if (ac.UB2D_enabled) { ub2d_managers[0].visit_new_depth(); }
                if (ac.PBF_enabled) { pbf_managers[0].visit_new_depth(); }
                if (ac.REC_enabled) { rec_managers[0].visit_new_depth(); }
                if (ac.UB2D_enabled || ac.PBF_enabled) { si_caches[0].visit_new_depth(); }

                // search through the candidates
                iterative_search();
            }

            // stop time
            time_needed = get_elapsed_seconds(time_sp, get_time_point());

            // write if specified
            if (ac.write_output) {
                write_to_JSON();
            }

            return best_s;
        };

        /**
         * Iteratively searches for a solution.
         */
        void iterative_search() {
            // init vars
            size_t depth = 0;
            size_t r = k - s_size;
            TypeSF score = scores[depth];
            TypeSF r_score = best_score - score;
            size_t depth_action = DEPTH_STAY;

            while (depth != 0 || depth_action != DEPTH_UP) {

                if (depth_action == DEPTH_STAY) {
                    bool close = false;
                    if constexpr (std::is_same<TypeSF, double>::value) {
                        close = double_eq(best_score, ac.score_threshold, 0.0001);
                    }
                    if (best_score > ac.score_threshold && !close) {
                        // we have found a set with more than the desired score
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    if (best_score == t.max_reachable_score) {
                        // we have found a set with the best score
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    // depth, s_size, and r assumed to have the correct value
                    if (c_managers[depth].offset + r > c_managers[depth].size) {
                        // if not enough elements go up
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    TypeSF sub_bound = std::numeric_limits<TypeSF>::max();
                    TypeSF ub2d_bound = std::numeric_limits<TypeSF>::max();
                    TypeSF pbf_bound = std::numeric_limits<TypeSF>::max();
                    TypeSF rec_bound = std::numeric_limits<TypeSF>::max();

                    // Simple Upper Bound
                    if (ac.SUB_enabled) {
                        sub_bound = SUB(depth);

                        if (sub_bound <= best_score) {
                            depth_action = DEPTH_UP;
                            continue;
                        }
                    }

                    // Reduction of Possible Candidates
                    if (ac.RPC_enabled) {
                        c_managers[depth].size = RPC(depth);
                        if (c_managers[depth].size - c_managers[depth].offset == r) {
                            r_remaining_candidates(depth);
                            depth_action = DEPTH_UP;
                            continue;
                        }
                    }

                    // Upper Bound 2D
                    if (ac.UB2D_enabled && ac.UB2D_depth_enabled[depth] && r_score >= (sub_bound - scores[depth]) * ac.UB2D_sub_bound_percentage_vec[depth]) {

                        if (ac.measure_oracle_time) {
                            sp_ub2d = get_time_point();
                        }

                        ub2d_bound = UB2D(depth);
                        if (ub2d_bound <= best_score) {

                            if (ac.measure_oracle_time) {
                                UB2D_success_time += get_elapsed_seconds(sp_ub2d, get_time_point());
                                UB2D_n_success += 1;
                            }

                            depth_action = DEPTH_UP;
                            continue;
                        }

                        UpperBoundManager<TypeSF> &ubm_ub2d = ub2d_managers[depth];
                        if (ac.UB2D_RPC_enabled && ubm_ub2d.has_r1_solution()) {
                            size_t largest_idx = ubm_ub2d.get_largest_idx_of_r1_solution();
                            TypeSF r1_si_bound = ubm_ub2d.get_r1_si_bound();

                            c_managers[depth].size = general_RPC(r1_si_bound, largest_idx, depth);
                            if (c_managers[depth].size - c_managers[depth].offset == r) {
                                r_remaining_candidates(depth);

                                if (ac.measure_oracle_time) {
                                    UB2D_success_time += get_elapsed_seconds(sp_ub2d, get_time_point());
                                    UB2D_n_success += 1;
                                }

                                depth_action = DEPTH_UP;
                                continue;
                            }
                        }

                        if (ac.measure_oracle_time) {
                            UB2D_failure_time += get_elapsed_seconds(sp_ub2d, get_time_point());
                            UB2D_n_failure += 1;
                        }
                    }

                    // Partial Brute Force
                    if (ac.PBF_enabled && ac.PBF_depth_enabled[depth] && r_score >= (sub_bound - scores[depth]) * ac.PBF_sub_bound_percentage_vec[depth]) {
                        if (ac.measure_oracle_time) {
                            sp_pbf = get_time_point();
                        }

                        pbf_bound = PBF(depth);
                        if (pbf_bound <= best_score) {

                            if (ac.measure_oracle_time) {
                                PBF_success_time += get_elapsed_seconds(sp_pbf, get_time_point());
                                PBF_n_success += 1;
                            }

                            depth_action = DEPTH_UP;
                            continue;
                        }

                        UpperBoundManager<TypeSF> &ubm_pbf = pbf_managers[depth];
                        if (ac.PBF_RPC_enabled && ubm_pbf.has_r1_solution()) {
                            size_t largest_idx = ubm_pbf.get_largest_idx_of_r1_solution();
                            TypeSF r1_si_bound = ubm_pbf.get_r1_si_bound();

                            c_managers[depth].size = general_RPC(r1_si_bound, largest_idx, depth);
                            if (c_managers[depth].size - c_managers[depth].offset == r) {
                                r_remaining_candidates(depth);

                                if (ac.measure_oracle_time) {
                                    PBF_success_time += get_elapsed_seconds(sp_pbf, get_time_point());
                                    PBF_n_success += 1;
                                }

                                depth_action = DEPTH_UP;
                                continue;
                            }
                        }

                        if (ac.measure_oracle_time) {
                            PBF_failure_time += get_elapsed_seconds(sp_pbf, get_time_point());
                            PBF_n_failure += 1;
                        }
                    }

                    // Recursive heuristic
                    if (ac.REC_enabled && ac.REC_depth_enabled[depth] && r_score >= (sub_bound - scores[depth]) * ac.REC_sub_bound_percentage_vec[depth]) {
                        if (ac.measure_oracle_time) {
                            sp_rec = get_time_point();
                        }

                        rec_bound = REC(depth);
                        if (rec_bound <= best_score) {

                            if (ac.measure_oracle_time) {
                                REC_success_time += get_elapsed_seconds(sp_rec, get_time_point());
                                REC_n_success += 1;
                            }

                            depth_action = DEPTH_UP;
                            continue;
                        }

                        UpperBoundManager<TypeSF> &ubm_rec = rec_managers[depth];
                        if (ac.REC_RPC_enabled && ubm_rec.has_r1_solution()) {
                            size_t largest_idx = ubm_rec.get_largest_idx_of_r1_solution();
                            TypeSF r1_si_bound = ubm_rec.get_r1_si_bound();

                            c_managers[depth].size = general_RPC(r1_si_bound, largest_idx, depth);
                            if (c_managers[depth].size - c_managers[depth].offset == r) {
                                r_remaining_candidates(depth);

                                if (ac.measure_oracle_time) {
                                    REC_success_time += get_elapsed_seconds(sp_rec, get_time_point());
                                    REC_n_success += 1;
                                }

                                depth_action = DEPTH_UP;
                                continue;
                            }
                        }

                        if (ac.measure_oracle_time) {
                            REC_failure_time += get_elapsed_seconds(sp_rec, get_time_point());
                            REC_n_failure += 1;
                        }
                    }

                    // no heuristic pruned, so go down the tree
                    s[s_size] = c_managers[depth].get_c(c_managers[depth].offset);
                    scores[depth + 1] = c_managers[depth].get_acc(c_managers[depth].offset) ? scores[depth] + c_managers[depth].get_si(c_managers[depth].offset) : t.evaluate_1D(s, s_size + 1);
                    c_managers[depth].offset += 1;

                    t.visit_new_depth(s, s_size + 1);
                    depth_action = DEPTH_DOWN;
                    continue;
                } else if (depth_action == DEPTH_UP) {
                    // we will go one depth up
                    // depth, s_size, and r are in the wrong state

                    depth -= 1;
                    s_size -= 1;
                    r = k - s_size;
                    score = scores[depth];
                    r_score = best_score - score;

                    t.return_from_last_depth();

                    if (has_time_exceeded()) {
                        return;
                    }

                    depth_action = DEPTH_STAY;
                    continue;
                } else {
                    // we will go one depth down
                    // depth, s_size, and r are in the wrong state

                    depth += 1;
                    s_size += 1;
                    r = k - s_size;
                    score = scores[depth];
                    r_score = best_score - score;

                    size_t p_c_size = c_managers[depth - 1].size;
                    size_t p_offset = c_managers[depth - 1].offset;

                    // special case
                    if (r == 1) {
                        one_element_remaining(depth);
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    // special case
                    if (p_c_size - p_offset == r) {
                        r_remaining_candidates(depth - 1);
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    // special case
                    if (p_c_size - p_offset == r + 1) {
                        one_candidate_too_much(depth);
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    // special case
                    if (p_c_size - p_offset <= ac.bf_threshold_n && r <= ac.bf_threshold_r) {
                        brute_force(depth);
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    // signal to structures that a new depth was reached
                    if (ac.UB2D_enabled) { ub2d_managers[depth].visit_new_depth(); }
                    if (ac.PBF_enabled) { pbf_managers[depth].visit_new_depth(); }
                    if (ac.REC_enabled) { rec_managers[depth].visit_new_depth(); }
                    if (ac.UB2D_enabled || ac.PBF_enabled) { si_caches[depth].visit_new_depth(); }

                    // dynamic candidate ordering
                    c_managers[depth].size = DCO(depth);

                    if (c_managers[depth].size < r) {
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    c_managers[depth].offset = 0;

                    depth_action = DEPTH_STAY;
                    continue;
                }
            }
        }

        /**
        * Special function to handle the case when only one more element has to be
        * chosen. It will look at the candidates from the parent and determine the
        * best one.
        *
        * @param p_c_size Size of the parent candidates.
        * @param p_offset Offset into the parent candidates.
        * @param s_size Current size of set S.
        * @param score Current score of set S.
        * @param depth The current depth of the search tree.
        */
        void one_element_remaining(const size_t depth) {
            const CandidateManager<TypeSF> &pc_manager = c_managers[depth - 1];
            const size_t p_size = pc_manager.size;
            const size_t p_offset = pc_manager.offset;

            const TypeSF score = scores[depth];
            TypeSF remaining_score = best_score - score;

            for (size_t i = p_offset; i < p_size; ++i) {
                if (pc_manager.get_si(i) >= remaining_score) {
                    s[s_size] = pc_manager.get_c(i);
                    TypeSF new_score = t.evaluate_1D(s, s_size + 1);

                    if (new_score > best_score) {
                        update_best(s, new_score);
                        remaining_score = best_score - score;
                    }
                } else {
                    break;
                }
            }
        };

        /**
         * Special function, that handles the case if the number of remaining
         * candidates is the same as the number of remaining elements in set S.
         *
         * @param offset Offset into the candidates.
         * @param s_size The current size of set S.
         * @param depth The current depth of the search tree.
         */
        void r_remaining_candidates(const size_t depth) {
            const CandidateManager<TypeSF> &c1d_manager = c_managers[depth];
            const size_t offset = c1d_manager.offset;
            const size_t r = k - s_size;

            for (size_t i = 0; i < r; ++i) {
                s[s_size + i] = c1d_manager.get_c(offset + i);
            }

            TypeSF score = t.evaluate_general(s, k);

            if (score > best_score) {
                update_best(s, score);
            }
        }

        /**
         * Special function, that handles the case if the number of remaining
         * candidates - 1 is the same as the number of remaining elements in set S.
         *
         * @param offset Offset into the candidates.
         * @param s_size The current size of set S.
         * @param depth The current depth of the search tree.
         */
        void one_candidate_too_much(const size_t depth) {
            const CandidateManager<TypeSF> &pc_manager = c_managers[depth - 1];
            const size_t p_offset = pc_manager.offset;
            const size_t r = k - s_size;

            for (size_t i = 0; i < r; ++i) {
                s[s_size + i] = pc_manager.get_c(p_offset + i);
            }
            TypeSF new_score = t.evaluate_general(s, s_size + r);
            if (new_score > best_score) {
                update_best(s, new_score);
            }

            for (size_t i = 0; i < r; ++i) {
                s[s_size + r - 1 - i] = pc_manager.get_c(p_offset + r - i);
                new_score = t.evaluate_general(s, s_size + r);

                if (new_score > best_score) {
                    update_best(s, new_score);
                }
            }
        }

        /**
         * Special function, that brute forces all possible combinations.
         *
         * @param p_c_size Size of the parent candidates.
         * @param p_offset Offset into the parent candidates.
         * @param s_size The current size of set S.
         * @param depth The current depth of the search tree.
         */
        void brute_force(const size_t depth) {
            const CandidateManager<TypeSF> &pc_manager = c_managers[depth - 1];
            const size_t p_size = pc_manager.size;
            const size_t p_offset = pc_manager.offset;

            const size_t n_remaining = p_size - p_offset;
            const size_t r = k - s_size;

            set.resize(r);
            std::iota(set.begin(), set.end(), 0);
            set[r - 1] -= 1;

            while (next_subset(set, r, n_remaining)) {
                for (size_t i = 0; i < r; ++i) {
                    s[s_size + i] = pc_manager.get_c(p_offset + set[i]);
                }
                TypeSF score = t.evaluate_general(s, k);

                if (score > best_score) {
                    update_best(s, score);
                }
            }
        }

        /**
         * Updates the best found element.
         *
         * @param temp New best set.
         * @param new_score New best score.
         * @param depth The depth at which the set was found.
         */
        void update_best(const std::vector<uint32_t> &temp, const TypeSF new_score) {
            std::copy(temp.begin(), temp.end(), best_s.begin());
            best_score = new_score;
        }

        /**
         * A simple greedy algorithm that tries to find a good starting solution. It
         * iteratively chooses the best element k times.
         */
        void simple_greedy() {
            CandidateManager<TypeSF> &c_manager = c_managers[0];

            for (size_t i = s_size; i < k; ++i) {
                uint32_t best_c = 0;

                for (size_t j = c_manager.offset; j < c_manager.size; ++j) {
                    uint32_t c = c_manager.get_c(j);
                    if (!contains(s, i, c)) {
                        s[i] = c;
                        TypeSF new_score = t.evaluate_general(s, i + 1);

                        if (new_score >= best_greedy_score) {
                            best_greedy_score = new_score;
                            best_c = c;
                        }

                        if (has_time_exceeded()) {
                            return;
                        }
                    }
                }
                s[i] = best_c;
            }

            std::copy(s.begin(), s.end(), best_greedy_s.begin());

            if (best_greedy_score > best_score) {
                best_score = best_greedy_score;
                std::copy(s.begin(), s.end(), best_s.begin());
            }
        }

        /**
         * A simple greedy algorithm that tries to find a good starting solution. It
         * iteratively chooses the best element k times.
         */
        void greedy_local_search() {
            CandidateManager<TypeSF> &c_manager = c_managers[0];

            for (size_t i = s_size; i < k; ++i) {
                uint32_t best_c = 0;

                for (size_t j = c_manager.offset; j < c_manager.size; ++j) {
                    uint32_t c = c_manager.get_c(j);
                    if (!contains(s, i, c)) {
                        s[i] = c;
                        TypeSF new_score = t.evaluate_general(s, i + 1);

                        if (new_score >= best_greedy_score) {
                            best_greedy_score = new_score;
                            best_c = c;
                        }

                        if (has_time_exceeded()) {
                            return;
                        }
                    }
                }
                s[i] = best_c;
            }

            std::copy(s.begin(), s.end(), best_greedy_s.begin());

            if (best_greedy_score > best_score) {
                best_score = best_greedy_score;
                std::copy(s.begin(), s.end(), best_s.begin());
            }

            // now perform a local search, with 5% of the time limit
            double time_limit = std::min(90.0, ac.time_limit * 0.1);
            auto sp = get_time_point();
            bool better_set_found = true;

            while (better_set_found && get_elapsed_seconds(sp, get_time_point()) < time_limit) {
                better_set_found = false;
                for (size_t i = s_size; i < k; ++i) {
                    uint32_t best_c = s[i];

                    for (size_t j = c_manager.offset; j < c_manager.size; ++j) {
                        uint32_t c = c_manager.get_c(j);
                        if (!contains(s, k, c)) {
                            s[i] = c;
                            TypeSF new_score = t.evaluate_general(s, k);

                            if (new_score > best_greedy_score) {
                                best_greedy_score = new_score;
                                best_c = c;
                                better_set_found = true;

                                std::copy(s.begin(), s.end(), best_greedy_s.begin());
                                if (best_greedy_score > best_score) {
                                    best_score = best_greedy_score;
                                    std::copy(s.begin(), s.end(), best_s.begin());
                                }

                            }

                            if (has_time_exceeded() || get_elapsed_seconds(sp, get_time_point()) > time_limit) {
                                return;
                            }
                        }
                    }
                    s[i] = best_c;
                }
            }
        }

        /**
         * Dynamic Candidate Ordering at depth 0. It will calculate the score
         * improvement for all candidates and sort them descending.
         *
         * @param score The current score.
         * @return The size of the candidate set.
         */
        size_t DCO_depth_0() {

            CandidateManager<TypeSF> &c_manger = c_managers[0];
            const TypeSF score = scores[0];

            for (size_t i = c_manger.offset; i < c_manger.size; ++i) {
                if (c_manger.get_acc(i) == 0) {
                    uint32_t c = c_manger.get_c(i);
                    s[s_size] = c;
                    TypeSF si = t.evaluate_1D(s, s_size + 1) - score;
                    c_manger.set_entry(i, c, si, 1);
                }
            }

            c_manger.sort();
            c_manger.calc_csum();

            return c_manger.size;
        };

        /**
         * Calculates the 1D score improvement for each candidate and then reorders
         * them so the improvements are ordered descending. It will iterate through
         * the parents candidates and decide for each candidate if the score
         * improvement should be recomputed. This can be controlled via the 'SLU'
         * options.
         *
         * @param p_c_size Size of the parent candidates.
         * @param p_offset Offset into the parent candidates.
         * @param s_size Current size of set S.
         * @param score Current score of set S.
         * @param depth The current depth of the search tree.
         */
        size_t DCO(const size_t depth) {
            const CandidateManager<TypeSF> &pc_manager = c_managers[depth - 1];
            const size_t p_c_size = pc_manager.size;
            const size_t p_offset = pc_manager.offset;

            const size_t r = k - s_size;
            const TypeSF score = scores[depth];
            const TypeSF r_score = best_score - score;
            const double r_score_avg = (double) r_score / (double) r;

            const double score_threshold = ac.determine_LE_score_threshold(r_score_avg);
            const size_t rank_threshold = ac.determine_LE_rank_threshold(n, k, p_c_size - p_offset, r);

            CandidateManager<TypeSF> &c_manager = c_managers[depth];
            c_manager.clear(r);

            // default: update all
            auto update_scheme = [&](size_t curr_rank, double curr_score) {
                if (ac.LE_mode == 0) {
                    return true;
                } else if (ac.LE_mode == 1) {
                    return curr_score >= score_threshold || curr_rank <= rank_threshold;
                } else if (ac.LE_mode == 2) {
                    return curr_score >= score_threshold && curr_rank <= rank_threshold;
                }
                return true;
            };

            if (ac.LE_mode == 0) {
                // process all elements
                for (size_t i = p_offset; i < p_c_size; ++i) {
                    // get old values
                    uint32_t c = pc_manager.get_c(i);
                    TypeSF si = pc_manager.get_si(i);
                    uint8_t accurate = 0;

                    // update and also update heap
                    s[s_size] = c;
                    si = t.evaluate_1D(s, s_size + 1) - score;
                    accurate = 1;

                    // insert into the manager
                    c_manager.add_entry(c, si, accurate);
                }
            } else {
                // initialize the heap
                for (size_t i = p_offset; i < p_offset + r; ++i) {
                    uint32_t c = pc_manager.get_c(i);
                    TypeSF si = pc_manager.get_si(i);
                    c_manager.SUB_heap_add(c, si);
                }

                // process first r candidates
                for (size_t i = p_offset; i < p_offset + r; ++i) {
                    // get old values
                    uint32_t c = pc_manager.get_c(i);
                    TypeSF si = pc_manager.get_si(i);
                    uint8_t accurate = 0;

                    // check if we can abort early
                    if (c_manager.get_SUB_heap_sum() <= r_score && c_manager.SUB_heap_min() >= si) {
                        return 0;
                    }

                    // check for update
                    bool update = update_scheme(i - p_offset, si);

                    if (update) {
                        // update and also update heap
                        s[s_size] = c;
                        si = t.evaluate_1D(s, s_size + 1) - score;
                        accurate = 1;

                        // update the heap
                        size_t idx = c_manager.SUB_heap_find(c);
                        c_manager.SUB_heap_update(idx, si);
                    }

                    // insert into the manager
                    c_manager.add_entry(c, si, accurate);
                }

                // process all other elements
                for (size_t i = p_offset + r; i < p_c_size; ++i) {
                    // get old values
                    uint32_t c = pc_manager.get_c(i);
                    TypeSF si = pc_manager.get_si(i);
                    uint8_t accurate = 0;

                    // check if we can abort early
                    if (c_manager.get_SUB_heap_sum() < r_score && c_manager.SUB_heap_min() >= si) {
                        return 0;
                    }

                    // check for update
                    bool update = update_scheme(i - p_offset, si);

                    if (update) {
                        // update and also update heap
                        s[s_size] = c;
                        si = t.evaluate_1D(s, s_size + 1) - score;
                        accurate = 1;
                    }

                    // update the heap, this will only insert if si is large enough
                    c_manager.SUB_heap_add(c, si);

                    // insert into the manager
                    c_manager.add_entry(c, si, accurate);
                }
            }

            c_manager.sort();
            c_manager.calc_csum();

            return c_manager.size;
        };

        /**
         * General RPC function that can be called by all heuristics.
         *
         * @param c_size Size of the candidates.
         * @param offset Offset into the candidates.
         * @param s_size Current size of set S.
         * @param score Current score of set S.
         * @param r1_si The score improvement of r - 1 elements.
         * @param largest_idx Largest Idx of the r-1 solution.
         * @param depth The current depth of the search tree.
         * @return The new number of candidates.
         */
        size_t general_RPC(const TypeSF r1_si, const size_t largest_idx, const size_t depth) {
            CandidateManager<TypeSF> &c_manager = c_managers[depth];
            const size_t c_size = c_manager.size;
            const size_t offset = c_manager.offset;
            const TypeSF score = scores[depth];
            size_t new_c_size = c_size;
            size_t r = k - s_size;

            while ((score + r1_si + c_manager.get_si(new_c_size - 1) <= best_score) && (new_c_size > largest_idx) && (offset + r < new_c_size)) {
                new_c_size -= 1;
            }

            return new_c_size;
        };

        /**
         * The 'Simple Upper Bound' heuristic. It computes the sum of the r
         * largest score improvements after the offset.
         *
         * @param c_size Size of the candidates.
         * @param offset Offset into the candidates.
         * @param s_size Current size of set S.
         * @param score Current score of set S.
         * @param depth The current depth of the search tree.
         * @return The bound.
         */
        TypeSF SUB(const size_t depth) {
            const CandidateManager<TypeSF> &c_manager = c_managers[depth];
            const size_t r = k - s_size;
            const TypeSF score = scores[depth];

            TypeSF bound = score + c_manager.get_partial_sum(c_manager.offset, r);

            return bound;
        };

        /**
         * The 'Reduction of Possible Candidates' heuristic. It checks if the sum
         * of the top (r-1) score improvements + the lowest score improvement is
         * smaller than the currently best achieved score. If so the last candidate
         * can be removed from the candidates.
         *
         * @param c_size Size of the candidates.
         * @param offset Offset into the candidates.
         * @param s_size Current size of set S.
         * @param score Current score of set S.
         * @param depth The current depth of the search tree.
         * @return The new number of candidates.
         */
        size_t RPC(const size_t depth) {
            const CandidateManager<TypeSF> &c_manager = c_managers[depth];
            const size_t c_size = c_manager.size;
            const size_t offset = c_manager.offset;
            const size_t r = k - s_size;
            const TypeSF score = scores[depth];

            size_t new_c_size = c_size;
            TypeSF top_r1_score = score + c_manager.get_partial_sum(offset, r - 1);

            while (top_r1_score + c_manager.get_si(new_c_size - 1) <= best_score && offset + r < new_c_size) {
                new_c_size -= 1;
            }

            return new_c_size;
        };

        /**
         * The 'Upper Bound 2D Heuristic'. It will calculate the pairwise score
         * improvements of the first l candidates and then choose a set of them so
         * r candidates are picked and the score is maximized. There are 2 versions
         * available. Either 'Greedy', the algorithm does not care that a candidate
         * will be picked multiple times or 'Matching', the algorithm fulfills the
         * constraint that each candidate can be picked at most one time. If we have
         * an uneven number of candidates to choose from the pairwise improvements,
         * we need to pick a solo candidate and its score improvement. There are two
         * ways how to accomplish this. Either we use option 'A' and pick the
         * largest score improvement first and then pick the pairs, or option 'B' we
         * first pick the pairs and after this the solo candidate with the greatest
         * score improvement remaining. Option 'AB' is also available and uses the
         * minimum of both options. In both variants the solo candidate is not
         * overlapping with any of the pairs. We need to use this procedure in l+1
         * cases. Case i is that i candidates are from the part of the pairs and
         * r-i candidates are from the part without the pairs. This is needed, since
         * we do not know the distribution of S^*.
         *
         * @param c_size Size of the candidates.
         * @param offset Offset into the candidates.
         * @param s_size Current size of set S.
         * @param score Current score of set S.
         * @param depth The current depth of the search tree.
         * @return The bound.
         */
        TypeSF UB2D(const size_t depth) {
            const CandidateManager<TypeSF> &c_manager = c_managers[depth];
            const size_t r = k - s_size;
            const TypeSF score = scores[depth];
            const size_t n_remaining = c_manager.size - c_manager.offset;
            const bool need_candidates = ac.UB2D_safe_skip_enabled;

            size_t l = std::min(ac.determine_UB2D_l(n, k, n_remaining, r), n_remaining);
            if (l < 2) {
                return std::numeric_limits<TypeSF>::max();
            }

            auto &ub2d_manager = ub2d_managers[depth];
            auto &si_cache = si_caches[depth];
            auto &ub2d_algorithm = ub2d_algorithms[depth];
            ub2d_algorithm.reinitialize(l, ac.UB2D_alg_type);

            // lazy skip
            if (ub2d_manager.do_lazy_skip()) {
                return std::numeric_limits<TypeSF>::max();
            }

            // safe skip
            if (ub2d_manager.safe_skip_enabled) {
                TypeSF inaccurate_bound = score + ub2d_manager.get_updated_r_si_bound(c_manager, si_cache);
                if (inaccurate_bound > best_score) {
                    return inaccurate_bound;
                }
            }
            ub2d_manager.clear();

            // Initialize the pairwise score improvements
            for (size_t i = 0; i < l; ++i) {
                size_t idx_1 = c_manager.offset + i;
                s[s_size] = c_manager.get_c(idx_1);

                for (size_t j = i + 1; j < l; ++j) {
                    size_t idx_2 = c_manager.offset + j;

                    // look in cache for value
                    size_t hash = si_cache.hash_2D(idx_1, idx_2);
                    TypeSF score_imp = si_cache.get_entry_2D(hash, idx_1, idx_2);
                    if (score_imp < 0) {
                        s[s_size + 1] = c_manager.get_c(idx_2);
                        score_imp = t.evaluate_2D(s, s_size + 2) - score;
                        si_cache.insert_entry_2D(hash, idx_1, idx_2, score_imp);
                    }
                    ub2d_algorithm.add_edge(i, j, score_imp);
                }
            }

            // calculate the upper bound
            TypeSF ub2d_r_si_bound = -std::numeric_limits<TypeSF>::max();
            TypeSF ub2d_r1_si_bound = -std::numeric_limits<TypeSF>::max();

            // clear all needed vectors
            p1_sets_a.clear();
            p1_sets_a.reserve(std::min(l, r) / 2);
            p1_sets_b.clear();
            p1_sets_b.reserve(std::min(l, r) / 2);
            p1_sets_even.clear();
            p1_sets_even.reserve(std::min(l, r) / 2);
            p1_si_a.clear();
            p1_si_a.reserve(std::min(l, r) / 2);
            p1_si_b.clear();
            p1_si_b.reserve(std::min(l, r) / 2);
            p1_si_even.clear();
            p1_si_even.reserve(std::min(l, r) / 2);
            VectorOfVectors<size_t> &p1_sets = p1_sets_a;
            std::vector<TypeSF> &p1_si = p1_si_a;

            // check all cases
            for (size_t i = l + r - std::min(n_remaining, l + r); i <= std::min(l, r); ++i) {
                size_t n_candidates_p1 = i;
                size_t n_candidates_p2 = r - i;

                TypeSF sum_even, sum_p1, sum_p2;

                if (n_candidates_p1 & 1) {
                    // odd case
                    size_t n_endpoints = n_candidates_p1 - 1;
                    TypeSF sum_a = std::numeric_limits<TypeSF>::max();
                    TypeSF sum_b = std::numeric_limits<TypeSF>::max();

                    // option A, top 1 and then pairs
                    if (ac.UB2D_odd_type == 1 || ac.UB2D_odd_type == 3) {
                        if (p1_sets_even.size != 0 && !p1_sets_even.contains(c_manager.offset)) {
                            p1_sets_a.copy_from(p1_sets_even);
                            copy(p1_si_a, p1_si_even);
                            sum_a = sum_even;
                        } else {
                            p1_sets_a.clear();
                            p1_si_a.clear();

                            sum_a = ub2d_algorithm.get_upper_bound_skip_first(p1_sets_a, p1_si_a, n_endpoints, c_manager.offset);
                        }
                        p1_sets_a.push_back(c_manager.offset);
                        p1_si_a.push_back(c_manager.get_si(c_manager.offset));
                        sum_a += c_manager.get_si(c_manager.offset);
                    }

                    // option B, pairs and then top 1
                    if (ac.UB2D_odd_type == 2 || ac.UB2D_odd_type == 3) {
                        if (p1_sets_even.size != 0) {
                            p1_sets_b.copy_from(p1_sets_even);
                            copy(p1_si_b, p1_si_even);
                            sum_b = sum_even;
                        } else {
                            p1_sets_b.clear();
                            p1_si_b.clear();

                            sum_b = ub2d_algorithm.get_upper_bound(p1_sets_b, p1_si_b, n_endpoints, c_manager.offset);
                        }
                        for (size_t j = c_manager.offset; j < c_manager.offset + l; ++j) {
                            if (!p1_sets_b.contains(j)) {
                                p1_sets_b.push_back(j);
                                p1_si_b.push_back(c_manager.get_si(j));
                                sum_b += c_manager.get_si(j);
                                break;
                            }
                        }
                    }

                    // update the reference
                    if (sum_a < sum_b) {
                        p1_sets = p1_sets_a;
                        p1_si = p1_si_a;
                        sum_p1 = sum_a;
                    } else {
                        p1_sets = p1_sets_b;
                        p1_si = p1_si_b;
                        sum_p1 = sum_b;
                    }
                } else {
                    // even case
                    p1_sets_even.clear();
                    p1_si_even.clear();

                    size_t n_endpoints = n_candidates_p1;
                    sum_even = ub2d_algorithm.get_upper_bound(p1_sets_even, p1_si_even, n_endpoints, c_manager.offset);

                    p1_sets = p1_sets_even;
                    p1_si = p1_si_even;
                    sum_p1 = sum_even;
                }

                // add top r - i from part2
                sum_p2 = c_manager.get_partial_sum(c_manager.offset + l, n_candidates_p2);

                // check if a new r solution was found
                if (ub2d_r_si_bound < sum_p1 + sum_p2) {
                    ub2d_r_si_bound = sum_p1 + sum_p2;
                    if (need_candidates) {
                        ub2d_manager.p1_sets.copy_from(p1_sets);
                        copy(ub2d_manager.p1_si, p1_si);

                        ub2d_manager.p2_c.resize(n_candidates_p2);
                        ub2d_manager.p2_si.resize(n_candidates_p2);
                        for (size_t j = 0; j < n_candidates_p2; ++j) {
                            ub2d_manager.p2_c[j] = c_manager.offset + l + j;
                            ub2d_manager.p2_si[j] = c_manager.get_si(c_manager.offset + l + j);
                        }

                        ub2d_manager.sum_p1 = sum_p1;
                        ub2d_manager.sum_p2 = sum_p2;
                    }
                }

                // check if a new r-1 solution was found
                if (n_candidates_p2 > 0) {
                    TypeSF sum_p2_r1 = c_manager.get_partial_sum(c_manager.offset + l, n_candidates_p2 - 1);
                    if (ub2d_r1_si_bound < sum_p1 + sum_p2_r1) {
                        ub2d_r1_si_bound = sum_p1 + sum_p2_r1;
                        if (need_candidates) {
                            ub2d_manager.p1_sets_r1.copy_from(p1_sets);
                            copy(ub2d_manager.p1_si_r1, p1_si);

                            ub2d_manager.p2_c_r1.resize(n_candidates_p2 - 1);
                            ub2d_manager.p2_si_r1.resize(n_candidates_p2 - 1);
                            for (size_t j = 0; j < n_candidates_p2 - 1; ++j) {
                                ub2d_manager.p2_c_r1[j] = c_manager.offset + l + j;
                                ub2d_manager.p2_si_r1[j] = c_manager.get_si(c_manager.offset + l + j);
                            }

                            ub2d_manager.sum_p1_r1 = sum_p1;
                            ub2d_manager.sum_p2_r1 = sum_p2_r1;
                        }
                    }
                }

                if (!ac.UB2D_RPC_enabled && best_score < score + ub2d_r_si_bound) {
                    // when no rpc, we can stop if we have found a true better solution
                    break;
                }

            }

            return score + ub2d_r_si_bound;
        };

        /**
         * The 'Partial Brute Force' heuristic.
         *
         * @param c_size Size of the candidates.
         * @param offset Offset into the candidates.
         * @param s_size Current size of set S.
         * @param score Current score of set S.
         * @param depth The current depth of the search tree.
         * @return The bound.
         */
        TypeSF PBF(const size_t depth) {
            const CandidateManager<TypeSF> &c_manager = c_managers[depth];
            const size_t r = k - s_size;
            const TypeSF score = scores[depth];
            const size_t n_remaining = c_manager.size - c_manager.offset;
            const bool need_candidates = ac.PBF_safe_skip_enabled;

            UpperBoundManager<TypeSF> &pbf_manager = pbf_managers[depth];
            SICache<TypeSF> &si_cache = si_caches[depth];

            // lazy skip
            if (pbf_manager.do_lazy_skip()) {
                return std::numeric_limits<TypeSF>::max();
            }

            // safe skip
            if (ac.PBF_safe_skip_enabled) {
                TypeSF inaccurate_bound = score + pbf_manager.get_updated_r_si_bound(c_manager, si_cache);
                if (inaccurate_bound > best_score) {
                    return inaccurate_bound;
                }
            }

            size_t n_b = ac.determine_PBF_n(n, k, n_remaining, r);
            size_t l_b = std::min(r, ac.determine_PBF_l(n, k, n_remaining, r)); // if r is smaller, there is no use to check larger sets

            while (n_b * l_b > n_remaining) {
                n_b -= 1;
            }
            if (n_b == 0) {
                return std::numeric_limits<TypeSF>::max();
            }

            const size_t l = n_b * l_b;

            PBFAlgorithm<TypeSF> &pbf_alg = pbf_algorithms[depth];
            pbf_alg.reinitialize(n_b, l_b, r, ac.PBF_alg_type, ac.PBF_safe_skip_enabled);

            // insert the 1d score improvements
            for (size_t b_id = 0; b_id < n_b; ++b_id) {
                for (size_t j = 0; j < l_b; ++j) {
                    size_t idx = c_manager.offset + (b_id * l_b) + j;
                    pbf_alg.add_set_1(b_id, idx, c_manager.get_si(idx));
                }
            }

            // insert the 2d score improvements
            for (size_t b_id = 0; b_id < n_b; ++b_id) {
                // Initialize the pairwise score improvements
                for (size_t i = 0; i < l_b; ++i) {
                    for (size_t j = i + 1; j < l_b; ++j) {
                        size_t idx_1 = c_manager.offset + (b_id * l_b) + i;
                        size_t idx_2 = c_manager.offset + (b_id * l_b) + j;

                        // look in cache for value
                        size_t hash = si_cache.hash_2D(idx_1, idx_2);
                        TypeSF score_imp = si_cache.get_entry_2D(hash, idx_1, idx_2);
                        if (score_imp < 0) {
                            s[s_size] = c_manager.get_c(idx_1);
                            s[s_size + 1] = c_manager.get_c(idx_2);
                            score_imp = t.evaluate_2D(s, s_size + 2) - score;
                            si_cache.insert_entry_2D(hash, idx_1, idx_2, score_imp);
                        }
                        pbf_alg.add_set_2(b_id, idx_1, idx_2, score_imp);
                    }
                }
            }

            // insert all other dimensions
            for (size_t b_id = 0; b_id < n_b; ++b_id) {
                size_t block_start = c_manager.offset + (b_id * l_b);
                // iterate over all dimensions
                for (size_t j = 3; j <= l_b; ++j) {
                    comb.resize(j);
                    std::iota(comb.begin(), comb.end(), 0);
                    comb[j - 1] -= 1;

                    set.resize(j);

                    // iterate over all possible sets
                    while (next_subset(comb, j, l_b)) {

                        // insert set
                        for (size_t m = 0; m < j; ++m) {
                            set[m] = block_start + comb[m];
                            s[s_size + m] = c_manager.get_c(block_start + comb[m]);
                        }

                        // look for score improvement in cache
                        size_t hash = si_cache.vector_hash(set);
                        TypeSF score_imp = si_cache.get_entry(hash, set);
                        if (score_imp < 0) {
                            score_imp = t.evaluate_XD(s, s_size + j) - score;
                            si_cache.insert_entry(hash, set, score_imp);
                        }
                        pbf_alg.add_set(b_id, set, j, score_imp);
                    }
                }
            }

            // calculate the upper bound
            TypeSF pbf_r_si_bound = -std::numeric_limits<TypeSF>::max();
            TypeSF pbf_r1_si_bound = -std::numeric_limits<TypeSF>::max();

            // choose from the pbf algorithm
            size_t loop_start = l + r - std::min(n_remaining, l + r);
            size_t loop_end = std::min(l, r) + 1;
            for (size_t i = loop_start; i < loop_end; ++i) {
                size_t n_candidates_p1 = i;
                size_t n_candidates_p2 = r - i;

                // get solution from p1
                p1_sets_even.clear();
                p1_si_even.clear();
                TypeSF sum_p1 = pbf_alg.get_solution(n_candidates_p1, p1_sets_even, p1_si_even);

                // get solution from p2
                TypeSF sum_p2 = c_manager.get_partial_sum(c_manager.offset + l, n_candidates_p2);

                // check if a new r solution was found
                if (pbf_r_si_bound < sum_p1 + sum_p2) {
                    pbf_r_si_bound = sum_p1 + sum_p2;
                    if (need_candidates) {
                        pbf_manager.p1_sets.copy_from(p1_sets_even);
                        copy(pbf_manager.p1_si, p1_si_even);

                        pbf_manager.p2_c.clear();
                        pbf_manager.p2_si.clear();
                        for (size_t j = 0; j < n_candidates_p2; ++j) {
                            pbf_manager.p2_c.push_back(c_manager.offset + l + j);
                            pbf_manager.p2_si.push_back(c_manager.get_si(c_manager.offset + l + j));
                        }

                        pbf_manager.sum_p1 = sum_p1;
                        pbf_manager.sum_p2 = sum_p2;
                    }
                }

                // check if a new r-1 solution was found
                if (n_candidates_p2 > 0) {
                    TypeSF sum_p2_r1 = c_manager.get_partial_sum(c_manager.offset + l, n_candidates_p2 - 1);
                    if (pbf_r1_si_bound < sum_p1 + sum_p2_r1) {
                        pbf_r1_si_bound = sum_p1 + sum_p2_r1;
                        if (need_candidates) {
                            pbf_manager.p1_sets_r1.copy_from(p1_sets_even);
                            copy(pbf_manager.p1_si_r1, p1_si_even);

                            pbf_manager.p2_c_r1.clear();
                            pbf_manager.p2_si_r1.clear();
                            for (size_t j = 0; j < n_candidates_p2 - 1; ++j) {
                                pbf_manager.p2_c_r1.push_back(c_manager.offset + l + j);
                                pbf_manager.p2_si_r1.push_back(c_manager.get_si(c_manager.offset + l + j));
                            }

                            pbf_manager.sum_p1_r1 = sum_p1;
                            pbf_manager.sum_p2_r1 = sum_p2_r1;
                        }
                    }
                }

                if (!ac.PBF_RPC_enabled && best_score < score + pbf_r_si_bound) {
                    // when no rpc, we can stop if we have found a true better solution
                    break;
                }

            }

            return score + pbf_r_si_bound;
        };

        /**
         * Recursive heuristic
         */
        TypeSF REC(const size_t depth) {
            const CandidateManager<TypeSF> &c_manager = c_managers[depth];
            const size_t r = k - s_size;
            const size_t n_remaining = c_manager.size - c_manager.offset;
            const TypeSF score = scores[depth];
            const TypeSF r_score = best_score - score;
            const bool need_candidates = ac.REC_safe_skip_enabled;

            auto rec_manager = rec_managers[depth];
            SICache<TypeSF> &si_cache = si_caches[depth];

            size_t l = std::min(ac.determine_REC_l(n, k, n_remaining, r), n_remaining);

            if (l == n_remaining) {
                return std::numeric_limits<TypeSF>::max();
            }

            // lazy skip
            if (rec_manager.do_lazy_skip()) {
                return std::numeric_limits<TypeSF>::max();
            }

            // safe skip
            if (ac.REC_safe_skip_enabled) {
                TypeSF inaccurate_bound = score + rec_manager.get_updated_r_si_bound(c_manager, si_cache);
                if (inaccurate_bound > best_score) {
                    return inaccurate_bound;
                }
            }

            rec_ac.partial_s.resize(s_size);
            for (size_t i = 0; i < s_size; ++i) {
                rec_ac.partial_s[i] = s[i];
            }

            rec_ac.init_candidates.resize(l);
            rec_ac.init_si.resize(l);
            rec_ac.init_acc.resize(l);
            for (size_t i = c_manager.offset; i < c_manager.offset + l; ++i) {
                rec_ac.init_candidates[i - c_manager.offset] = c_manager.get_c(i);
                rec_ac.init_si[i - c_manager.offset] = c_manager.get_si(i);
                rec_ac.init_acc[i - c_manager.offset] = c_manager.get_acc(i);
            }

            TypeSF rec_r_si_bound = -std::numeric_limits<TypeSF>::max();
            TypeSF rec_r1_si_bound = -std::numeric_limits<TypeSF>::max();

            p1_sets_a.clear();
            p1_sets_a.reserve(1);
            p1_si_a.clear();
            p1_si_a.reserve(1);
            VectorOfVectors<size_t> &p1_sets = p1_sets_a;
            std::vector<TypeSF> &p1_si = p1_si_a;

            TypeSF sum_p1, sum_p2;
            std::vector<uint32_t> rec_sol;
            std::vector<size_t> rec_sol_index;

            for (size_t i = 0; i <= r; ++i) {
                size_t n_candidates_p1 = i;
                size_t n_candidates_p2 = r - i;

                if (n_candidates_p1 > l || n_candidates_p2 > n_remaining - l) {
                    continue;
                }

                sum_p2 = c_manager.get_partial_sum(c_manager.offset + l, n_candidates_p2);

                rec_sol_index.clear();

                if (n_candidates_p1 == 0) {
                    sum_p1 = 0;
                } else if (n_candidates_p1 == 1) {
                    sum_p1 = c_manager.get_si(c_manager.offset);
                    rec_sol.resize(1);
                    rec_sol[0] = c_manager.get_c(c_manager.offset);
                } else {
                    rec_ac.k = s_size + n_candidates_p1;
                    rec_ac.score_threshold = score + (r_score - sum_p2);
                    rec_ac.finalize();
                    auto tree_search = TreeSearch<T, TypeSF>(t, rec_ac);
                    auto res = tree_search.search();

                    sum_p1 = tree_search.best_score - score;
                    rec_sol = tree_search.best_s;
                }

                if (rec_r_si_bound < sum_p1 + sum_p2) {
                    rec_r_si_bound = sum_p1 + sum_p2;
                    if (need_candidates) {
                        // fill the index
                        for (auto c: rec_sol) {
                            for (size_t j = c_manager.offset; j < c_manager.offset + l; ++j) {
                                if (c == c_manager.get_c(j)) {
                                    rec_sol_index.push_back(j);
                                    break;
                                }
                            }
                        }

                        p1_sets.push_back(rec_sol_index);
                        p1_si.push_back(sum_p1);

                        rec_manager.p1_sets.copy_from(p1_sets);
                        copy(rec_manager.p1_si, p1_si);

                        rec_manager.p2_c.resize(n_candidates_p2);
                        rec_manager.p2_si.resize(n_candidates_p2);
                        for (size_t j = 0; j < n_candidates_p2; ++j) {
                            rec_manager.p2_c[j] = c_manager.offset + l + j;
                            rec_manager.p2_si[j] = c_manager.get_si(c_manager.offset + l + j);
                        }

                        rec_manager.sum_p1 = sum_p1;
                        rec_manager.sum_p2 = sum_p2;
                    }
                }

                if (n_candidates_p2 > 0) {
                    TypeSF sum_p2_r1 = c_manager.get_partial_sum(c_manager.offset + l, n_candidates_p2 - 1);
                    if (rec_r1_si_bound < sum_p1 + sum_p2_r1) {
                        rec_r1_si_bound = sum_p1 + sum_p2_r1;
                        if (need_candidates) {

                            if (rec_sol_index.empty()) {
                                // fill the index
                                for (auto c: rec_sol) {
                                    for (size_t j = c_manager.offset; j < c_manager.offset + l; ++j) {
                                        if (c == c_manager.get_c(j)) {
                                            rec_sol_index.push_back(j);
                                            break;
                                        }
                                    }
                                }
                            }

                            p1_sets.push_back(rec_sol_index);
                            p1_si.push_back(sum_p1);

                            rec_manager.p1_sets.copy_from(p1_sets);
                            copy(rec_manager.p1_si, p1_si);

                            rec_manager.p2_c.resize(n_candidates_p2 - 1);
                            rec_manager.p2_si.resize(n_candidates_p2 - 1);
                            for (size_t j = 0; j < n_candidates_p2 - 1; ++j) {
                                rec_manager.p2_c[j] = c_manager.offset + l + j;
                                rec_manager.p2_si[j] = c_manager.get_si(c_manager.offset + l + j);
                            }

                            rec_manager.sum_p1 = sum_p1;
                            rec_manager.sum_p2 = sum_p2;
                        }
                    }
                }

                if (!ac.REC_RPC_enabled && rec_r_si_bound + score > best_score) {
                    break;
                }
            }
            return rec_r_si_bound + score;
        }

        /**
         * Check if the execution time has exceeded the specified time limit.
         *
         * @return `true` if the execution time has exceeded the time limit, otherwise `false`.
         */
        bool has_time_exceeded() {
            if (n_tried_calls > 100) {
                double seconds = get_elapsed_seconds(time_sp, get_time_point());
                time_exceeded = seconds > ac.time_limit;
                n_tried_calls = 0;
                return time_exceeded;
            }
            n_tried_calls += 1;
            return false;
        };

        /**
         * Serialize object state to a JSON-formatted string.
         *
         * This function serializes the state of an object into a JSON-formatted string. It constructs
         * a JSON object with key-value pairs representing the object's attributes. Each key is a string
         * representing an attribute name, and the corresponding value is the JSON representation of
         * the attribute. The function returns the resulting JSON string.
         *
         * @return A JSON-formatted string representing the serialized object state.
         */
        std::string parse_to_JSON() {
            std::string content = "{\n";
            content += "\"n\" : " + to_JSON_value(n) + ",\n";
            content += "\"k\" : " + to_JSON_value(k) + ",\n";
            content += "\"initial_s\" : " + to_JSON(initial_s) + ",\n";
            content += "\"initial_score\" : " + to_JSON_value(initial_score) + ",\n";
            content += "\"best_greedy_s\" : " + to_JSON(best_greedy_s) + ",\n";
            content += "\"best_greedy_score\" : " + to_JSON_value(best_greedy_score) + ",\n";
            content += "\"best_s\" : " + to_JSON(best_s) + ",\n";
            content += "\"best_score\" : " + to_JSON_value(best_score) + ",\n";
            content += "\"time_needed\" : " + to_JSON_value(time_needed) + ",\n";
            content += "\"time_limit_exceeded\" : " + to_JSON_value(time_exceeded) + ",\n";
            content += "\"UB2D_success_time\" : " + to_JSON_value(UB2D_success_time) + ",\n";
            content += "\"UB2D_failure_time\" : " + to_JSON_value(UB2D_failure_time) + ",\n";
            content += "\"UB2D_n_success\" : " + to_JSON_value(UB2D_n_success) + ",\n";
            content += "\"UB2D_n_failure\" : " + to_JSON_value(UB2D_n_failure) + ",\n";
            content += "\"PBF_success_time\" : " + to_JSON_value(PBF_success_time) + ",\n";
            content += "\"PBF_failure_time\" : " + to_JSON_value(PBF_failure_time) + ",\n";
            content += "\"PBF_n_success\" : " + to_JSON_value(PBF_n_success) + ",\n";
            content += "\"PBF_n_failure\" : " + to_JSON_value(PBF_n_failure) + ",\n";
            content += "\"REC_success_time\" : " + to_JSON_value(REC_success_time) + ",\n";
            content += "\"REC_failure_time\" : " + to_JSON_value(REC_failure_time) + ",\n";
            content += "\"REC_n_success\" : " + to_JSON_value(REC_n_success) + ",\n";
            content += "\"REC_n_failure\" : " + to_JSON_value(REC_n_failure) + ",\n";
            content += "\"program_options\" : {\n" + ac.to_JSON() + "\n}\n";
            content += "}";
            return content;
        };

        /**
         * Write serialized object state to a JSON file.
         *
         * This function serializes the state of an object into a JSON-formatted string using the
         * `parse_to_JSON` function and writes it to a JSON file specified by the `output_file_path`
         * attribute of the object. The existing content of the file, if any, will be overwritten.
         * If the file does not exist, it will be created.
         */
        void write_to_JSON() {
            std::string content = parse_to_JSON();
            std::ofstream file;
            file.open(ac.output_file_path);
            file << content;
            file.close();
        };
    };

}

#endif //SUBSETOPTIMIZATION_TREESEARCH_H
