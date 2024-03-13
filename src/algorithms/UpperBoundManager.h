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

#ifndef SUBSETOPTIMIZATION_UPPERBOUNDMANAGER_H
#define SUBSETOPTIMIZATION_UPPERBOUNDMANAGER_H

#include <cstddef>

#include "../utility/Util.h"
#include "SICache.h"
#include "../utility/VectorOfVectors.h"
#include "CandidateManager.h"

namespace SMSM {

/**
 * Class that manages the results from heuristics.
 */
    template<typename TypeSF>
    class UpperBoundManager {
    public:
        // solution for r
        VectorOfVectors<size_t> p1_sets;
        std::vector<TypeSF> p1_si;
        std::vector<size_t> p2_c;
        std::vector<TypeSF> p2_si;
        TypeSF sum_p1;
        TypeSF sum_p2;

        // solution for r-1
        VectorOfVectors<size_t> p1_sets_r1;
        std::vector<TypeSF> p1_si_r1;
        std::vector<size_t> p2_c_r1;
        std::vector<TypeSF> p2_si_r1;
        TypeSF sum_p1_r1;
        TypeSF sum_p2_r1;

        std::vector<size_t> valid_r_configurations;

        // values for lazy and safe skipping
        double lazy_skip_value;
        double lazy_skip_start_value;
        double lazy_skip_add;
        bool safe_skip_enabled;

        /**
        * Constructs an `UpperBoundManager` instance with initial skip settings.
        *
        * @param init_skip_value The initial skip value.
        * @param init_skip_value_add The value to be added to lazy_skip_value.
        * @param init_safe_skip_enabled A flag indicating whether safe skip is enabled.
        */
        inline UpperBoundManager(double init_skip_value, double init_skip_value_add, bool init_safe_skip_enabled) {
            lazy_skip_value = init_skip_value;
            lazy_skip_start_value = init_skip_value;
            lazy_skip_add = init_skip_value_add;
            safe_skip_enabled = init_safe_skip_enabled;
        }

        /**
         * Reinitializes the UpperBoundManager.
         *
         * @param init_skip_value The initial skip value.
         * @param init_skip_value_add The value to be added to lazy_skip_value.
         * @param init_safe_skip_enabled A flag indicating whether safe skip is enabled.
         */
        inline void reinitialize(double init_skip_value, double init_skip_value_add, bool init_safe_skip_enabled) {
            lazy_skip_value = init_skip_value;
            lazy_skip_add = init_skip_value_add;
            safe_skip_enabled = init_safe_skip_enabled;
        }

        /**
         * Perform lazy skip by incrementing the skip value.
         *
         * @return `true` if we can lazy skip, `false` else.
         */
        inline bool do_lazy_skip() {
            lazy_skip_value += lazy_skip_add;
            if (lazy_skip_value >= 1.0) {
                lazy_skip_value -= 1.0;
                return false;
            }
            return true;
        }

        /**
         * Clears all stored results in the manager.
         */
        inline void clear() {
            p1_sets.clear();
            p1_si.clear();
            p2_c.clear();
            p2_si.clear();

            p1_sets_r1.clear();
            p1_si_r1.clear();
            p2_c_r1.clear();
            p2_si_r1.clear();
        }

        /**
         * Checks if the manager has a solution for r elements.
         *
         * @return `true` if it has a solution, `false` otherwise.
         */
        inline bool has_r_solution() const {
            return (p1_si.size() + p2_si.size()) > 0;
        }

        /**
         * Checks if the manager has a solution for r-1 elements.
         *
         * @return `true` if it has a solution, `false` otherwise.
         */
        inline bool has_r1_solution() const {
            return (p1_sets_r1.size + p2_c_r1.size()) > 0;
        }

        /**
         * Gets the current (inaccurate) si bound for r elements.
         *
         * @return The sum of the score improvements.
         */
        inline TypeSF get_r_si_bound() const {
            return sum_p1 + sum_p2;
        }

        /**
         * Gets the current (inaccurate) si bound for r-1 elements.
         *
         * @return The sum of the score improvements.
         */
        inline TypeSF get_r1_si_bound() const {
            return sum_p1_r1 + sum_p2_r1;
        }

        /**
         * Determines the largest index in the r solution.
         *
         * @return The largest index or `std::numeric_limits<size_t>::max()` if no solution exists.
         */
        inline size_t get_largest_idx_of_r_solution() {
            if (!has_r_solution()) {
                return std::numeric_limits<size_t>::max();
            }

            size_t max_idx = 0;
            if (!p2_c.empty()) {
                for (auto idx: p2_c) {
                    max_idx = std::max(max_idx, idx);
                }
                return max_idx;
            }

            for (size_t i = 0; i < p1_sets.size; ++i) {
                std::vector<size_t> &vec = p1_sets.vec[i];
                for (auto idx: vec) {
                    max_idx = std::max(max_idx, idx);
                }
            }

            return max_idx;
        }

        /**
         * Determines the largest index in the r solution.
         *
         * @return The largest index or `std::numeric_limits<size_t>::max()` if no solution exists.
         */
        inline size_t get_largest_idx_of_r1_solution() {
            if (!has_r_solution()) {
                return std::numeric_limits<size_t>::max();
            }

            size_t max_idx = 0;
            if (!p2_c_r1.empty()) {
                for (auto idx: p2_c_r1) {
                    max_idx = std::max(max_idx, idx);
                }
                return max_idx;
            }

            for (size_t i = 0; i < p1_sets_r1.size; ++i) {
                std::vector<size_t> &vec = p1_sets_r1.vec[i];
                for (auto idx: vec) {
                    max_idx = std::max(max_idx, idx);
                }
            }

            return max_idx;
        }

        /**
         * Calculates the updated si bound for r elements.
         *
         * @param c_size The size of the candidate set.
         * @param offset The offset for indexing candidates.
         * @param c1d_manager The candidate manager.
         * @param si_cache The SI cache.
         * @return The updated si bound.
         */
        inline TypeSF get_updated_r_si_bound(const CandidateManager<TypeSF> &c1d_manager,
                                             SICache<TypeSF> &si_cache) {
            TypeSF sum = 0.0;
            for (size_t j = 0; j < p1_sets.size; ++j) {
                std::vector<size_t> &set = p1_sets.vec[j];
                size_t original_size = set.size();
                for (size_t i = 0; i < set.size(); ++i) {
                    size_t idx = set.size() - i - 1;
                    if (set[idx] < c1d_manager.offset || set[idx] >= c1d_manager.size) {
                        std::swap(set[idx], set[set.size() - 1]);
                        set.pop_back();
                    }
                }

                if (original_size == set.size()) {
                    sum += p1_si[j];
                } else if (set.size() == 1) {
                    sum += c1d_manager.get_si(set[0]);
                } else if (set.size() == 2) {
                    std::sort(set.begin(), set.end());
                    size_t hash = si_cache.hash_2D(set[0], set[1]);
                    TypeSF score_imp = si_cache.get_entry_2D(hash, set[0], set[1]);
                    if (score_imp >= 0) {
                        sum += score_imp;
                        // we do not calculate the improvement if we don't know the set
                    }
                } else if (set.size() > 2) {
                    std::sort(set.begin(), set.end());
                    size_t hash = si_cache.vector_hash(set);
                    TypeSF score_imp = si_cache.get_entry(hash, set);
                    if (score_imp >= 0) {
                        sum += score_imp;
                        // we do not calculate the improvement if we don't know the set
                    }
                }
            }
            return sum;
        }

        /**
         * Notifies the manager of visiting a new depth, clearing previous results.
         */
        inline void visit_new_depth() {
            lazy_skip_value = lazy_skip_start_value;
            clear();
        }
    };

}

#endif //SUBSETOPTIMIZATION_UPPERBOUNDMANAGER_H
