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

#ifndef SUBSETOPTIMIZATION_CANDIDATEMANAGER_H
#define SUBSETOPTIMIZATION_CANDIDATEMANAGER_H

#include <cstdint>
#include <algorithm>

#include "../utility/Util.h"

namespace SMSM {

/**
 * Represents a candidate with score improvement and accuracy information.
 *
 * @tparam TypeSF The type of the 'score_improvement' attribute (e.g., int, float, double).
 */
    template<typename TypeSF>
    struct CandidateEntry {
        TypeSF score_improvement;
        uint32_t candidate;
        uint8_t accurate;
    };

/**
 * Manages a collection of candidate entries for optimization.
 *
 * @tparam TypeSF The numeric type used for score improvement.
 */
    template<typename TypeSF>
    class CandidateManager {
    private:
        std::vector<CandidateEntry<TypeSF>> candidates;
        std::vector<TypeSF> csum_si;

        std::vector<CandidateEntry<TypeSF>> SUB_heap; // Simple Upper-Bound Heap
        size_t SUB_heap_size;
        size_t SUB_heap_max_size;
        TypeSF SUB_heap_sum;

    public:
        size_t offset; // current offset into the candidates
        size_t size; // current size of the candidates
        size_t max_size; // maximum size of the candidates

        /**
         * Constructor.
         *
         * @param n Number of candidates.
         */
        explicit CandidateManager(size_t n) {
            candidates.resize(n);
            csum_si.resize(n + 1);

            offset = 0;
            size = 0;
            max_size = n;

            SUB_heap_size = 0;
            SUB_heap_max_size = 0;
            SUB_heap_sum = 0;
        };

        /**
         * Reinitializes the CandidateManager.
         *
         * @param n New maximum number of candidates.
         */
        void reinitialize(size_t n) {
            candidates.resize(n);
            csum_si.resize(n + 1);

            offset = 0;
            size = 0;
            max_size = n;

            SUB_heap_size = 0;
            SUB_heap_max_size = 0;
            SUB_heap_sum = 0;
        }

        /**
         * Fills the candidate array with values 0, ..., n-1.
         *
         * @param n Number of elements.
         */
        void fill_candidates(size_t n,
                             const std::vector<uint32_t> &s_part,
                             const std::vector<uint32_t> &init_candidates,
                             const std::vector<double> &init_si,
                             const std::vector<uint8_t> &init_acc) {
            size_t idx = 0;
            if (init_candidates.empty()) {
                for (uint32_t c = 0; c < n; ++c) {
                    if (!contains(s_part, c)) {
                        candidates[idx].candidate = c;
                        candidates[idx].accurate = 0;
                        idx += 1;
                    }
                }
            } else {
                if (init_candidates.size() == init_si.size()) {
                    for (size_t i = 0; i < init_candidates.size(); ++i) {
                        uint32_t c = init_candidates[i];
                        TypeSF si = init_si[i];
                        uint8_t acc = init_acc[i];
                        if (!contains(s_part, c)) {
                            candidates[idx].candidate = c;
                            candidates[idx].score_improvement = si;
                            candidates[idx].accurate = acc;
                            idx += 1;
                        }
                    }
                } else {
                    for (unsigned int c: init_candidates) {
                        if (!contains(s_part, c)) {
                            candidates[idx].candidate = c;
                            candidates[idx].accurate = 0;
                            idx += 1;
                        }
                    }
                }
            }
            offset = 0;
            size = idx;
        };

        /**
         * Fills the CandidateManager, based on another CandidateManager.
         *
         * @param ref_manager The other CandidateManager.
         */
        void fill_from_candidate_manager(const CandidateManager<TypeSF> &ref_manager) {
            for (size_t i = ref_manager.offset; i < ref_manager.size; ++i) {
                candidates[i - ref_manager.offset].candidate = ref_manager.get_c(i);
            }
            offset = 0;
            size = ref_manager.size - ref_manager.offset;
        }

        /**
        * Add a candidate entry to the manager.
        *
        * @param candidate_id An identifier or index for the candidate.
        * @param score_improvement The score improvement.
        * @param is_accurate A flag indicating the accuracy.
        */
        void add_entry(uint32_t candidate_id, TypeSF score_improvement, uint8_t is_accurate) {
            candidates[size].candidate = candidate_id;
            candidates[size].score_improvement = score_improvement;
            candidates[size].accurate = is_accurate;
            size += 1;
        }

        /**
         * Sets an entry in the CandidateManager.
         *
         * @param idx Index of the entry.
         * @param candidate_id New id.
         * @param score_improvement New score improvement.
         * @param is_accurate New accuracy.
         */
        void set_entry(size_t idx, uint32_t candidate_id, TypeSF score_improvement, uint8_t is_accurate) {
            candidates[idx].candidate = candidate_id;
            candidates[idx].score_improvement = score_improvement;
            candidates[idx].accurate = is_accurate;
        }

        /**
        * Sort the candidates by score improvement in descending order.
        */
        void sort() {
            auto first = candidates.begin();

            auto last = candidates.begin();
            std::advance(last, size);

            std::sort(first, last, [](const CandidateEntry<TypeSF> &a, const CandidateEntry<TypeSF> &b) {
                return a.score_improvement > b.score_improvement;
            });
        }

        /**
        * Calculate the cumulative sum of score improvements for candidates.
        */
        void calc_csum() {
            TypeSF cumulative_sum = 0;
            csum_si[0] = 0; // Initialize the first element.

            for (size_t i = 0; i < size; ++i) {
                cumulative_sum += candidates[i].score_improvement;
                csum_si[i + 1] = cumulative_sum;
            }
        }

        /**
         * Clears the manager.
         */
        void clear(size_t r) {
            SUB_heap.resize(r);
            SUB_heap_size = 0;
            SUB_heap_max_size = r;
            SUB_heap_sum = 0;

            size = 0;
        }

        /**
         * Adds a candidate onto the heap.
         *
         * @param candidate The candidate.
         * @param score_improvement The score improvement.
         */
        void SUB_heap_add(uint32_t candidate, TypeSF score_improvement) {
            if (SUB_heap_size < SUB_heap_max_size) {
                // add the element to the heap if it has less than r elements
                SUB_heap_sum += score_improvement;
                SUB_heap[SUB_heap_size].candidate = candidate;
                SUB_heap[SUB_heap_size].score_improvement = score_improvement;
                SUB_heap_size += 1;

                if (SUB_heap_size == SUB_heap_max_size) {
                    // if the heap is full, then heapify
                    std::make_heap(SUB_heap.begin(), SUB_heap.end(), [](const CandidateEntry<TypeSF> &a, const CandidateEntry<TypeSF> &b) { return a.score_improvement > b.score_improvement; });
                }
            } else {
                if (score_improvement > SUB_heap[0].score_improvement) {
                    // we have to push the score improvement on the heap
                    SUB_heap_sum = SUB_heap_sum - SUB_heap[0].score_improvement + score_improvement;
                    SUB_heap[0].candidate = candidate;
                    SUB_heap[0].score_improvement = score_improvement;

                    std::make_heap(SUB_heap.begin(), SUB_heap.end(), [](const CandidateEntry<TypeSF> &a, const CandidateEntry<TypeSF> &b) { return a.score_improvement > b.score_improvement; });
                }
            }
        }

        /**
         * Finds a candidate in the heap and returns its index.
         *
         * @param candidate The candidate.
         * @return The index of the candidate.
         */
        size_t SUB_heap_find(uint32_t candidate) {
            for (size_t i = 0; i < SUB_heap_size; ++i) {
                if (SUB_heap[i].candidate == candidate) {
                    return i;
                }
            }
            std::cout << "Error: Could not find candidate " << candidate << std::endl;
            abort();
        }

        /**
         * Updates the score improvement at the current index of the heap.
         *
         * @param idx The index.
         * @param score_improvement The new score improvement.
         */
        void SUB_heap_update(size_t idx, TypeSF score_improvement) {
            SUB_heap_sum = SUB_heap_sum - SUB_heap[idx].score_improvement + score_improvement;

            SUB_heap[idx].score_improvement = score_improvement;
            std::make_heap(SUB_heap.begin(), SUB_heap.end(), [](const CandidateEntry<TypeSF> &a, const CandidateEntry<TypeSF> &b) { return a.score_improvement > b.score_improvement; });
        }

        /**
         * Returns the minimum score improvement the heap.
         *
         * @return The score improvement.
         */
        TypeSF SUB_heap_min() {
            return SUB_heap[0].score_improvement;
        }

        /**
         * Returns the sum of the heap.
         *
         * @return Sum of the heap.
         */
        TypeSF get_SUB_heap_sum() {
            return SUB_heap_sum;
        }

        /**
        * Calculate the sum of score improvements for a range of candidates.
        *
        * @param start_idx The index of the first candidate in the range.
        * @param count The number of candidates to include in the sum.
        * @return The sum of score improvements for the specified range of candidates.
        */
        TypeSF get_partial_sum(size_t start_idx, size_t count) const {
            // Calculate the sum of score improvements for the specified range.
            return csum_si[start_idx + count] - csum_si[start_idx];
        }

        /**
         * Returns the score improvement of the candidate at the specified index.
         *
         * @param idx The index.
         * @return The score improvement.
         */
        TypeSF get_si(size_t idx) const {
            return candidates[idx].score_improvement;
        }

        /**
         * Returns the candidate at the specified index.
         *
         * @param idx The index.
         * @return The candidate.
         */
        uint32_t get_c(size_t idx) const {
            return candidates[idx].candidate;
        }

        /**
         * Returns the accuracy of score improvement of the candidate at the specified index.
         *
         * @param idx The index.
         * @return The accuracy.
         */
        uint8_t get_acc(size_t idx) const {
            return candidates[idx].accurate;
        }
    };

}

#endif //SUBSETOPTIMIZATION_CANDIDATEMANAGER_H
