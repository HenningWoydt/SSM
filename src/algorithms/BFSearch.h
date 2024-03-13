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

#ifndef SUBSETOPTIMIZATION_BFSEARCH_H
#define SUBSETOPTIMIZATION_BFSEARCH_H

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <vector>
#include <cstring>

#include "../utility/Util.h"

namespace SMSM {

    template<typename T, typename TypeSF>
    class BFSearch {
    private:
        T &t;
        size_t n;
        size_t k;
        size_t sf_evaluated = 0;

        std::vector<uint32_t> best_s;
        std::vector<uint32_t> s;

        TypeSF best_score;

    public:
        /**
         * Constructor for the Brute Force Search algorithm.
         *
         * @param t The data structure to be optimized.
         * @param k Number of elements in the solution.
         */
        BFSearch(T &t,
                 size_t k) : t(t), k(k) {
            n = t.get_n();

            best_s.resize(k);
            s.resize(k);

            best_score = -std::numeric_limits<TypeSF>::max();
        }

        /**
         * Reinitializes the Brute-Force Search algorithm.
         *
         * @param t_new The new data structure.
         * @param k_new The new number of elements.
         */
        void reinitialize(T &t_new,
                          size_t k_new) {
            t = t_new;
            n = t.get_n();
            k = k_new;

            best_s.resize(k);
            s.resize(k);

            best_score = -std::numeric_limits<TypeSF>::max();
        }

        /**
         * Starts the search.
         *
         * @return The solution.
         */
        std::vector<uint32_t> search() {
            // initial greedy solution
            simple_greedy();
            best_score = -std::numeric_limits<TypeSF>::max();

            // special case at the top of the tree
            for (uint32_t c = 0; c < n; ++c) {
                s[0] = c;
                recursive_search(1);
            }
            return best_s;
        };

    private:
        /**
         * Recursively searches for a solution in a tree-like manner.
         *
         * @param depth Current depth of the search algorithm.
         */
        void recursive_search(size_t depth) {
            if (depth == k) {
                // if the size of S is k then evaluate
                sf_evaluated += 1;
                TypeSF score = t.evaluate_general(s, depth);
                if (score > best_score) {
                    best_score = score;
                    for (size_t i = 0; i < k; ++i) {
                        best_s[i] = s[i];
                    }
                }
            } else {
                // determine which element to add to the set
                uint32_t max_c = s[depth - 1];
                for (uint32_t c = max_c + 1; c < n; ++c) {
                    s[depth] = c;
                    recursive_search(depth + 1);
                }
            }
        }

        /**
         * Simple greedy algorithm to compute an initial solution.
         */
        void simple_greedy() {
            TypeSF score = -std::numeric_limits<TypeSF>::max();

            for (size_t i = 0; i < k; ++i) {
                uint32_t best_c = 0;

                for (uint32_t c = 0; c < n; ++c) {
                    if (!contains(best_s, i, c)) {
                        best_s[i] = c;
                        TypeSF new_score = t.evaluate_general(best_s, i);

                        if (new_score >= score) {
                            score = new_score;
                            best_c = c;
                        }
                    }
                }
                s[i] = best_c;
            }
        }
    };
}

#endif //SUBSETOPTIMIZATION_BFSEARCH_H
