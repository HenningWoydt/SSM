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

#ifndef SUBSETOPTIMIZATION_GRAPHNEGATIVEGROUPFARNESS_H
#define SUBSETOPTIMIZATION_GRAPHNEGATIVEGROUPFARNESS_H

#include "../utility/Util.h"
#include "Graph.h"

namespace SMSM {

/**
 * Graph structure to optimize for negative group farness.
 */
    template<typename TypeSF>
    class GraphNegativeGroupFarness final : public Graph<TypeSF> {
    public:
        using Graph<TypeSF>::Graph;
        size_t padded_n_nodes = round_up(Graph<TypeSF>::n_nodes, (size_t) 64);
        size_t extra_nodes = padded_n_nodes - Graph<TypeSF>::n_nodes;

        // distance matrix
        std::vector<std::vector<TypeSF>> dist_mtx;

        // structures to speed up score function evaluation
        size_t depth = 0;
        std::vector<std::vector<TypeSF>> min_dist;
        std::vector<TypeSF> temp_min;

        inline TypeSF evaluate_empty_set() override {
            return -(TypeSF) (Graph<TypeSF>::n_nodes * Graph<TypeSF>::n_nodes);
        };

        inline TypeSF evaluate_1D(const std::vector<uint32_t> &s, const size_t s_size) override {
            TypeSF score = sum_of_min(min_dist[depth], dist_mtx[s[s_size - 1]]);
            ASSERT(evaluate_general(s, s_size) == -score);
            return -score;
        };

        inline TypeSF evaluate_2D(const std::vector<uint32_t> &s, const size_t s_size) override {
            TypeSF score = sum_of_min(min_dist[depth], dist_mtx[s[s_size - 2]], dist_mtx[s[s_size - 1]]);
            ASSERT(evaluate_general(s, s_size) == -score);
            return -score;
        };

        inline TypeSF evaluate_XD(const std::vector<uint32_t> &s, size_t s_size) override {
            size_t n_new_elements = s_size - depth;

            min(temp_min, min_dist[depth], dist_mtx[s[depth]]);
            for (size_t j = 1; j < n_new_elements - 1; ++j) {
                min_in_place(temp_min, dist_mtx[s[depth + j]]);
            }
            TypeSF score = sum_of_min(temp_min, dist_mtx[s[depth + n_new_elements - 1]]);
            ASSERT(evaluate_general(s, s_size) == -score);
            return -score;
        };

        inline TypeSF evaluate_general(const std::vector<uint32_t> &s, size_t s_size) override {
            if (s_size == 0) {
                return -(TypeSF) (Graph<TypeSF>::n_nodes * Graph<TypeSF>::n_nodes);
            } else if (s_size == 1) {
                return -sum(dist_mtx[s[0]]);
            } else if (s_size == 2) {
                return -sum_of_min(dist_mtx[s[0]], dist_mtx[s[1]]);
            }

            min(temp_min, dist_mtx[s[0]], dist_mtx[s[1]]);
            for (size_t j = 2; j < s_size - 1; ++j) {
                min_in_place(temp_min, dist_mtx[s[j]]);
            }
            TypeSF score = sum_of_min(temp_min, dist_mtx[s[s_size - 1]]);
            return -score;
        };

        inline void finalize() override {
            Graph<TypeSF>::sort_unique_neighbours();
            temp_min.resize(padded_n_nodes);
            initialize_dist_mtx();

            Graph<TypeSF>::max_reachable_score = 0;
        };

        inline void initialize_helping_structures(size_t k) override {
            min_dist.clear();
            min_dist.resize((k + 1), std::vector<TypeSF>(padded_n_nodes, std::numeric_limits<TypeSF>::max()));

            for (size_t i = 0; i < k + 1; ++i) {
                for (size_t j = 0; j < extra_nodes; ++j) {
                    min_dist[i][padded_n_nodes - extra_nodes + j] = 0;
                }
            }

            depth = 0;
        };

        inline void visit_new_depth(const std::vector<uint32_t> &s, size_t s_size) override {
            depth += 1;

            min(min_dist[depth], min_dist[(depth - 1)], dist_mtx[s[s_size - 1]]);
        };

        inline void return_from_last_depth() override {
            depth -= 1;
        };

        /**
         * Initializes the distance matrix.
         */
        inline void initialize_dist_mtx() {
            dist_mtx.resize(Graph<TypeSF>::n_nodes, std::vector<TypeSF>(padded_n_nodes, 0));

            // initialize arrays to help
            std::vector<uint8_t> bool_arr(Graph<TypeSF>::n_nodes);
            std::vector<uint32_t> stack1(Graph<TypeSF>::n_nodes);
            std::vector<uint32_t> stack2(Graph<TypeSF>::n_nodes);

            for (uint32_t i = 0; i < Graph<TypeSF>::n_nodes; ++i) {
                // for every node
                int curr_distance = 0;
                size_t stack1_size = 0;
                size_t stack2_size = 0;
                std::fill(bool_arr.begin(), bool_arr.end(), 0);
                stack1[stack1_size++] = i;

                while (stack1_size != 0) {
                    // set distance for current stack
                    for (size_t j = 0; j < stack1_size; ++j) {
                        uint32_t node = stack1[j];
                        dist_mtx[i][node] = curr_distance;
                        bool_arr[node] = 1;
                    }

                    // get the next stack
                    for (size_t j = 0; j < stack1_size; ++j) {
                        uint32_t node = stack1[j];

                        for (uint32_t neighbour: Graph<TypeSF>::adj_list[node]) {
                            if (!bool_arr[neighbour]) {
                                bool_arr[neighbour] = 1;
                                stack2[stack2_size++] = neighbour;
                            }
                        }
                    }
                    curr_distance++;

                    stack1.swap(stack2);
                    std::swap(stack1_size, stack2_size);
                    stack2_size = 0;
                }
            }
        };
    };

}

#endif //SUBSETOPTIMIZATION_GRAPHNEGATIVEGROUPFARNESS_H
