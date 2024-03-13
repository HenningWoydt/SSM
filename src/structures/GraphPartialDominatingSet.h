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

#ifndef SUBSETOPTIMIZATION_GRAPHPARTIALDOMINATINGSET_H
#define SUBSETOPTIMIZATION_GRAPHPARTIALDOMINATINGSET_H

#include "Graph.h"

namespace SMSM {

/**
 * Graph structure to optimize for negative group farness.
 */
    template<typename TypeSF>
    class GraphPartialDominatingSet final : public Graph<TypeSF> {
    public:
        using Graph<TypeSF>::Graph;

        // structures to speed up score function evaluation
        size_t depth = 0;
        std::vector<TypeSF> scores;
        std::vector<std::vector<uint8_t>> vertex_sets;
        std::vector<uint8_t> temp;

        inline TypeSF evaluate_empty_set() override {
            return 0;
        };

        inline TypeSF evaluate_1D(const std::vector<uint32_t> &s, const size_t s_size) override {
            TypeSF score = scores[depth];

            score += !vertex_sets[depth][s[s_size - 1]];

            for (uint32_t neighbour: Graph<TypeSF>::adj_list[s[s_size - 1]]) {
                score += !vertex_sets[depth][neighbour];
            }

            ASSERT(evaluate_general(s, s_size) == score);
            return score;
        };

        inline TypeSF evaluate_2D(const std::vector<uint32_t> &s, const size_t s_size) override {
            overwrite(temp, vertex_sets[depth]);
            TypeSF score = scores[depth];

            score += !temp[s[s_size - 2]];
            temp[s[s_size - 2]] = 1;
            score += !temp[s[s_size - 1]];
            temp[s[s_size - 1]] = 1;

            for (uint32_t neighbour: Graph<TypeSF>::adj_list[s[s_size - 2]]) {
                score += !temp[neighbour];
                temp[neighbour] = 1;
            }
            for (uint32_t neighbour: Graph<TypeSF>::adj_list[s[s_size - 1]]) {
                score += !temp[neighbour];
                temp[neighbour] = 1;
            }

            ASSERT(evaluate_general(s, s_size) == score);
            return score;
        };

        inline TypeSF evaluate_XD(const std::vector<uint32_t> &s, size_t s_size) override {
            size_t n_new_elements = s_size - depth;

            overwrite(temp, vertex_sets[depth]);
            TypeSF score = scores[depth];

            for (size_t i = 0; i < n_new_elements; ++i) {

                score += !temp[s[depth + i]];
                temp[s[depth + i]] = 1;

                for (uint32_t neighbour: Graph<TypeSF>::adj_list[s[depth + i]]) {
                    score += !temp[neighbour];
                    temp[neighbour] = 1;
                }
            }

            ASSERT(evaluate_general(s, s_size) == score);
            return score;
        };

        inline TypeSF evaluate_general(const std::vector<uint32_t> &s, size_t s_size) override {
            TypeSF score = 0;

            std::fill(temp.begin(), temp.end(), 0);

            for (size_t i = 0; i < s_size; ++i) {

                score += !temp[s[i]];
                temp[s[i]] = 1;

                for (uint32_t neighbour: Graph<TypeSF>::adj_list[s[i]]) {
                    score += !temp[neighbour];
                    temp[neighbour] = 1;
                }
            }

            return score;
        };

        inline void finalize() override {
            Graph<TypeSF>::sort_unique_neighbours();
            temp.resize(Graph<TypeSF>::n_nodes);

            Graph<TypeSF>::max_reachable_score = Graph<TypeSF>::n_nodes;
        };

        inline void initialize_helping_structures(size_t k) override {
            vertex_sets.clear();
            vertex_sets.resize(k + 1, std::vector<uint8_t>(Graph<TypeSF>::n_nodes, 0));

            scores.clear();
            scores.resize(k + 1, 0);

            depth = 0;
        };

        inline void visit_new_depth(const std::vector<uint32_t> &s, size_t s_size) override {
            depth += 1;

            scores[depth] = scores[depth - 1];
            overwrite(vertex_sets[depth], vertex_sets[depth - 1]);

            scores[depth] += !vertex_sets[depth][s[s_size - 1]];
            vertex_sets[depth][s[s_size - 1]] = 1;

            for (uint32_t neighbour: Graph<TypeSF>::adj_list[s[s_size - 1]]) {
                scores[depth] += !vertex_sets[depth][neighbour];
                vertex_sets[depth][neighbour] = 1;
            }
        };

        inline void return_from_last_depth() override {
            depth -= 1;
        };
    };

}

#endif //SUBSETOPTIMIZATION_GRAPHPARTIALDOMINATINGSET_H
