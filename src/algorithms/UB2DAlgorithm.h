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

#ifndef SUBSETOPTIMIZATION_UB2DALGORITHM_H
#define SUBSETOPTIMIZATION_UB2DALGORITHM_H

#include <cstdint>
#include <vector>
#include <cstddef>
#include <algorithm>

#include "../utility/Util.h"
#include "UB2DStorage.h"
#include "../../3rd_party_tools/blossom5-v2.05.src/PerfectMatching.h"
#include "../utility/VectorOfVectors.h"

namespace SMSM {

/**
 * Represents an edge in a graph with associated weight.
 *
 * This structure defines an edge in a graph, where 'n1' and 'n2' are the indices
 * of the two nodes connected by the edge, and 'weight' represents the weight or
 * cost associated with the edge.
 *
 * @tparam TypeSF The type of the edge weight (e.g. int, float, double).
 */
    template<typename TypeSF>
    struct Edge {
        size_t n1;
        size_t n2;
        TypeSF weight;
    };

/**
 * This class serves as a superclass for all methods that return an upper bound
 * with 2D score improvements.
 */
    template<typename TypeSF>
    class UB2DAlgorithm {
    private:
        // graph structure
        size_t n_nodes;
        size_t n_edges;
        std::vector<Edge<TypeSF>> edges;
        size_t edges_size;
        TypeSF max_weight = -std::numeric_limits<TypeSF>::max();

        // variable to hold the algorithm type
        size_t algorithm_type;

        // needed for greedy
        std::vector<Edge<TypeSF>> greedy_edges;
        bool greedy_is_sorted = false;

        // needed for matching
        PerfectMatching::Options options;
        PerfectMatching *pm = nullptr;
        size_t pm_last_n_endpoints = 0;
        PerfectMatching *pm_skip_first = nullptr;
        size_t pm_skip_first_last_n_endpoints = 0;

        // needed for brute force iterative
        std::vector<size_t> set;
        std::vector<size_t> best_set;
        std::vector<size_t> count;

        // needed for dynamic programming
        std::vector<uint8_t> solution_ready;
        std::vector<TypeSF> cache_si;
        std::vector<size_t> cache_bitset_t1;
        std::vector<size_t> cache_bitset_t2;
        std::vector<size_t> node_set;

    public:
        /**
         * Makes space for n vertices.
         *
         * @param n Number of vertices.
         * @param m Number of edges.
         */
        explicit UB2DAlgorithm(size_t n, size_t algo_type) {
            n_nodes = n;
            n_edges = (n * (n - 1)) / 2;

            algorithm_type = algo_type;

            edges_size = 0;

            initialize(n, algo_type);
        }

        /**
         * De-constructor.
         */
        ~UB2DAlgorithm() {
            delete pm;
            delete pm_skip_first;
        }

        /**
         * Clears all structures.
         */
        void clear() {
            n_nodes = 0;
            n_edges = 0;
            edges_size = 0;
            max_weight = -std::numeric_limits<TypeSF>::max();

            if (algorithm_type == 1) {
                greedy_is_sorted = false;
            }

            if (algorithm_type == 2) {
                delete pm;
                pm = nullptr;
                pm_last_n_endpoints = 0;

                delete pm_skip_first;
                pm_skip_first = nullptr;
                pm_skip_first_last_n_endpoints = 0;
            }

            if (algorithm_type == 4) {
                solution_ready.clear();
                cache_si.clear();
            }
        }

        /**
         * Reinitialize the graph with a new number of nodes and edges.
         *
         * @param n The new number of nodes in the graph.
         * @param m The new number of edges in the graph.
         */
        void reinitialize(size_t n, size_t algo_type) {
            clear();
            initialize(n, algo_type);
        }

        /**
         * Add an edge with a specified weight to the graph.
         *
         * @param u The index of the first node connected by the edge.
         * @param v The index of the second node connected by the edge.
         * @param weight The weight or cost associated with the edge.
         */
        void add_edge(size_t u, size_t v, TypeSF weight) {
            edges[edges_size].n1 = u;
            edges[edges_size].n2 = v;
            edges[edges_size].weight = weight;
            edges_size += 1;

            max_weight = std::max(max_weight, weight);
        }

        /**
         * Calculate an upper bound for a subset optimization problem based on the specified algorithm type.
         *
         * This function calculates an upper bound for a subset optimization problem using one of several
         * algorithms based on the provided `algorithm_type`:
         * - If `algorithm_type` is 1, it uses a greedy algorithm.
         * - If `algorithm_type` is 2, it uses a matching algorithm.
         * - If `algorithm_type` is 3, it uses a brute-force algorithm.
         * - If `algorithm_type` is 4, it uses a dynamic programming algorithm.
         *
         * @param v_s The vector of vectors containing subsets of endpoints.
         * @param w_s The vector of weights corresponding to the subsets.
         * @param n_endpoints The number of endpoints.
         * @param offset An offset value used in the calculation.
         * @return The calculated upper bound based on the selected algorithm.
         * @throws std::runtime_error if the specified `algorithm_type` is not recognized.
         */
        TypeSF get_upper_bound(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());
            ASSERT(n_endpoints <= n_nodes);

            ASSERT(n_endpoints == 0 || n_endpoints == 2 || n_endpoints == 4 || algorithm_type == 1 || algorithm_type == 2 || algorithm_type == 3 || algorithm_type == 4);

            TypeSF result;
            if (n_endpoints == 0) {
                result = 0;
            } else if (n_endpoints == 2) {
                result = two_endpoints(v_s, w_s, offset);
            } else if (n_endpoints == 4) {
                if (n_nodes <= 10) {
                    result = four_endpoints_10(v_s, w_s, offset);
                } else {
                    result = four_endpoints(v_s, w_s, offset);
                }
            } else if (algorithm_type == 1) {
                result = get_upper_bound_greedy(v_s, w_s, n_endpoints, offset);
            } else if (algorithm_type == 2) {
                result = get_upper_bound_matching(v_s, w_s, n_endpoints, offset);
            } else if (algorithm_type == 3) {
                if (n_nodes <= 10) {
                    result = get_upper_bound_bforce_lookup(v_s, w_s, n_endpoints, offset);
                } else {
                    result = get_upper_bound_bforce_iterative(v_s, w_s, n_endpoints, offset);
                }
            } else if (algorithm_type == 4) {
                result = get_upper_bound_dynamic(v_s, w_s, n_endpoints, offset);
            } else {
                std::cout << "Dont recognize algorithm type " << algorithm_type << " for UB2D!" << std::endl;
                exit(EXIT_FAILURE);
            }

            return result;
        }

        /**
         * Calculate an upper bound for a subset optimization problem with the first endpoint skipped, based on the specified algorithm type.
         *
         * This function calculates an upper bound for a subset optimization problem with the first endpoint skipped
         * using one of several algorithms based on the provided `algorithm_type`:
         * - If `algorithm_type` is 1, it uses a greedy algorithm.
         * - If `algorithm_type` is 2, it uses a matching algorithm.
         * - If `algorithm_type` is 3, it uses a brute-force algorithm.
         * - If `algorithm_type` is 4, it uses a dynamic programming algorithm.
         *
         * @param v_s The vector of vectors containing subsets of endpoints.
         * @param w_s The vector of weights corresponding to the subsets.
         * @param n_endpoints The number of endpoints.
         * @param offset An offset value used in the calculation.
         * @return The calculated upper bound based on the selected algorithm with the first endpoint skipped.
         * @throws std::runtime_error if the specified `algorithm_type` is not recognized.
         */
        TypeSF get_upper_bound_skip_first(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());
            ASSERT(n_endpoints <= n_nodes);

            ASSERT(n_endpoints == 0 || n_endpoints == 2 || n_endpoints == 4 || algorithm_type == 1 || algorithm_type == 2 || algorithm_type == 3 || algorithm_type == 4);

            TypeSF result;
            if (n_endpoints == 0) {
                result = 0;
            } else if (n_endpoints == 2) {
                result = two_endpoints_skip_first(v_s, w_s, offset);
            } else if (n_endpoints == 4) {
                if (n_nodes <= 10) {
                    result = four_endpoints_10_skip_first(v_s, w_s, offset);
                } else {
                    result = four_endpoints_skip_first(v_s, w_s, offset);
                }
            } else if (algorithm_type == 1) {
                result = get_upper_bound_skip_first_greedy(v_s, w_s, n_endpoints, offset);
            } else if (algorithm_type == 2) {
                result = get_upper_bound_skip_first_matching(v_s, w_s, n_endpoints, offset);
            } else if (algorithm_type == 3) {
                if (n_nodes <= 10) {
                    result = get_upper_bound_skip_first_bforce_lookup(v_s, w_s, n_endpoints, offset);
                } else {
                    result = get_upper_bound_skip_first_bforce_iterative(v_s, w_s, n_endpoints, offset);
                }
            } else if (algorithm_type == 4) {
                result = get_upper_bound_skip_first_dynamic(v_s, w_s, n_endpoints, offset);
            } else {
                std::cout << "Dont recognize algorithm type " << algorithm_type << " for UB2D!" << std::endl;
                exit(EXIT_FAILURE);
            }

            return result;
        }

    private:
        /**
         * Reserves memory for all needed data structures.
         *
         * @param n The number of nodes in the graph.
         * @param m The number of edges in the graph.
         */
        void initialize(size_t n, size_t algo_type) {
            n_nodes = n;
            n_edges = (n * (n - 1)) / 2;
            max_weight = -std::numeric_limits<TypeSF>::max();

            algorithm_type = algo_type;

            edges.resize(n_edges);
            edges_size = 0;

            if (algorithm_type == 2) {
                options.fractional_jumpstart = true; // true is better than false
                options.dual_greedy_update_option = 0; // 0, 1 better than 2
                options.dual_LP_threshold = 0.00;
                options.update_duals_before = false;
                options.update_duals_after = false;
                options.single_tree_threshold = 1.00;
                options.verbose = false;
            }

            if (algorithm_type == 4) {
                size_t c = ceil(n, 2);
                size_t cache_size = pow_2(n);
                solution_ready.resize(c);
                std::fill(solution_ready.begin(), solution_ready.end(), 0);

                cache_si.resize(cache_size);
                std::fill(cache_si.begin(), cache_si.end(), std::numeric_limits<TypeSF>::max());

                cache_bitset_t1.resize(cache_size);
                cache_bitset_t2.resize(cache_size);
            }
        }

        /**
         * Adds the edge with the maximum weight.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF two_endpoints(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            size_t best_i = 0;
            for (size_t i = 1; i < edges_size; ++i) {
                if (edges[i].weight > edges[best_i].weight) {
                    best_i = i;
                }
            }

            v_s.push_back(edges[best_i].n1 + offset, edges[best_i].n2 + offset);
            w_s.push_back(edges[best_i].weight);

            return edges[best_i].weight;
        }

#define UPDATE(i, j) if (edges[i].weight + edges[j].weight > edges[best_i].weight + edges[best_j].weight){ best_i = i; best_j = j; } ((void)0)
#define NO_OVERLAP(i, j) (edges[i].n1 != edges[j].n1 && edges[i].n1 != edges[j].n2 && edges[i].n2 != edges[j].n1 && edges[i].n2 != edges[j].n2)

        /**
         * Adds two edges that have the combined maximum weight.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF four_endpoints(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            size_t best_i = 0, best_j = n_nodes + n_nodes - 1;
            TypeSF best_sum = edges[best_i].weight + edges[best_j].weight;

            for (size_t i = 0; i < n_edges; ++i) {
                for (size_t j = i + 1; j < n_edges; ++j) {
                    if (edges[i].weight + edges[j].weight > best_sum && NO_OVERLAP(i, j)) {
                        best_i = i;
                        best_j = j;
                        best_sum = edges[i].weight + edges[j].weight;
                    }
                }
            }

            v_s.push_back(edges[best_i].n1 + offset, edges[best_i].n2 + offset);
            v_s.push_back(edges[best_j].n1 + offset, edges[best_j].n2 + offset);
            w_s.push_back(edges[best_i].weight);
            w_s.push_back(edges[best_j].weight);

            return best_sum;
        }

        TypeSF four_endpoints_10(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());
            ASSERT(n_nodes >= 4);
            ASSERT(n_nodes <= 10);

            size_t best_i = 0, best_j = n_edges - 1;

            if (n_nodes == 4) {
// 0
                UPDATE(0, 5);
// 1
                UPDATE(1, 4);
// 2
                UPDATE(2, 3);
// 3
// 4
// 5
            } else if (n_nodes == 5) {
// 0
                UPDATE(0, 7);
                UPDATE(0, 8);
                UPDATE(0, 9);
// 1
                UPDATE(1, 5);
                UPDATE(1, 6);
                UPDATE(1, 9);
// 2
                UPDATE(2, 4);
                UPDATE(2, 6);
                UPDATE(2, 8);
// 3
                UPDATE(3, 4);
                UPDATE(3, 5);
                UPDATE(3, 7);
// 4
                UPDATE(4, 9);
// 5
                UPDATE(5, 8);
// 6
                UPDATE(6, 7);
// 7
// 8
// 9
            } else if (n_nodes == 6) {
// 0
                UPDATE(0, 9);
                UPDATE(0, 10);
                UPDATE(0, 11);
                UPDATE(0, 12);
                UPDATE(0, 13);
                UPDATE(0, 14);
// 1
                UPDATE(1, 6);
                UPDATE(1, 7);
                UPDATE(1, 8);
                UPDATE(1, 12);
                UPDATE(1, 13);
                UPDATE(1, 14);
// 2
                UPDATE(2, 5);
                UPDATE(2, 7);
                UPDATE(2, 8);
                UPDATE(2, 10);
                UPDATE(2, 11);
                UPDATE(2, 14);
// 3
                UPDATE(3, 5);
                UPDATE(3, 6);
                UPDATE(3, 8);
                UPDATE(3, 9);
                UPDATE(3, 11);
                UPDATE(3, 13);
// 4
                UPDATE(4, 5);
                UPDATE(4, 6);
                UPDATE(4, 7);
                UPDATE(4, 9);
                UPDATE(4, 10);
                UPDATE(4, 12);
// 5
                UPDATE(5, 12);
                UPDATE(5, 13);
                UPDATE(5, 14);
// 6
                UPDATE(6, 10);
                UPDATE(6, 11);
                UPDATE(6, 14);
// 7
                UPDATE(7, 9);
                UPDATE(7, 11);
                UPDATE(7, 13);
// 8
                UPDATE(8, 9);
                UPDATE(8, 10);
                UPDATE(8, 12);
// 9
                UPDATE(9, 14);
// 10
                UPDATE(10, 13);
// 11
                UPDATE(11, 12);
// 12
// 13
// 14
            } else if (n_nodes == 7) {
// 0
                UPDATE(0, 11);
                UPDATE(0, 12);
                UPDATE(0, 13);
                UPDATE(0, 14);
                UPDATE(0, 15);
                UPDATE(0, 16);
                UPDATE(0, 17);
                UPDATE(0, 18);
                UPDATE(0, 19);
                UPDATE(0, 20);
// 1
                UPDATE(1, 7);
                UPDATE(1, 8);
                UPDATE(1, 9);
                UPDATE(1, 10);
                UPDATE(1, 15);
                UPDATE(1, 16);
                UPDATE(1, 17);
                UPDATE(1, 18);
                UPDATE(1, 19);
                UPDATE(1, 20);
// 2
                UPDATE(2, 6);
                UPDATE(2, 8);
                UPDATE(2, 9);
                UPDATE(2, 10);
                UPDATE(2, 12);
                UPDATE(2, 13);
                UPDATE(2, 14);
                UPDATE(2, 18);
                UPDATE(2, 19);
                UPDATE(2, 20);
// 3
                UPDATE(3, 6);
                UPDATE(3, 7);
                UPDATE(3, 9);
                UPDATE(3, 10);
                UPDATE(3, 11);
                UPDATE(3, 13);
                UPDATE(3, 14);
                UPDATE(3, 16);
                UPDATE(3, 17);
                UPDATE(3, 20);
// 4
                UPDATE(4, 6);
                UPDATE(4, 7);
                UPDATE(4, 8);
                UPDATE(4, 10);
                UPDATE(4, 11);
                UPDATE(4, 12);
                UPDATE(4, 14);
                UPDATE(4, 15);
                UPDATE(4, 17);
                UPDATE(4, 19);
// 5
                UPDATE(5, 6);
                UPDATE(5, 7);
                UPDATE(5, 8);
                UPDATE(5, 9);
                UPDATE(5, 11);
                UPDATE(5, 12);
                UPDATE(5, 13);
                UPDATE(5, 15);
                UPDATE(5, 16);
                UPDATE(5, 18);
// 6
                UPDATE(6, 15);
                UPDATE(6, 16);
                UPDATE(6, 17);
                UPDATE(6, 18);
                UPDATE(6, 19);
                UPDATE(6, 20);
// 7
                UPDATE(7, 12);
                UPDATE(7, 13);
                UPDATE(7, 14);
                UPDATE(7, 18);
                UPDATE(7, 19);
                UPDATE(7, 20);
// 8
                UPDATE(8, 11);
                UPDATE(8, 13);
                UPDATE(8, 14);
                UPDATE(8, 16);
                UPDATE(8, 17);
                UPDATE(8, 20);
// 9
                UPDATE(9, 11);
                UPDATE(9, 12);
                UPDATE(9, 14);
                UPDATE(9, 15);
                UPDATE(9, 17);
                UPDATE(9, 19);
// 10
                UPDATE(10, 11);
                UPDATE(10, 12);
                UPDATE(10, 13);
                UPDATE(10, 15);
                UPDATE(10, 16);
                UPDATE(10, 18);
// 11
                UPDATE(11, 18);
                UPDATE(11, 19);
                UPDATE(11, 20);
// 12
                UPDATE(12, 16);
                UPDATE(12, 17);
                UPDATE(12, 20);
// 13
                UPDATE(13, 15);
                UPDATE(13, 17);
                UPDATE(13, 19);
// 14
                UPDATE(14, 15);
                UPDATE(14, 16);
                UPDATE(14, 18);
// 15
                UPDATE(15, 20);
// 16
                UPDATE(16, 19);
// 17
                UPDATE(17, 18);
// 18
// 19
// 20
            } else if (n_nodes == 8) {
// 0
                UPDATE(0, 13);
                UPDATE(0, 14);
                UPDATE(0, 15);
                UPDATE(0, 16);
                UPDATE(0, 17);
                UPDATE(0, 18);
                UPDATE(0, 19);
                UPDATE(0, 20);
                UPDATE(0, 21);
                UPDATE(0, 22);
                UPDATE(0, 23);
                UPDATE(0, 24);
                UPDATE(0, 25);
                UPDATE(0, 26);
                UPDATE(0, 27);
// 1
                UPDATE(1, 8);
                UPDATE(1, 9);
                UPDATE(1, 10);
                UPDATE(1, 11);
                UPDATE(1, 12);
                UPDATE(1, 18);
                UPDATE(1, 19);
                UPDATE(1, 20);
                UPDATE(1, 21);
                UPDATE(1, 22);
                UPDATE(1, 23);
                UPDATE(1, 24);
                UPDATE(1, 25);
                UPDATE(1, 26);
                UPDATE(1, 27);
// 2
                UPDATE(2, 7);
                UPDATE(2, 9);
                UPDATE(2, 10);
                UPDATE(2, 11);
                UPDATE(2, 12);
                UPDATE(2, 14);
                UPDATE(2, 15);
                UPDATE(2, 16);
                UPDATE(2, 17);
                UPDATE(2, 22);
                UPDATE(2, 23);
                UPDATE(2, 24);
                UPDATE(2, 25);
                UPDATE(2, 26);
                UPDATE(2, 27);
// 3
                UPDATE(3, 7);
                UPDATE(3, 8);
                UPDATE(3, 10);
                UPDATE(3, 11);
                UPDATE(3, 12);
                UPDATE(3, 13);
                UPDATE(3, 15);
                UPDATE(3, 16);
                UPDATE(3, 17);
                UPDATE(3, 19);
                UPDATE(3, 20);
                UPDATE(3, 21);
                UPDATE(3, 25);
                UPDATE(3, 26);
                UPDATE(3, 27);
// 4
                UPDATE(4, 7);
                UPDATE(4, 8);
                UPDATE(4, 9);
                UPDATE(4, 11);
                UPDATE(4, 12);
                UPDATE(4, 13);
                UPDATE(4, 14);
                UPDATE(4, 16);
                UPDATE(4, 17);
                UPDATE(4, 18);
                UPDATE(4, 20);
                UPDATE(4, 21);
                UPDATE(4, 23);
                UPDATE(4, 24);
                UPDATE(4, 27);
// 5
                UPDATE(5, 7);
                UPDATE(5, 8);
                UPDATE(5, 9);
                UPDATE(5, 10);
                UPDATE(5, 12);
                UPDATE(5, 13);
                UPDATE(5, 14);
                UPDATE(5, 15);
                UPDATE(5, 17);
                UPDATE(5, 18);
                UPDATE(5, 19);
                UPDATE(5, 21);
                UPDATE(5, 22);
                UPDATE(5, 24);
                UPDATE(5, 26);
// 6
                UPDATE(6, 7);
                UPDATE(6, 8);
                UPDATE(6, 9);
                UPDATE(6, 10);
                UPDATE(6, 11);
                UPDATE(6, 13);
                UPDATE(6, 14);
                UPDATE(6, 15);
                UPDATE(6, 16);
                UPDATE(6, 18);
                UPDATE(6, 19);
                UPDATE(6, 20);
                UPDATE(6, 22);
                UPDATE(6, 23);
                UPDATE(6, 25);
// 7
                UPDATE(7, 18);
                UPDATE(7, 19);
                UPDATE(7, 20);
                UPDATE(7, 21);
                UPDATE(7, 22);
                UPDATE(7, 23);
                UPDATE(7, 24);
                UPDATE(7, 25);
                UPDATE(7, 26);
                UPDATE(7, 27);
// 8
                UPDATE(8, 14);
                UPDATE(8, 15);
                UPDATE(8, 16);
                UPDATE(8, 17);
                UPDATE(8, 22);
                UPDATE(8, 23);
                UPDATE(8, 24);
                UPDATE(8, 25);
                UPDATE(8, 26);
                UPDATE(8, 27);
// 9
                UPDATE(9, 13);
                UPDATE(9, 15);
                UPDATE(9, 16);
                UPDATE(9, 17);
                UPDATE(9, 19);
                UPDATE(9, 20);
                UPDATE(9, 21);
                UPDATE(9, 25);
                UPDATE(9, 26);
                UPDATE(9, 27);
// 10
                UPDATE(10, 13);
                UPDATE(10, 14);
                UPDATE(10, 16);
                UPDATE(10, 17);
                UPDATE(10, 18);
                UPDATE(10, 20);
                UPDATE(10, 21);
                UPDATE(10, 23);
                UPDATE(10, 24);
                UPDATE(10, 27);
// 11
                UPDATE(11, 13);
                UPDATE(11, 14);
                UPDATE(11, 15);
                UPDATE(11, 17);
                UPDATE(11, 18);
                UPDATE(11, 19);
                UPDATE(11, 21);
                UPDATE(11, 22);
                UPDATE(11, 24);
                UPDATE(11, 26);
// 12
                UPDATE(12, 13);
                UPDATE(12, 14);
                UPDATE(12, 15);
                UPDATE(12, 16);
                UPDATE(12, 18);
                UPDATE(12, 19);
                UPDATE(12, 20);
                UPDATE(12, 22);
                UPDATE(12, 23);
                UPDATE(12, 25);
// 13
                UPDATE(13, 22);
                UPDATE(13, 23);
                UPDATE(13, 24);
                UPDATE(13, 25);
                UPDATE(13, 26);
                UPDATE(13, 27);
// 14
                UPDATE(14, 19);
                UPDATE(14, 20);
                UPDATE(14, 21);
                UPDATE(14, 25);
                UPDATE(14, 26);
                UPDATE(14, 27);
// 15
                UPDATE(15, 18);
                UPDATE(15, 20);
                UPDATE(15, 21);
                UPDATE(15, 23);
                UPDATE(15, 24);
                UPDATE(15, 27);
// 16
                UPDATE(16, 18);
                UPDATE(16, 19);
                UPDATE(16, 21);
                UPDATE(16, 22);
                UPDATE(16, 24);
                UPDATE(16, 26);
// 17
                UPDATE(17, 18);
                UPDATE(17, 19);
                UPDATE(17, 20);
                UPDATE(17, 22);
                UPDATE(17, 23);
                UPDATE(17, 25);
// 18
                UPDATE(18, 25);
                UPDATE(18, 26);
                UPDATE(18, 27);
// 19
                UPDATE(19, 23);
                UPDATE(19, 24);
                UPDATE(19, 27);
// 20
                UPDATE(20, 22);
                UPDATE(20, 24);
                UPDATE(20, 26);
// 21
                UPDATE(21, 22);
                UPDATE(21, 23);
                UPDATE(21, 25);
// 22
                UPDATE(22, 27);
// 23
                UPDATE(23, 26);
// 24
                UPDATE(24, 25);
// 25
// 26
// 27
            } else if (n_nodes == 9) {
// 0
                UPDATE(0, 15);
                UPDATE(0, 16);
                UPDATE(0, 17);
                UPDATE(0, 18);
                UPDATE(0, 19);
                UPDATE(0, 20);
                UPDATE(0, 21);
                UPDATE(0, 22);
                UPDATE(0, 23);
                UPDATE(0, 24);
                UPDATE(0, 25);
                UPDATE(0, 26);
                UPDATE(0, 27);
                UPDATE(0, 28);
                UPDATE(0, 29);
                UPDATE(0, 30);
                UPDATE(0, 31);
                UPDATE(0, 32);
                UPDATE(0, 33);
                UPDATE(0, 34);
                UPDATE(0, 35);
// 1
                UPDATE(1, 9);
                UPDATE(1, 10);
                UPDATE(1, 11);
                UPDATE(1, 12);
                UPDATE(1, 13);
                UPDATE(1, 14);
                UPDATE(1, 21);
                UPDATE(1, 22);
                UPDATE(1, 23);
                UPDATE(1, 24);
                UPDATE(1, 25);
                UPDATE(1, 26);
                UPDATE(1, 27);
                UPDATE(1, 28);
                UPDATE(1, 29);
                UPDATE(1, 30);
                UPDATE(1, 31);
                UPDATE(1, 32);
                UPDATE(1, 33);
                UPDATE(1, 34);
                UPDATE(1, 35);
// 2
                UPDATE(2, 8);
                UPDATE(2, 10);
                UPDATE(2, 11);
                UPDATE(2, 12);
                UPDATE(2, 13);
                UPDATE(2, 14);
                UPDATE(2, 16);
                UPDATE(2, 17);
                UPDATE(2, 18);
                UPDATE(2, 19);
                UPDATE(2, 20);
                UPDATE(2, 26);
                UPDATE(2, 27);
                UPDATE(2, 28);
                UPDATE(2, 29);
                UPDATE(2, 30);
                UPDATE(2, 31);
                UPDATE(2, 32);
                UPDATE(2, 33);
                UPDATE(2, 34);
                UPDATE(2, 35);
// 3
                UPDATE(3, 8);
                UPDATE(3, 9);
                UPDATE(3, 11);
                UPDATE(3, 12);
                UPDATE(3, 13);
                UPDATE(3, 14);
                UPDATE(3, 15);
                UPDATE(3, 17);
                UPDATE(3, 18);
                UPDATE(3, 19);
                UPDATE(3, 20);
                UPDATE(3, 22);
                UPDATE(3, 23);
                UPDATE(3, 24);
                UPDATE(3, 25);
                UPDATE(3, 30);
                UPDATE(3, 31);
                UPDATE(3, 32);
                UPDATE(3, 33);
                UPDATE(3, 34);
                UPDATE(3, 35);
// 4
                UPDATE(4, 8);
                UPDATE(4, 9);
                UPDATE(4, 10);
                UPDATE(4, 12);
                UPDATE(4, 13);
                UPDATE(4, 14);
                UPDATE(4, 15);
                UPDATE(4, 16);
                UPDATE(4, 18);
                UPDATE(4, 19);
                UPDATE(4, 20);
                UPDATE(4, 21);
                UPDATE(4, 23);
                UPDATE(4, 24);
                UPDATE(4, 25);
                UPDATE(4, 27);
                UPDATE(4, 28);
                UPDATE(4, 29);
                UPDATE(4, 33);
                UPDATE(4, 34);
                UPDATE(4, 35);
// 5
                UPDATE(5, 8);
                UPDATE(5, 9);
                UPDATE(5, 10);
                UPDATE(5, 11);
                UPDATE(5, 13);
                UPDATE(5, 14);
                UPDATE(5, 15);
                UPDATE(5, 16);
                UPDATE(5, 17);
                UPDATE(5, 19);
                UPDATE(5, 20);
                UPDATE(5, 21);
                UPDATE(5, 22);
                UPDATE(5, 24);
                UPDATE(5, 25);
                UPDATE(5, 26);
                UPDATE(5, 28);
                UPDATE(5, 29);
                UPDATE(5, 31);
                UPDATE(5, 32);
                UPDATE(5, 35);
// 6
                UPDATE(6, 8);
                UPDATE(6, 9);
                UPDATE(6, 10);
                UPDATE(6, 11);
                UPDATE(6, 12);
                UPDATE(6, 14);
                UPDATE(6, 15);
                UPDATE(6, 16);
                UPDATE(6, 17);
                UPDATE(6, 18);
                UPDATE(6, 20);
                UPDATE(6, 21);
                UPDATE(6, 22);
                UPDATE(6, 23);
                UPDATE(6, 25);
                UPDATE(6, 26);
                UPDATE(6, 27);
                UPDATE(6, 29);
                UPDATE(6, 30);
                UPDATE(6, 32);
                UPDATE(6, 34);
// 7
                UPDATE(7, 8);
                UPDATE(7, 9);
                UPDATE(7, 10);
                UPDATE(7, 11);
                UPDATE(7, 12);
                UPDATE(7, 13);
                UPDATE(7, 15);
                UPDATE(7, 16);
                UPDATE(7, 17);
                UPDATE(7, 18);
                UPDATE(7, 19);
                UPDATE(7, 21);
                UPDATE(7, 22);
                UPDATE(7, 23);
                UPDATE(7, 24);
                UPDATE(7, 26);
                UPDATE(7, 27);
                UPDATE(7, 28);
                UPDATE(7, 30);
                UPDATE(7, 31);
                UPDATE(7, 33);
// 8
                UPDATE(8, 21);
                UPDATE(8, 22);
                UPDATE(8, 23);
                UPDATE(8, 24);
                UPDATE(8, 25);
                UPDATE(8, 26);
                UPDATE(8, 27);
                UPDATE(8, 28);
                UPDATE(8, 29);
                UPDATE(8, 30);
                UPDATE(8, 31);
                UPDATE(8, 32);
                UPDATE(8, 33);
                UPDATE(8, 34);
                UPDATE(8, 35);
// 9
                UPDATE(9, 16);
                UPDATE(9, 17);
                UPDATE(9, 18);
                UPDATE(9, 19);
                UPDATE(9, 20);
                UPDATE(9, 26);
                UPDATE(9, 27);
                UPDATE(9, 28);
                UPDATE(9, 29);
                UPDATE(9, 30);
                UPDATE(9, 31);
                UPDATE(9, 32);
                UPDATE(9, 33);
                UPDATE(9, 34);
                UPDATE(9, 35);
// 10
                UPDATE(10, 15);
                UPDATE(10, 17);
                UPDATE(10, 18);
                UPDATE(10, 19);
                UPDATE(10, 20);
                UPDATE(10, 22);
                UPDATE(10, 23);
                UPDATE(10, 24);
                UPDATE(10, 25);
                UPDATE(10, 30);
                UPDATE(10, 31);
                UPDATE(10, 32);
                UPDATE(10, 33);
                UPDATE(10, 34);
                UPDATE(10, 35);
// 11
                UPDATE(11, 15);
                UPDATE(11, 16);
                UPDATE(11, 18);
                UPDATE(11, 19);
                UPDATE(11, 20);
                UPDATE(11, 21);
                UPDATE(11, 23);
                UPDATE(11, 24);
                UPDATE(11, 25);
                UPDATE(11, 27);
                UPDATE(11, 28);
                UPDATE(11, 29);
                UPDATE(11, 33);
                UPDATE(11, 34);
                UPDATE(11, 35);
// 12
                UPDATE(12, 15);
                UPDATE(12, 16);
                UPDATE(12, 17);
                UPDATE(12, 19);
                UPDATE(12, 20);
                UPDATE(12, 21);
                UPDATE(12, 22);
                UPDATE(12, 24);
                UPDATE(12, 25);
                UPDATE(12, 26);
                UPDATE(12, 28);
                UPDATE(12, 29);
                UPDATE(12, 31);
                UPDATE(12, 32);
                UPDATE(12, 35);
// 13
                UPDATE(13, 15);
                UPDATE(13, 16);
                UPDATE(13, 17);
                UPDATE(13, 18);
                UPDATE(13, 20);
                UPDATE(13, 21);
                UPDATE(13, 22);
                UPDATE(13, 23);
                UPDATE(13, 25);
                UPDATE(13, 26);
                UPDATE(13, 27);
                UPDATE(13, 29);
                UPDATE(13, 30);
                UPDATE(13, 32);
                UPDATE(13, 34);
// 14
                UPDATE(14, 15);
                UPDATE(14, 16);
                UPDATE(14, 17);
                UPDATE(14, 18);
                UPDATE(14, 19);
                UPDATE(14, 21);
                UPDATE(14, 22);
                UPDATE(14, 23);
                UPDATE(14, 24);
                UPDATE(14, 26);
                UPDATE(14, 27);
                UPDATE(14, 28);
                UPDATE(14, 30);
                UPDATE(14, 31);
                UPDATE(14, 33);
// 15
                UPDATE(15, 26);
                UPDATE(15, 27);
                UPDATE(15, 28);
                UPDATE(15, 29);
                UPDATE(15, 30);
                UPDATE(15, 31);
                UPDATE(15, 32);
                UPDATE(15, 33);
                UPDATE(15, 34);
                UPDATE(15, 35);
// 16
                UPDATE(16, 22);
                UPDATE(16, 23);
                UPDATE(16, 24);
                UPDATE(16, 25);
                UPDATE(16, 30);
                UPDATE(16, 31);
                UPDATE(16, 32);
                UPDATE(16, 33);
                UPDATE(16, 34);
                UPDATE(16, 35);
// 17
                UPDATE(17, 21);
                UPDATE(17, 23);
                UPDATE(17, 24);
                UPDATE(17, 25);
                UPDATE(17, 27);
                UPDATE(17, 28);
                UPDATE(17, 29);
                UPDATE(17, 33);
                UPDATE(17, 34);
                UPDATE(17, 35);
// 18
                UPDATE(18, 21);
                UPDATE(18, 22);
                UPDATE(18, 24);
                UPDATE(18, 25);
                UPDATE(18, 26);
                UPDATE(18, 28);
                UPDATE(18, 29);
                UPDATE(18, 31);
                UPDATE(18, 32);
                UPDATE(18, 35);
// 19
                UPDATE(19, 21);
                UPDATE(19, 22);
                UPDATE(19, 23);
                UPDATE(19, 25);
                UPDATE(19, 26);
                UPDATE(19, 27);
                UPDATE(19, 29);
                UPDATE(19, 30);
                UPDATE(19, 32);
                UPDATE(19, 34);
// 20
                UPDATE(20, 21);
                UPDATE(20, 22);
                UPDATE(20, 23);
                UPDATE(20, 24);
                UPDATE(20, 26);
                UPDATE(20, 27);
                UPDATE(20, 28);
                UPDATE(20, 30);
                UPDATE(20, 31);
                UPDATE(20, 33);
// 21
                UPDATE(21, 30);
                UPDATE(21, 31);
                UPDATE(21, 32);
                UPDATE(21, 33);
                UPDATE(21, 34);
                UPDATE(21, 35);
// 22
                UPDATE(22, 27);
                UPDATE(22, 28);
                UPDATE(22, 29);
                UPDATE(22, 33);
                UPDATE(22, 34);
                UPDATE(22, 35);
// 23
                UPDATE(23, 26);
                UPDATE(23, 28);
                UPDATE(23, 29);
                UPDATE(23, 31);
                UPDATE(23, 32);
                UPDATE(23, 35);
// 24
                UPDATE(24, 26);
                UPDATE(24, 27);
                UPDATE(24, 29);
                UPDATE(24, 30);
                UPDATE(24, 32);
                UPDATE(24, 34);
// 25
                UPDATE(25, 26);
                UPDATE(25, 27);
                UPDATE(25, 28);
                UPDATE(25, 30);
                UPDATE(25, 31);
                UPDATE(25, 33);
// 26
                UPDATE(26, 33);
                UPDATE(26, 34);
                UPDATE(26, 35);
// 27
                UPDATE(27, 31);
                UPDATE(27, 32);
                UPDATE(27, 35);
// 28
                UPDATE(28, 30);
                UPDATE(28, 32);
                UPDATE(28, 34);
// 29
                UPDATE(29, 30);
                UPDATE(29, 31);
                UPDATE(29, 33);
// 30
                UPDATE(30, 35);
// 31
                UPDATE(31, 34);
// 32
                UPDATE(32, 33);
// 33
// 34
// 35
            } else if (n_nodes == 10) {
// 0
                UPDATE(0, 17);
                UPDATE(0, 18);
                UPDATE(0, 19);
                UPDATE(0, 20);
                UPDATE(0, 21);
                UPDATE(0, 22);
                UPDATE(0, 23);
                UPDATE(0, 24);
                UPDATE(0, 25);
                UPDATE(0, 26);
                UPDATE(0, 27);
                UPDATE(0, 28);
                UPDATE(0, 29);
                UPDATE(0, 30);
                UPDATE(0, 31);
                UPDATE(0, 32);
                UPDATE(0, 33);
                UPDATE(0, 34);
                UPDATE(0, 35);
                UPDATE(0, 36);
                UPDATE(0, 37);
                UPDATE(0, 38);
                UPDATE(0, 39);
                UPDATE(0, 40);
                UPDATE(0, 41);
                UPDATE(0, 42);
                UPDATE(0, 43);
                UPDATE(0, 44);
// 1
                UPDATE(1, 10);
                UPDATE(1, 11);
                UPDATE(1, 12);
                UPDATE(1, 13);
                UPDATE(1, 14);
                UPDATE(1, 15);
                UPDATE(1, 16);
                UPDATE(1, 24);
                UPDATE(1, 25);
                UPDATE(1, 26);
                UPDATE(1, 27);
                UPDATE(1, 28);
                UPDATE(1, 29);
                UPDATE(1, 30);
                UPDATE(1, 31);
                UPDATE(1, 32);
                UPDATE(1, 33);
                UPDATE(1, 34);
                UPDATE(1, 35);
                UPDATE(1, 36);
                UPDATE(1, 37);
                UPDATE(1, 38);
                UPDATE(1, 39);
                UPDATE(1, 40);
                UPDATE(1, 41);
                UPDATE(1, 42);
                UPDATE(1, 43);
                UPDATE(1, 44);
// 2
                UPDATE(2, 9);
                UPDATE(2, 11);
                UPDATE(2, 12);
                UPDATE(2, 13);
                UPDATE(2, 14);
                UPDATE(2, 15);
                UPDATE(2, 16);
                UPDATE(2, 18);
                UPDATE(2, 19);
                UPDATE(2, 20);
                UPDATE(2, 21);
                UPDATE(2, 22);
                UPDATE(2, 23);
                UPDATE(2, 30);
                UPDATE(2, 31);
                UPDATE(2, 32);
                UPDATE(2, 33);
                UPDATE(2, 34);
                UPDATE(2, 35);
                UPDATE(2, 36);
                UPDATE(2, 37);
                UPDATE(2, 38);
                UPDATE(2, 39);
                UPDATE(2, 40);
                UPDATE(2, 41);
                UPDATE(2, 42);
                UPDATE(2, 43);
                UPDATE(2, 44);
// 3
                UPDATE(3, 9);
                UPDATE(3, 10);
                UPDATE(3, 12);
                UPDATE(3, 13);
                UPDATE(3, 14);
                UPDATE(3, 15);
                UPDATE(3, 16);
                UPDATE(3, 17);
                UPDATE(3, 19);
                UPDATE(3, 20);
                UPDATE(3, 21);
                UPDATE(3, 22);
                UPDATE(3, 23);
                UPDATE(3, 25);
                UPDATE(3, 26);
                UPDATE(3, 27);
                UPDATE(3, 28);
                UPDATE(3, 29);
                UPDATE(3, 35);
                UPDATE(3, 36);
                UPDATE(3, 37);
                UPDATE(3, 38);
                UPDATE(3, 39);
                UPDATE(3, 40);
                UPDATE(3, 41);
                UPDATE(3, 42);
                UPDATE(3, 43);
                UPDATE(3, 44);
// 4
                UPDATE(4, 9);
                UPDATE(4, 10);
                UPDATE(4, 11);
                UPDATE(4, 13);
                UPDATE(4, 14);
                UPDATE(4, 15);
                UPDATE(4, 16);
                UPDATE(4, 17);
                UPDATE(4, 18);
                UPDATE(4, 20);
                UPDATE(4, 21);
                UPDATE(4, 22);
                UPDATE(4, 23);
                UPDATE(4, 24);
                UPDATE(4, 26);
                UPDATE(4, 27);
                UPDATE(4, 28);
                UPDATE(4, 29);
                UPDATE(4, 31);
                UPDATE(4, 32);
                UPDATE(4, 33);
                UPDATE(4, 34);
                UPDATE(4, 39);
                UPDATE(4, 40);
                UPDATE(4, 41);
                UPDATE(4, 42);
                UPDATE(4, 43);
                UPDATE(4, 44);
// 5
                UPDATE(5, 9);
                UPDATE(5, 10);
                UPDATE(5, 11);
                UPDATE(5, 12);
                UPDATE(5, 14);
                UPDATE(5, 15);
                UPDATE(5, 16);
                UPDATE(5, 17);
                UPDATE(5, 18);
                UPDATE(5, 19);
                UPDATE(5, 21);
                UPDATE(5, 22);
                UPDATE(5, 23);
                UPDATE(5, 24);
                UPDATE(5, 25);
                UPDATE(5, 27);
                UPDATE(5, 28);
                UPDATE(5, 29);
                UPDATE(5, 30);
                UPDATE(5, 32);
                UPDATE(5, 33);
                UPDATE(5, 34);
                UPDATE(5, 36);
                UPDATE(5, 37);
                UPDATE(5, 38);
                UPDATE(5, 42);
                UPDATE(5, 43);
                UPDATE(5, 44);
// 6
                UPDATE(6, 9);
                UPDATE(6, 10);
                UPDATE(6, 11);
                UPDATE(6, 12);
                UPDATE(6, 13);
                UPDATE(6, 15);
                UPDATE(6, 16);
                UPDATE(6, 17);
                UPDATE(6, 18);
                UPDATE(6, 19);
                UPDATE(6, 20);
                UPDATE(6, 22);
                UPDATE(6, 23);
                UPDATE(6, 24);
                UPDATE(6, 25);
                UPDATE(6, 26);
                UPDATE(6, 28);
                UPDATE(6, 29);
                UPDATE(6, 30);
                UPDATE(6, 31);
                UPDATE(6, 33);
                UPDATE(6, 34);
                UPDATE(6, 35);
                UPDATE(6, 37);
                UPDATE(6, 38);
                UPDATE(6, 40);
                UPDATE(6, 41);
                UPDATE(6, 44);
// 7
                UPDATE(7, 9);
                UPDATE(7, 10);
                UPDATE(7, 11);
                UPDATE(7, 12);
                UPDATE(7, 13);
                UPDATE(7, 14);
                UPDATE(7, 16);
                UPDATE(7, 17);
                UPDATE(7, 18);
                UPDATE(7, 19);
                UPDATE(7, 20);
                UPDATE(7, 21);
                UPDATE(7, 23);
                UPDATE(7, 24);
                UPDATE(7, 25);
                UPDATE(7, 26);
                UPDATE(7, 27);
                UPDATE(7, 29);
                UPDATE(7, 30);
                UPDATE(7, 31);
                UPDATE(7, 32);
                UPDATE(7, 34);
                UPDATE(7, 35);
                UPDATE(7, 36);
                UPDATE(7, 38);
                UPDATE(7, 39);
                UPDATE(7, 41);
                UPDATE(7, 43);
// 8
                UPDATE(8, 9);
                UPDATE(8, 10);
                UPDATE(8, 11);
                UPDATE(8, 12);
                UPDATE(8, 13);
                UPDATE(8, 14);
                UPDATE(8, 15);
                UPDATE(8, 17);
                UPDATE(8, 18);
                UPDATE(8, 19);
                UPDATE(8, 20);
                UPDATE(8, 21);
                UPDATE(8, 22);
                UPDATE(8, 24);
                UPDATE(8, 25);
                UPDATE(8, 26);
                UPDATE(8, 27);
                UPDATE(8, 28);
                UPDATE(8, 30);
                UPDATE(8, 31);
                UPDATE(8, 32);
                UPDATE(8, 33);
                UPDATE(8, 35);
                UPDATE(8, 36);
                UPDATE(8, 37);
                UPDATE(8, 39);
                UPDATE(8, 40);
                UPDATE(8, 42);
// 9
                UPDATE(9, 24);
                UPDATE(9, 25);
                UPDATE(9, 26);
                UPDATE(9, 27);
                UPDATE(9, 28);
                UPDATE(9, 29);
                UPDATE(9, 30);
                UPDATE(9, 31);
                UPDATE(9, 32);
                UPDATE(9, 33);
                UPDATE(9, 34);
                UPDATE(9, 35);
                UPDATE(9, 36);
                UPDATE(9, 37);
                UPDATE(9, 38);
                UPDATE(9, 39);
                UPDATE(9, 40);
                UPDATE(9, 41);
                UPDATE(9, 42);
                UPDATE(9, 43);
                UPDATE(9, 44);
// 10
                UPDATE(10, 18);
                UPDATE(10, 19);
                UPDATE(10, 20);
                UPDATE(10, 21);
                UPDATE(10, 22);
                UPDATE(10, 23);
                UPDATE(10, 30);
                UPDATE(10, 31);
                UPDATE(10, 32);
                UPDATE(10, 33);
                UPDATE(10, 34);
                UPDATE(10, 35);
                UPDATE(10, 36);
                UPDATE(10, 37);
                UPDATE(10, 38);
                UPDATE(10, 39);
                UPDATE(10, 40);
                UPDATE(10, 41);
                UPDATE(10, 42);
                UPDATE(10, 43);
                UPDATE(10, 44);
// 11
                UPDATE(11, 17);
                UPDATE(11, 19);
                UPDATE(11, 20);
                UPDATE(11, 21);
                UPDATE(11, 22);
                UPDATE(11, 23);
                UPDATE(11, 25);
                UPDATE(11, 26);
                UPDATE(11, 27);
                UPDATE(11, 28);
                UPDATE(11, 29);
                UPDATE(11, 35);
                UPDATE(11, 36);
                UPDATE(11, 37);
                UPDATE(11, 38);
                UPDATE(11, 39);
                UPDATE(11, 40);
                UPDATE(11, 41);
                UPDATE(11, 42);
                UPDATE(11, 43);
                UPDATE(11, 44);
// 12
                UPDATE(12, 17);
                UPDATE(12, 18);
                UPDATE(12, 20);
                UPDATE(12, 21);
                UPDATE(12, 22);
                UPDATE(12, 23);
                UPDATE(12, 24);
                UPDATE(12, 26);
                UPDATE(12, 27);
                UPDATE(12, 28);
                UPDATE(12, 29);
                UPDATE(12, 31);
                UPDATE(12, 32);
                UPDATE(12, 33);
                UPDATE(12, 34);
                UPDATE(12, 39);
                UPDATE(12, 40);
                UPDATE(12, 41);
                UPDATE(12, 42);
                UPDATE(12, 43);
                UPDATE(12, 44);
// 13
                UPDATE(13, 17);
                UPDATE(13, 18);
                UPDATE(13, 19);
                UPDATE(13, 21);
                UPDATE(13, 22);
                UPDATE(13, 23);
                UPDATE(13, 24);
                UPDATE(13, 25);
                UPDATE(13, 27);
                UPDATE(13, 28);
                UPDATE(13, 29);
                UPDATE(13, 30);
                UPDATE(13, 32);
                UPDATE(13, 33);
                UPDATE(13, 34);
                UPDATE(13, 36);
                UPDATE(13, 37);
                UPDATE(13, 38);
                UPDATE(13, 42);
                UPDATE(13, 43);
                UPDATE(13, 44);
// 14
                UPDATE(14, 17);
                UPDATE(14, 18);
                UPDATE(14, 19);
                UPDATE(14, 20);
                UPDATE(14, 22);
                UPDATE(14, 23);
                UPDATE(14, 24);
                UPDATE(14, 25);
                UPDATE(14, 26);
                UPDATE(14, 28);
                UPDATE(14, 29);
                UPDATE(14, 30);
                UPDATE(14, 31);
                UPDATE(14, 33);
                UPDATE(14, 34);
                UPDATE(14, 35);
                UPDATE(14, 37);
                UPDATE(14, 38);
                UPDATE(14, 40);
                UPDATE(14, 41);
                UPDATE(14, 44);
// 15
                UPDATE(15, 17);
                UPDATE(15, 18);
                UPDATE(15, 19);
                UPDATE(15, 20);
                UPDATE(15, 21);
                UPDATE(15, 23);
                UPDATE(15, 24);
                UPDATE(15, 25);
                UPDATE(15, 26);
                UPDATE(15, 27);
                UPDATE(15, 29);
                UPDATE(15, 30);
                UPDATE(15, 31);
                UPDATE(15, 32);
                UPDATE(15, 34);
                UPDATE(15, 35);
                UPDATE(15, 36);
                UPDATE(15, 38);
                UPDATE(15, 39);
                UPDATE(15, 41);
                UPDATE(15, 43);
// 16
                UPDATE(16, 17);
                UPDATE(16, 18);
                UPDATE(16, 19);
                UPDATE(16, 20);
                UPDATE(16, 21);
                UPDATE(16, 22);
                UPDATE(16, 24);
                UPDATE(16, 25);
                UPDATE(16, 26);
                UPDATE(16, 27);
                UPDATE(16, 28);
                UPDATE(16, 30);
                UPDATE(16, 31);
                UPDATE(16, 32);
                UPDATE(16, 33);
                UPDATE(16, 35);
                UPDATE(16, 36);
                UPDATE(16, 37);
                UPDATE(16, 39);
                UPDATE(16, 40);
                UPDATE(16, 42);
// 17
                UPDATE(17, 30);
                UPDATE(17, 31);
                UPDATE(17, 32);
                UPDATE(17, 33);
                UPDATE(17, 34);
                UPDATE(17, 35);
                UPDATE(17, 36);
                UPDATE(17, 37);
                UPDATE(17, 38);
                UPDATE(17, 39);
                UPDATE(17, 40);
                UPDATE(17, 41);
                UPDATE(17, 42);
                UPDATE(17, 43);
                UPDATE(17, 44);
// 18
                UPDATE(18, 25);
                UPDATE(18, 26);
                UPDATE(18, 27);
                UPDATE(18, 28);
                UPDATE(18, 29);
                UPDATE(18, 35);
                UPDATE(18, 36);
                UPDATE(18, 37);
                UPDATE(18, 38);
                UPDATE(18, 39);
                UPDATE(18, 40);
                UPDATE(18, 41);
                UPDATE(18, 42);
                UPDATE(18, 43);
                UPDATE(18, 44);
// 19
                UPDATE(19, 24);
                UPDATE(19, 26);
                UPDATE(19, 27);
                UPDATE(19, 28);
                UPDATE(19, 29);
                UPDATE(19, 31);
                UPDATE(19, 32);
                UPDATE(19, 33);
                UPDATE(19, 34);
                UPDATE(19, 39);
                UPDATE(19, 40);
                UPDATE(19, 41);
                UPDATE(19, 42);
                UPDATE(19, 43);
                UPDATE(19, 44);
// 20
                UPDATE(20, 24);
                UPDATE(20, 25);
                UPDATE(20, 27);
                UPDATE(20, 28);
                UPDATE(20, 29);
                UPDATE(20, 30);
                UPDATE(20, 32);
                UPDATE(20, 33);
                UPDATE(20, 34);
                UPDATE(20, 36);
                UPDATE(20, 37);
                UPDATE(20, 38);
                UPDATE(20, 42);
                UPDATE(20, 43);
                UPDATE(20, 44);
// 21
                UPDATE(21, 24);
                UPDATE(21, 25);
                UPDATE(21, 26);
                UPDATE(21, 28);
                UPDATE(21, 29);
                UPDATE(21, 30);
                UPDATE(21, 31);
                UPDATE(21, 33);
                UPDATE(21, 34);
                UPDATE(21, 35);
                UPDATE(21, 37);
                UPDATE(21, 38);
                UPDATE(21, 40);
                UPDATE(21, 41);
                UPDATE(21, 44);
// 22
                UPDATE(22, 24);
                UPDATE(22, 25);
                UPDATE(22, 26);
                UPDATE(22, 27);
                UPDATE(22, 29);
                UPDATE(22, 30);
                UPDATE(22, 31);
                UPDATE(22, 32);
                UPDATE(22, 34);
                UPDATE(22, 35);
                UPDATE(22, 36);
                UPDATE(22, 38);
                UPDATE(22, 39);
                UPDATE(22, 41);
                UPDATE(22, 43);
// 23
                UPDATE(23, 24);
                UPDATE(23, 25);
                UPDATE(23, 26);
                UPDATE(23, 27);
                UPDATE(23, 28);
                UPDATE(23, 30);
                UPDATE(23, 31);
                UPDATE(23, 32);
                UPDATE(23, 33);
                UPDATE(23, 35);
                UPDATE(23, 36);
                UPDATE(23, 37);
                UPDATE(23, 39);
                UPDATE(23, 40);
                UPDATE(23, 42);
// 24
                UPDATE(24, 35);
                UPDATE(24, 36);
                UPDATE(24, 37);
                UPDATE(24, 38);
                UPDATE(24, 39);
                UPDATE(24, 40);
                UPDATE(24, 41);
                UPDATE(24, 42);
                UPDATE(24, 43);
                UPDATE(24, 44);
// 25
                UPDATE(25, 31);
                UPDATE(25, 32);
                UPDATE(25, 33);
                UPDATE(25, 34);
                UPDATE(25, 39);
                UPDATE(25, 40);
                UPDATE(25, 41);
                UPDATE(25, 42);
                UPDATE(25, 43);
                UPDATE(25, 44);
// 26
                UPDATE(26, 30);
                UPDATE(26, 32);
                UPDATE(26, 33);
                UPDATE(26, 34);
                UPDATE(26, 36);
                UPDATE(26, 37);
                UPDATE(26, 38);
                UPDATE(26, 42);
                UPDATE(26, 43);
                UPDATE(26, 44);
// 27
                UPDATE(27, 30);
                UPDATE(27, 31);
                UPDATE(27, 33);
                UPDATE(27, 34);
                UPDATE(27, 35);
                UPDATE(27, 37);
                UPDATE(27, 38);
                UPDATE(27, 40);
                UPDATE(27, 41);
                UPDATE(27, 44);
// 28
                UPDATE(28, 30);
                UPDATE(28, 31);
                UPDATE(28, 32);
                UPDATE(28, 34);
                UPDATE(28, 35);
                UPDATE(28, 36);
                UPDATE(28, 38);
                UPDATE(28, 39);
                UPDATE(28, 41);
                UPDATE(28, 43);
// 29
                UPDATE(29, 30);
                UPDATE(29, 31);
                UPDATE(29, 32);
                UPDATE(29, 33);
                UPDATE(29, 35);
                UPDATE(29, 36);
                UPDATE(29, 37);
                UPDATE(29, 39);
                UPDATE(29, 40);
                UPDATE(29, 42);
// 30
                UPDATE(30, 39);
                UPDATE(30, 40);
                UPDATE(30, 41);
                UPDATE(30, 42);
                UPDATE(30, 43);
                UPDATE(30, 44);
// 31
                UPDATE(31, 36);
                UPDATE(31, 37);
                UPDATE(31, 38);
                UPDATE(31, 42);
                UPDATE(31, 43);
                UPDATE(31, 44);
// 32
                UPDATE(32, 35);
                UPDATE(32, 37);
                UPDATE(32, 38);
                UPDATE(32, 40);
                UPDATE(32, 41);
                UPDATE(32, 44);
// 33
                UPDATE(33, 35);
                UPDATE(33, 36);
                UPDATE(33, 38);
                UPDATE(33, 39);
                UPDATE(33, 41);
                UPDATE(33, 43);
// 34
                UPDATE(34, 35);
                UPDATE(34, 36);
                UPDATE(34, 37);
                UPDATE(34, 39);
                UPDATE(34, 40);
                UPDATE(34, 42);
// 35
                UPDATE(35, 42);
                UPDATE(35, 43);
                UPDATE(35, 44);
// 36
                UPDATE(36, 40);
                UPDATE(36, 41);
                UPDATE(36, 44);
// 37
                UPDATE(37, 39);
                UPDATE(37, 41);
                UPDATE(37, 43);
// 38
                UPDATE(38, 39);
                UPDATE(38, 40);
                UPDATE(38, 42);
// 39
                UPDATE(39, 44);
// 40
                UPDATE(40, 43);
// 41
                UPDATE(41, 42);
// 42
// 43
// 44
            }

            v_s.push_back(edges[best_i].n1 + offset, edges[best_i].n2 + offset);
            v_s.push_back(edges[best_j].n1 + offset, edges[best_j].n2 + offset);
            w_s.push_back(edges[best_i].weight);
            w_s.push_back(edges[best_j].weight);

            return edges[best_i].weight + edges[best_j].weight;
        }

        /**
         * Adds the edge with the maximum weight, but only edges, that do not
         * connect to the 0 node..
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF two_endpoints_skip_first(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            TypeSF max_w = -std::numeric_limits<TypeSF>::max();
            size_t best_i = std::numeric_limits<size_t>::max();

            for (size_t i = 0; i < edges_size; ++i) {
                if (edges[i].weight > max_w && edges[i].n1 != 0 && edges[i].n2 != 0) {
                    ASSERT(edges[i].n1 != 0);
                    ASSERT(edges[i].n2 != 0);

                    max_w = edges[i].weight;
                    best_i = i;
                }
            }

            ASSERT(best_i != std::numeric_limits<size_t>::max());

            v_s.push_back(edges[best_i].n1 + offset, edges[best_i].n2 + offset);
            w_s.push_back(edges[best_i].weight);

            return edges[best_i].weight;

        }

        /**
         * Adds two edges that have the combined maximum weight.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF four_endpoints_skip_first(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            size_t best_i = n_nodes, best_j = n_edges - 1;
            TypeSF best_sum = edges[best_i].weight + edges[best_j].weight;

            for (size_t i = n_nodes; i < n_edges; ++i) {
                for (size_t j = i + 1; j < n_edges; ++j) {
                    if (edges[i].weight + edges[j].weight > best_sum && NO_OVERLAP(i, j)) {
                        best_i = i;
                        best_j = j;
                        best_sum = edges[i].weight + edges[j].weight;
                    }
                }
            }

            v_s.push_back(edges[best_i].n1 + offset, edges[best_i].n2 + offset);
            v_s.push_back(edges[best_j].n1 + offset, edges[best_j].n2 + offset);
            w_s.push_back(edges[best_i].weight);
            w_s.push_back(edges[best_j].weight);

            return best_sum;
        }

        /**
         * Adds two edges that have the combined maximum weight.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF four_endpoints_10_skip_first(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());
            ASSERT(n_nodes >= 4);
            ASSERT(n_nodes <= 10);

            size_t best_i = n_nodes, best_j = n_edges - 1;

            if (n_nodes == 5) {
// 0
// 1
// 2
// 3
// 4
                UPDATE(4, 9);
// 5
                UPDATE(5, 8);
// 6
                UPDATE(6, 7);
// 7
// 8
// 9
            } else if (n_nodes == 6) {
// 0
// 1
// 2
// 3
// 4
// 5
                UPDATE(5, 12);
                UPDATE(5, 13);
                UPDATE(5, 14);
// 6
                UPDATE(6, 10);
                UPDATE(6, 11);
                UPDATE(6, 14);
// 7
                UPDATE(7, 9);
                UPDATE(7, 11);
                UPDATE(7, 13);
// 8
                UPDATE(8, 9);
                UPDATE(8, 10);
                UPDATE(8, 12);
// 9
                UPDATE(9, 14);
// 10
                UPDATE(10, 13);
// 11
                UPDATE(11, 12);
// 12
// 13
// 14
            } else if (n_nodes == 7) {
// 0
// 1
// 2
// 3
// 4
// 5
// 6
                UPDATE(6, 15);
                UPDATE(6, 16);
                UPDATE(6, 17);
                UPDATE(6, 18);
                UPDATE(6, 19);
                UPDATE(6, 20);
// 7
                UPDATE(7, 12);
                UPDATE(7, 13);
                UPDATE(7, 14);
                UPDATE(7, 18);
                UPDATE(7, 19);
                UPDATE(7, 20);
// 8
                UPDATE(8, 11);
                UPDATE(8, 13);
                UPDATE(8, 14);
                UPDATE(8, 16);
                UPDATE(8, 17);
                UPDATE(8, 20);
// 9
                UPDATE(9, 11);
                UPDATE(9, 12);
                UPDATE(9, 14);
                UPDATE(9, 15);
                UPDATE(9, 17);
                UPDATE(9, 19);
// 10
                UPDATE(10, 11);
                UPDATE(10, 12);
                UPDATE(10, 13);
                UPDATE(10, 15);
                UPDATE(10, 16);
                UPDATE(10, 18);
// 11
                UPDATE(11, 18);
                UPDATE(11, 19);
                UPDATE(11, 20);
// 12
                UPDATE(12, 16);
                UPDATE(12, 17);
                UPDATE(12, 20);
// 13
                UPDATE(13, 15);
                UPDATE(13, 17);
                UPDATE(13, 19);
// 14
                UPDATE(14, 15);
                UPDATE(14, 16);
                UPDATE(14, 18);
// 15
                UPDATE(15, 20);
// 16
                UPDATE(16, 19);
// 17
                UPDATE(17, 18);
// 18
// 19
// 20
            } else if (n_nodes == 8) {
// 0
// 1
// 2
// 3
// 4
// 5
// 6
// 7
                UPDATE(7, 18);
                UPDATE(7, 19);
                UPDATE(7, 20);
                UPDATE(7, 21);
                UPDATE(7, 22);
                UPDATE(7, 23);
                UPDATE(7, 24);
                UPDATE(7, 25);
                UPDATE(7, 26);
                UPDATE(7, 27);
// 8
                UPDATE(8, 14);
                UPDATE(8, 15);
                UPDATE(8, 16);
                UPDATE(8, 17);
                UPDATE(8, 22);
                UPDATE(8, 23);
                UPDATE(8, 24);
                UPDATE(8, 25);
                UPDATE(8, 26);
                UPDATE(8, 27);
// 9
                UPDATE(9, 13);
                UPDATE(9, 15);
                UPDATE(9, 16);
                UPDATE(9, 17);
                UPDATE(9, 19);
                UPDATE(9, 20);
                UPDATE(9, 21);
                UPDATE(9, 25);
                UPDATE(9, 26);
                UPDATE(9, 27);
// 10
                UPDATE(10, 13);
                UPDATE(10, 14);
                UPDATE(10, 16);
                UPDATE(10, 17);
                UPDATE(10, 18);
                UPDATE(10, 20);
                UPDATE(10, 21);
                UPDATE(10, 23);
                UPDATE(10, 24);
                UPDATE(10, 27);
// 11
                UPDATE(11, 13);
                UPDATE(11, 14);
                UPDATE(11, 15);
                UPDATE(11, 17);
                UPDATE(11, 18);
                UPDATE(11, 19);
                UPDATE(11, 21);
                UPDATE(11, 22);
                UPDATE(11, 24);
                UPDATE(11, 26);
// 12
                UPDATE(12, 13);
                UPDATE(12, 14);
                UPDATE(12, 15);
                UPDATE(12, 16);
                UPDATE(12, 18);
                UPDATE(12, 19);
                UPDATE(12, 20);
                UPDATE(12, 22);
                UPDATE(12, 23);
                UPDATE(12, 25);
// 13
                UPDATE(13, 22);
                UPDATE(13, 23);
                UPDATE(13, 24);
                UPDATE(13, 25);
                UPDATE(13, 26);
                UPDATE(13, 27);
// 14
                UPDATE(14, 19);
                UPDATE(14, 20);
                UPDATE(14, 21);
                UPDATE(14, 25);
                UPDATE(14, 26);
                UPDATE(14, 27);
// 15
                UPDATE(15, 18);
                UPDATE(15, 20);
                UPDATE(15, 21);
                UPDATE(15, 23);
                UPDATE(15, 24);
                UPDATE(15, 27);
// 16
                UPDATE(16, 18);
                UPDATE(16, 19);
                UPDATE(16, 21);
                UPDATE(16, 22);
                UPDATE(16, 24);
                UPDATE(16, 26);
// 17
                UPDATE(17, 18);
                UPDATE(17, 19);
                UPDATE(17, 20);
                UPDATE(17, 22);
                UPDATE(17, 23);
                UPDATE(17, 25);
// 18
                UPDATE(18, 25);
                UPDATE(18, 26);
                UPDATE(18, 27);
// 19
                UPDATE(19, 23);
                UPDATE(19, 24);
                UPDATE(19, 27);
// 20
                UPDATE(20, 22);
                UPDATE(20, 24);
                UPDATE(20, 26);
// 21
                UPDATE(21, 22);
                UPDATE(21, 23);
                UPDATE(21, 25);
// 22
                UPDATE(22, 27);
// 23
                UPDATE(23, 26);
// 24
                UPDATE(24, 25);
// 25
// 26
// 27
            } else if (n_nodes == 9) {
// 0
// 1
// 2
// 3
// 4
// 5
// 6
// 7
// 8
                UPDATE(8, 21);
                UPDATE(8, 22);
                UPDATE(8, 23);
                UPDATE(8, 24);
                UPDATE(8, 25);
                UPDATE(8, 26);
                UPDATE(8, 27);
                UPDATE(8, 28);
                UPDATE(8, 29);
                UPDATE(8, 30);
                UPDATE(8, 31);
                UPDATE(8, 32);
                UPDATE(8, 33);
                UPDATE(8, 34);
                UPDATE(8, 35);
// 9
                UPDATE(9, 16);
                UPDATE(9, 17);
                UPDATE(9, 18);
                UPDATE(9, 19);
                UPDATE(9, 20);
                UPDATE(9, 26);
                UPDATE(9, 27);
                UPDATE(9, 28);
                UPDATE(9, 29);
                UPDATE(9, 30);
                UPDATE(9, 31);
                UPDATE(9, 32);
                UPDATE(9, 33);
                UPDATE(9, 34);
                UPDATE(9, 35);
// 10
                UPDATE(10, 15);
                UPDATE(10, 17);
                UPDATE(10, 18);
                UPDATE(10, 19);
                UPDATE(10, 20);
                UPDATE(10, 22);
                UPDATE(10, 23);
                UPDATE(10, 24);
                UPDATE(10, 25);
                UPDATE(10, 30);
                UPDATE(10, 31);
                UPDATE(10, 32);
                UPDATE(10, 33);
                UPDATE(10, 34);
                UPDATE(10, 35);
// 11
                UPDATE(11, 15);
                UPDATE(11, 16);
                UPDATE(11, 18);
                UPDATE(11, 19);
                UPDATE(11, 20);
                UPDATE(11, 21);
                UPDATE(11, 23);
                UPDATE(11, 24);
                UPDATE(11, 25);
                UPDATE(11, 27);
                UPDATE(11, 28);
                UPDATE(11, 29);
                UPDATE(11, 33);
                UPDATE(11, 34);
                UPDATE(11, 35);
// 12
                UPDATE(12, 15);
                UPDATE(12, 16);
                UPDATE(12, 17);
                UPDATE(12, 19);
                UPDATE(12, 20);
                UPDATE(12, 21);
                UPDATE(12, 22);
                UPDATE(12, 24);
                UPDATE(12, 25);
                UPDATE(12, 26);
                UPDATE(12, 28);
                UPDATE(12, 29);
                UPDATE(12, 31);
                UPDATE(12, 32);
                UPDATE(12, 35);
// 13
                UPDATE(13, 15);
                UPDATE(13, 16);
                UPDATE(13, 17);
                UPDATE(13, 18);
                UPDATE(13, 20);
                UPDATE(13, 21);
                UPDATE(13, 22);
                UPDATE(13, 23);
                UPDATE(13, 25);
                UPDATE(13, 26);
                UPDATE(13, 27);
                UPDATE(13, 29);
                UPDATE(13, 30);
                UPDATE(13, 32);
                UPDATE(13, 34);
// 14
                UPDATE(14, 15);
                UPDATE(14, 16);
                UPDATE(14, 17);
                UPDATE(14, 18);
                UPDATE(14, 19);
                UPDATE(14, 21);
                UPDATE(14, 22);
                UPDATE(14, 23);
                UPDATE(14, 24);
                UPDATE(14, 26);
                UPDATE(14, 27);
                UPDATE(14, 28);
                UPDATE(14, 30);
                UPDATE(14, 31);
                UPDATE(14, 33);
// 15
                UPDATE(15, 26);
                UPDATE(15, 27);
                UPDATE(15, 28);
                UPDATE(15, 29);
                UPDATE(15, 30);
                UPDATE(15, 31);
                UPDATE(15, 32);
                UPDATE(15, 33);
                UPDATE(15, 34);
                UPDATE(15, 35);
// 16
                UPDATE(16, 22);
                UPDATE(16, 23);
                UPDATE(16, 24);
                UPDATE(16, 25);
                UPDATE(16, 30);
                UPDATE(16, 31);
                UPDATE(16, 32);
                UPDATE(16, 33);
                UPDATE(16, 34);
                UPDATE(16, 35);
// 17
                UPDATE(17, 21);
                UPDATE(17, 23);
                UPDATE(17, 24);
                UPDATE(17, 25);
                UPDATE(17, 27);
                UPDATE(17, 28);
                UPDATE(17, 29);
                UPDATE(17, 33);
                UPDATE(17, 34);
                UPDATE(17, 35);
// 18
                UPDATE(18, 21);
                UPDATE(18, 22);
                UPDATE(18, 24);
                UPDATE(18, 25);
                UPDATE(18, 26);
                UPDATE(18, 28);
                UPDATE(18, 29);
                UPDATE(18, 31);
                UPDATE(18, 32);
                UPDATE(18, 35);
// 19
                UPDATE(19, 21);
                UPDATE(19, 22);
                UPDATE(19, 23);
                UPDATE(19, 25);
                UPDATE(19, 26);
                UPDATE(19, 27);
                UPDATE(19, 29);
                UPDATE(19, 30);
                UPDATE(19, 32);
                UPDATE(19, 34);
// 20
                UPDATE(20, 21);
                UPDATE(20, 22);
                UPDATE(20, 23);
                UPDATE(20, 24);
                UPDATE(20, 26);
                UPDATE(20, 27);
                UPDATE(20, 28);
                UPDATE(20, 30);
                UPDATE(20, 31);
                UPDATE(20, 33);
// 21
                UPDATE(21, 30);
                UPDATE(21, 31);
                UPDATE(21, 32);
                UPDATE(21, 33);
                UPDATE(21, 34);
                UPDATE(21, 35);
// 22
                UPDATE(22, 27);
                UPDATE(22, 28);
                UPDATE(22, 29);
                UPDATE(22, 33);
                UPDATE(22, 34);
                UPDATE(22, 35);
// 23
                UPDATE(23, 26);
                UPDATE(23, 28);
                UPDATE(23, 29);
                UPDATE(23, 31);
                UPDATE(23, 32);
                UPDATE(23, 35);
// 24
                UPDATE(24, 26);
                UPDATE(24, 27);
                UPDATE(24, 29);
                UPDATE(24, 30);
                UPDATE(24, 32);
                UPDATE(24, 34);
// 25
                UPDATE(25, 26);
                UPDATE(25, 27);
                UPDATE(25, 28);
                UPDATE(25, 30);
                UPDATE(25, 31);
                UPDATE(25, 33);
// 26
                UPDATE(26, 33);
                UPDATE(26, 34);
                UPDATE(26, 35);
// 27
                UPDATE(27, 31);
                UPDATE(27, 32);
                UPDATE(27, 35);
// 28
                UPDATE(28, 30);
                UPDATE(28, 32);
                UPDATE(28, 34);
// 29
                UPDATE(29, 30);
                UPDATE(29, 31);
                UPDATE(29, 33);
// 30
                UPDATE(30, 35);
// 31
                UPDATE(31, 34);
// 32
                UPDATE(32, 33);
// 33
// 34
// 35
            } else if (n_nodes == 10) {
// 0
// 1
// 2
// 3
// 4
// 5
// 6
// 7
// 8
// 9
                UPDATE(9, 24);
                UPDATE(9, 25);
                UPDATE(9, 26);
                UPDATE(9, 27);
                UPDATE(9, 28);
                UPDATE(9, 29);
                UPDATE(9, 30);
                UPDATE(9, 31);
                UPDATE(9, 32);
                UPDATE(9, 33);
                UPDATE(9, 34);
                UPDATE(9, 35);
                UPDATE(9, 36);
                UPDATE(9, 37);
                UPDATE(9, 38);
                UPDATE(9, 39);
                UPDATE(9, 40);
                UPDATE(9, 41);
                UPDATE(9, 42);
                UPDATE(9, 43);
                UPDATE(9, 44);
// 10
                UPDATE(10, 18);
                UPDATE(10, 19);
                UPDATE(10, 20);
                UPDATE(10, 21);
                UPDATE(10, 22);
                UPDATE(10, 23);
                UPDATE(10, 30);
                UPDATE(10, 31);
                UPDATE(10, 32);
                UPDATE(10, 33);
                UPDATE(10, 34);
                UPDATE(10, 35);
                UPDATE(10, 36);
                UPDATE(10, 37);
                UPDATE(10, 38);
                UPDATE(10, 39);
                UPDATE(10, 40);
                UPDATE(10, 41);
                UPDATE(10, 42);
                UPDATE(10, 43);
                UPDATE(10, 44);
// 11
                UPDATE(11, 17);
                UPDATE(11, 19);
                UPDATE(11, 20);
                UPDATE(11, 21);
                UPDATE(11, 22);
                UPDATE(11, 23);
                UPDATE(11, 25);
                UPDATE(11, 26);
                UPDATE(11, 27);
                UPDATE(11, 28);
                UPDATE(11, 29);
                UPDATE(11, 35);
                UPDATE(11, 36);
                UPDATE(11, 37);
                UPDATE(11, 38);
                UPDATE(11, 39);
                UPDATE(11, 40);
                UPDATE(11, 41);
                UPDATE(11, 42);
                UPDATE(11, 43);
                UPDATE(11, 44);
// 12
                UPDATE(12, 17);
                UPDATE(12, 18);
                UPDATE(12, 20);
                UPDATE(12, 21);
                UPDATE(12, 22);
                UPDATE(12, 23);
                UPDATE(12, 24);
                UPDATE(12, 26);
                UPDATE(12, 27);
                UPDATE(12, 28);
                UPDATE(12, 29);
                UPDATE(12, 31);
                UPDATE(12, 32);
                UPDATE(12, 33);
                UPDATE(12, 34);
                UPDATE(12, 39);
                UPDATE(12, 40);
                UPDATE(12, 41);
                UPDATE(12, 42);
                UPDATE(12, 43);
                UPDATE(12, 44);
// 13
                UPDATE(13, 17);
                UPDATE(13, 18);
                UPDATE(13, 19);
                UPDATE(13, 21);
                UPDATE(13, 22);
                UPDATE(13, 23);
                UPDATE(13, 24);
                UPDATE(13, 25);
                UPDATE(13, 27);
                UPDATE(13, 28);
                UPDATE(13, 29);
                UPDATE(13, 30);
                UPDATE(13, 32);
                UPDATE(13, 33);
                UPDATE(13, 34);
                UPDATE(13, 36);
                UPDATE(13, 37);
                UPDATE(13, 38);
                UPDATE(13, 42);
                UPDATE(13, 43);
                UPDATE(13, 44);
// 14
                UPDATE(14, 17);
                UPDATE(14, 18);
                UPDATE(14, 19);
                UPDATE(14, 20);
                UPDATE(14, 22);
                UPDATE(14, 23);
                UPDATE(14, 24);
                UPDATE(14, 25);
                UPDATE(14, 26);
                UPDATE(14, 28);
                UPDATE(14, 29);
                UPDATE(14, 30);
                UPDATE(14, 31);
                UPDATE(14, 33);
                UPDATE(14, 34);
                UPDATE(14, 35);
                UPDATE(14, 37);
                UPDATE(14, 38);
                UPDATE(14, 40);
                UPDATE(14, 41);
                UPDATE(14, 44);
// 15
                UPDATE(15, 17);
                UPDATE(15, 18);
                UPDATE(15, 19);
                UPDATE(15, 20);
                UPDATE(15, 21);
                UPDATE(15, 23);
                UPDATE(15, 24);
                UPDATE(15, 25);
                UPDATE(15, 26);
                UPDATE(15, 27);
                UPDATE(15, 29);
                UPDATE(15, 30);
                UPDATE(15, 31);
                UPDATE(15, 32);
                UPDATE(15, 34);
                UPDATE(15, 35);
                UPDATE(15, 36);
                UPDATE(15, 38);
                UPDATE(15, 39);
                UPDATE(15, 41);
                UPDATE(15, 43);
// 16
                UPDATE(16, 17);
                UPDATE(16, 18);
                UPDATE(16, 19);
                UPDATE(16, 20);
                UPDATE(16, 21);
                UPDATE(16, 22);
                UPDATE(16, 24);
                UPDATE(16, 25);
                UPDATE(16, 26);
                UPDATE(16, 27);
                UPDATE(16, 28);
                UPDATE(16, 30);
                UPDATE(16, 31);
                UPDATE(16, 32);
                UPDATE(16, 33);
                UPDATE(16, 35);
                UPDATE(16, 36);
                UPDATE(16, 37);
                UPDATE(16, 39);
                UPDATE(16, 40);
                UPDATE(16, 42);
// 17
                UPDATE(17, 30);
                UPDATE(17, 31);
                UPDATE(17, 32);
                UPDATE(17, 33);
                UPDATE(17, 34);
                UPDATE(17, 35);
                UPDATE(17, 36);
                UPDATE(17, 37);
                UPDATE(17, 38);
                UPDATE(17, 39);
                UPDATE(17, 40);
                UPDATE(17, 41);
                UPDATE(17, 42);
                UPDATE(17, 43);
                UPDATE(17, 44);
// 18
                UPDATE(18, 25);
                UPDATE(18, 26);
                UPDATE(18, 27);
                UPDATE(18, 28);
                UPDATE(18, 29);
                UPDATE(18, 35);
                UPDATE(18, 36);
                UPDATE(18, 37);
                UPDATE(18, 38);
                UPDATE(18, 39);
                UPDATE(18, 40);
                UPDATE(18, 41);
                UPDATE(18, 42);
                UPDATE(18, 43);
                UPDATE(18, 44);
// 19
                UPDATE(19, 24);
                UPDATE(19, 26);
                UPDATE(19, 27);
                UPDATE(19, 28);
                UPDATE(19, 29);
                UPDATE(19, 31);
                UPDATE(19, 32);
                UPDATE(19, 33);
                UPDATE(19, 34);
                UPDATE(19, 39);
                UPDATE(19, 40);
                UPDATE(19, 41);
                UPDATE(19, 42);
                UPDATE(19, 43);
                UPDATE(19, 44);
// 20
                UPDATE(20, 24);
                UPDATE(20, 25);
                UPDATE(20, 27);
                UPDATE(20, 28);
                UPDATE(20, 29);
                UPDATE(20, 30);
                UPDATE(20, 32);
                UPDATE(20, 33);
                UPDATE(20, 34);
                UPDATE(20, 36);
                UPDATE(20, 37);
                UPDATE(20, 38);
                UPDATE(20, 42);
                UPDATE(20, 43);
                UPDATE(20, 44);
// 21
                UPDATE(21, 24);
                UPDATE(21, 25);
                UPDATE(21, 26);
                UPDATE(21, 28);
                UPDATE(21, 29);
                UPDATE(21, 30);
                UPDATE(21, 31);
                UPDATE(21, 33);
                UPDATE(21, 34);
                UPDATE(21, 35);
                UPDATE(21, 37);
                UPDATE(21, 38);
                UPDATE(21, 40);
                UPDATE(21, 41);
                UPDATE(21, 44);
// 22
                UPDATE(22, 24);
                UPDATE(22, 25);
                UPDATE(22, 26);
                UPDATE(22, 27);
                UPDATE(22, 29);
                UPDATE(22, 30);
                UPDATE(22, 31);
                UPDATE(22, 32);
                UPDATE(22, 34);
                UPDATE(22, 35);
                UPDATE(22, 36);
                UPDATE(22, 38);
                UPDATE(22, 39);
                UPDATE(22, 41);
                UPDATE(22, 43);
// 23
                UPDATE(23, 24);
                UPDATE(23, 25);
                UPDATE(23, 26);
                UPDATE(23, 27);
                UPDATE(23, 28);
                UPDATE(23, 30);
                UPDATE(23, 31);
                UPDATE(23, 32);
                UPDATE(23, 33);
                UPDATE(23, 35);
                UPDATE(23, 36);
                UPDATE(23, 37);
                UPDATE(23, 39);
                UPDATE(23, 40);
                UPDATE(23, 42);
// 24
                UPDATE(24, 35);
                UPDATE(24, 36);
                UPDATE(24, 37);
                UPDATE(24, 38);
                UPDATE(24, 39);
                UPDATE(24, 40);
                UPDATE(24, 41);
                UPDATE(24, 42);
                UPDATE(24, 43);
                UPDATE(24, 44);
// 25
                UPDATE(25, 31);
                UPDATE(25, 32);
                UPDATE(25, 33);
                UPDATE(25, 34);
                UPDATE(25, 39);
                UPDATE(25, 40);
                UPDATE(25, 41);
                UPDATE(25, 42);
                UPDATE(25, 43);
                UPDATE(25, 44);
// 26
                UPDATE(26, 30);
                UPDATE(26, 32);
                UPDATE(26, 33);
                UPDATE(26, 34);
                UPDATE(26, 36);
                UPDATE(26, 37);
                UPDATE(26, 38);
                UPDATE(26, 42);
                UPDATE(26, 43);
                UPDATE(26, 44);
// 27
                UPDATE(27, 30);
                UPDATE(27, 31);
                UPDATE(27, 33);
                UPDATE(27, 34);
                UPDATE(27, 35);
                UPDATE(27, 37);
                UPDATE(27, 38);
                UPDATE(27, 40);
                UPDATE(27, 41);
                UPDATE(27, 44);
// 28
                UPDATE(28, 30);
                UPDATE(28, 31);
                UPDATE(28, 32);
                UPDATE(28, 34);
                UPDATE(28, 35);
                UPDATE(28, 36);
                UPDATE(28, 38);
                UPDATE(28, 39);
                UPDATE(28, 41);
                UPDATE(28, 43);
// 29
                UPDATE(29, 30);
                UPDATE(29, 31);
                UPDATE(29, 32);
                UPDATE(29, 33);
                UPDATE(29, 35);
                UPDATE(29, 36);
                UPDATE(29, 37);
                UPDATE(29, 39);
                UPDATE(29, 40);
                UPDATE(29, 42);
// 30
                UPDATE(30, 39);
                UPDATE(30, 40);
                UPDATE(30, 41);
                UPDATE(30, 42);
                UPDATE(30, 43);
                UPDATE(30, 44);
// 31
                UPDATE(31, 36);
                UPDATE(31, 37);
                UPDATE(31, 38);
                UPDATE(31, 42);
                UPDATE(31, 43);
                UPDATE(31, 44);
// 32
                UPDATE(32, 35);
                UPDATE(32, 37);
                UPDATE(32, 38);
                UPDATE(32, 40);
                UPDATE(32, 41);
                UPDATE(32, 44);
// 33
                UPDATE(33, 35);
                UPDATE(33, 36);
                UPDATE(33, 38);
                UPDATE(33, 39);
                UPDATE(33, 41);
                UPDATE(33, 43);
// 34
                UPDATE(34, 35);
                UPDATE(34, 36);
                UPDATE(34, 37);
                UPDATE(34, 39);
                UPDATE(34, 40);
                UPDATE(34, 42);
// 35
                UPDATE(35, 42);
                UPDATE(35, 43);
                UPDATE(35, 44);
// 36
                UPDATE(36, 40);
                UPDATE(36, 41);
                UPDATE(36, 44);
// 37
                UPDATE(37, 39);
                UPDATE(37, 41);
                UPDATE(37, 43);
// 38
                UPDATE(38, 39);
                UPDATE(38, 40);
                UPDATE(38, 42);
// 39
                UPDATE(39, 44);
// 40
                UPDATE(40, 43);
// 41
                UPDATE(41, 42);
// 42
// 43
// 44
            }

            v_s.push_back(edges[best_i].n1 + offset, edges[best_i].n2 + offset);
            v_s.push_back(edges[best_j].n1 + offset, edges[best_j].n2 + offset);
            w_s.push_back(edges[best_i].weight);
            w_s.push_back(edges[best_j].weight);

            return edges[best_i].weight + edges[best_j].weight;
        }

        /**
         * Calculate an upper bound for a subset optimization problem using a greedy algorithm.
         *
         * This function calculates an upper bound for a subset optimization problem using a greedy algorithm.
         * It sorts the edges based on their weights in descending order and adds edges to the solution until
         * the number of included edges reaches half of the total number of endpoints.
         *
         * @param v_s The vector of vectors containing subsets of endpoints.
         * @param w_s The vector of weights corresponding to the subsets.
         * @param n_endpoints The number of endpoints.
         * @param offset An offset value used in the calculation.
         * @return The calculated upper bound using the greedy algorithm.
         */
        TypeSF get_upper_bound_greedy(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            if (!greedy_is_sorted) {
                greedy_edges.resize(n_edges);
                for (size_t i = 0; i < n_edges; ++i) {
                    greedy_edges[i] = edges[i];
                }

                std::sort(greedy_edges.begin(), greedy_edges.end(), [](const Edge<TypeSF> &x, const Edge<TypeSF> &y) { return x.weight > y.weight; });
                greedy_is_sorted = true;
            }

            size_t n_pairs_needed = n_endpoints / 2;
            TypeSF sum = 0;
            for (size_t i = 0; i < n_pairs_needed; ++i) {
                v_s.push_back(greedy_edges[i].n1 + offset, greedy_edges[i].n2 + offset);
                w_s.push_back(greedy_edges[i].weight);
                sum += greedy_edges[i].weight;
            }

            return sum;
        }

        /**
         * Computes an upper bound for a maximum weight matching, excluding edges connected to vertex with ID 0,
         * using a greedy algorithm.
         *
         * This function sorts the edges in descending order of weight and selects edges until reaching
         * n_endpoints/2, adding them to v_s and w_s. Edges connected to vertex with ID 0 are skipped.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of weights of selected edges.
         */
        TypeSF get_upper_bound_skip_first_greedy(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            if (!greedy_is_sorted) {
                greedy_edges.resize(n_edges);
                for (size_t i = 0; i < n_edges; ++i) {
                    greedy_edges[i] = edges[i];
                }

                std::sort(greedy_edges.begin(), greedy_edges.end(), [](const Edge<TypeSF> &x, const Edge<TypeSF> &y) { return x.weight > y.weight; });
                greedy_is_sorted = true;
            }

            size_t n_pairs_needed = n_endpoints / 2;
            TypeSF sum = 0;
            size_t n_pairs_chosen = 0;
            size_t i = 0;
            while (n_pairs_chosen < n_pairs_needed) {
                if (greedy_edges[i].n1 == 0 || greedy_edges[i].n2 == 0) {
                    ++i;
                    continue;
                }

                v_s.push_back(greedy_edges[i].n1 + offset, greedy_edges[i].n2 + offset);
                w_s.push_back(greedy_edges[i].weight);
                sum += greedy_edges[i].weight;

                ++i;
                ++n_pairs_chosen;
            }

            return sum;
        }

        /**
         * Gets a maximum weight matching.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights.
         */
        TypeSF get_upper_bound_matching(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            size_t dummy_n_nodes = n_nodes;
            size_t dummy_n_edges = dummy_n_nodes * n_nodes + ((1 + dummy_n_nodes) / 2);
            size_t total_n_nodes = n_nodes + dummy_n_nodes;
            size_t total_n_edges = n_edges + dummy_n_edges;

            // create the graph if not already created
            if (pm == nullptr) {
                pm_last_n_endpoints = 0;

                // create the structure
                pm = new PerfectMatching((int) total_n_nodes, (int) total_n_edges);
                pm->options = options;

                // add the clique
                for (size_t i = 0; i < n_edges; ++i) {
                    pm->AddEdge(edges[i].n1, edges[i].n2, -(double) edges[i].weight);
                }

                // add all fake edges
                for (size_t i = 0; i < dummy_n_nodes; ++i) {
                    for (size_t j = 0; j < n_nodes; ++j) {
                        pm->AddEdge(n_nodes + i, (int) j, -(double) (max_weight + 1));
                    }
                }

                // add the edges normally
                for (size_t i = 0; i < n_endpoints / 2; ++i) {
                    pm->AddEdge(n_nodes + (2 * i), n_nodes + (2 * i) + 1, -(double) (2 * max_weight + 3));
                }
            } else {
                // update the graph
                pm->StartUpdate();
                for (size_t i = pm_last_n_endpoints / 2; i < n_endpoints / 2; ++i) {
                    pm->AddNewEdge(n_nodes + (2 * i), n_nodes + (2 * i) + 1, -(double) (2 * max_weight + 3), false);
                }
                pm->FinishUpdate();
            }
            pm->Solve();
            pm_last_n_endpoints = n_endpoints;

            // get solution
            TypeSF sum = 0;
            size_t j = 0;
            for (size_t i = 0; i < n_nodes; ++i) {
                size_t match_i = pm->GetMatch((int) i);
                if (match_i >= n_nodes || match_i < i) {
                    continue;
                }

                while (edges[j].n1 != i) {
                    ++j;
                }
                while (edges[j].n2 != match_i) {
                    ++j;
                }

                v_s.push_back(i + offset, match_i + offset);
                w_s.push_back(edges[j].weight);

                sum += edges[j].weight;
            }

            return sum;
        }

        /**
         * Gets a maximum weight matching, but skips vertex with id 0 and all its
         * edges.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights.
         */
        TypeSF get_upper_bound_skip_first_matching(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            size_t clique_n_nodes = n_nodes - 1;
            size_t clique_n_edges = (clique_n_nodes * (clique_n_nodes - 1)) / 2;
            size_t dummy_n_nodes = clique_n_nodes;
            size_t dummy_n_edges = dummy_n_nodes * clique_n_nodes + ((1 + dummy_n_nodes) / 2);
            size_t total_n_nodes = clique_n_nodes + dummy_n_nodes;
            size_t total_n_edges = clique_n_edges + dummy_n_edges;

            if (pm_skip_first == nullptr) {
                pm_skip_first_last_n_endpoints = 0;

                pm_skip_first = new PerfectMatching((int) total_n_nodes, (int) total_n_edges);
                pm_skip_first->options = options;

                // add normal edges to boost graph
                for (size_t i = 0; i < n_edges; ++i) {
                    size_t node_1 = edges[i].n1;
                    size_t node_2 = edges[i].n2;

                    if (node_1 == 0 || node_2 == 0) {
                        // skip edges with node_id 0
                        continue;
                    }

                    // ID 0 in this structure is not used, but we have to use it in the
                    // 3rd party structure. We shift each node 1 to the left, but later
                    // have to shift them back to normal
                    size_t node_1_shifted = node_1 - 1;
                    size_t node_2_shifted = node_2 - 1;
                    pm_skip_first->AddEdge((int) node_1_shifted, (int) node_2_shifted, -(double) edges[i].weight);
                }

                // add fake edges
                for (size_t i = 0; i < dummy_n_nodes; ++i) {
                    for (size_t j = 0; j < n_nodes; ++j) {
                        size_t node_1 = clique_n_nodes + i; // fake nodes after clique
                        size_t node_2 = j;

                        if (node_2 == 0) {
                            // skip edges with node_id 0
                            continue;
                        }

                        size_t node_1_shifted = node_1; // no shift, because fake node
                        size_t node_2_shifted = node_2 - 1;
                        pm_skip_first->AddEdge((int) node_1_shifted, (int) node_2_shifted, -(double) (max_weight + 1));
                    }
                }

                // add the edges normally
                for (size_t i = 0; i < n_endpoints / 2; ++i) {
                    pm_skip_first->AddEdge((int) clique_n_nodes + (2 * i), (int) clique_n_nodes + (2 * i) + 1, -(double) (2 * max_weight + 3));
                }
            } else {
                // update the graph
                pm_skip_first->StartUpdate();
                for (size_t i = pm_skip_first_last_n_endpoints / 2; i < n_endpoints / 2; ++i) {
                    pm_skip_first->AddNewEdge((int) clique_n_nodes + (2 * i), (int) clique_n_nodes + (2 * i) + 1, -(double) (2 * max_weight + 3), false);
                }
                pm_skip_first->FinishUpdate();
            }

            pm_skip_first->Solve();
            pm_skip_first_last_n_endpoints = n_endpoints;

            // get solution
            TypeSF sum = 0;
            size_t j = 0;
            for (size_t i = 0; i < clique_n_nodes; ++i) {
                size_t match_i = pm_skip_first->GetMatch((int) i);
                if (match_i >= clique_n_nodes || match_i < i) {
                    continue;
                }

                while (edges[j].n1 != i + 1) {
                    ++j;
                }
                while (edges[j].n2 != match_i + 1) {
                    ++j;
                }
                v_s.push_back(i + 1 + offset, match_i + 1 + offset);
                w_s.push_back(edges[j].weight);
                sum += edges[j].weight;
            }

            return sum;
        }

        /**
         * Calculates an upper bound for a given set of vertices and edges using a brute-force
         * iterative approach.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF get_upper_bound_bforce_iterative(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            TypeSF best_sum = -std::numeric_limits<TypeSF>::max();
            size_t n_to_choose = n_endpoints / 2;

            best_set.resize(n_to_choose);

            set.resize(n_to_choose);
            std::iota(set.begin(), set.end(), 0);
            set[n_to_choose - 1] -= 1;

            count.resize(n_nodes, 0);

            // brute force iterate over all possible sets
            while (next_subset(set, n_to_choose, n_edges)) {
                if (valid_set(set, count, 0)) {
                    TypeSF sum = 0;
                    for (size_t i = 0; i < n_to_choose; ++i) {
                        sum += edges[set[i]].weight;
                    }

                    if (sum > best_sum) {
                        best_sum = sum;
                        std::copy(set.begin(), set.end(), best_set.begin());
                    }
                }
            }

            // Add the best-set edges to v_s and w_s
            for (size_t i = 0; i < n_to_choose; ++i) {
                size_t idx1 = offset + edges[best_set[i]].n1;
                size_t idx2 = offset + edges[best_set[i]].n2;
                v_s.push_back(idx1, idx2);
                w_s.push_back(edges[best_set[i]].weight);
            }

            return best_sum;
        }

        /**
         * Calculates an upper bound for a maximum weight matching by iterating through
         * precomputed valid sets of edges. This version of the function uses a lookup
         * table for valid sets and selects the subset with the maximum total weight.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF get_upper_bound_bforce_lookup(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            TypeSF best_sum = -std::numeric_limits<TypeSF>::max();
            size_t n_to_choose = n_endpoints / 2;

            best_set.resize(n_to_choose);
            // Lookup precomputed valid sets of edges
            const std::vector<size_t> &sets = valid_sets[n_nodes - 2][n_to_choose - 1];
            const size_t n_sets = sets.size() / n_to_choose;

            // Find the subset with the maximum total weight
            for (size_t i = 0; i < n_sets; ++i) {
                TypeSF sum = 0;
                for (size_t j = 0; j < n_to_choose; ++j) {
                    sum += edges[sets[i * n_to_choose + j]].weight;
                }

                if (sum > best_sum) {
                    best_sum = sum;
                    for (size_t j = 0; j < n_to_choose; ++j) {
                        best_set[j] = sets[i * n_to_choose + j];
                    }
                }
            }

            // Populate v_s and w_s with the selected edges
            for (size_t i = 0; i < n_to_choose; ++i) {
                size_t idx1 = offset + edges[best_set[i]].n1;
                size_t idx2 = offset + edges[best_set[i]].n2;
                v_s.push_back(idx1, idx2);
                w_s.push_back(edges[best_set[i]].weight);
            }

            return best_sum;
        }

        /**
         * Calculates an upper bound for a maximum weight matching by selecting the best
         * subset of edges iteratively while skipping the first vertex and its edges.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF get_upper_bound_skip_first_bforce_iterative(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            TypeSF best_sum = -std::numeric_limits<TypeSF>::max();
            size_t n_to_choose = n_endpoints / 2;

            best_set.resize(n_to_choose);

            set.resize(n_to_choose);
            std::iota(set.begin(), set.end(), 0);
            set[n_to_choose - 1] -= 1;

            count.clear();
            count.resize(n_nodes, 0);

            // Brute force iterate over all possible sets, excluding the first vertex and its edges
            while (next_subset(set, n_to_choose, n_edges - n_nodes - 1)) {
                if (valid_set(set, count, n_nodes - 1)) {
                    TypeSF sum = 0;
                    for (size_t i = 0; i < n_to_choose; ++i) {
                        sum += edges[n_nodes - 1 + set[i]].weight;
                    }

                    if (sum > best_sum) {
                        best_sum = sum;
                        std::copy(set.begin(), set.end(), best_set.begin());
                    }
                }
            }

            for (size_t i = 0; i < n_to_choose; ++i) {
                size_t idx1 = offset + edges[n_nodes - 1 + best_set[i]].n1;
                size_t idx2 = offset + edges[n_nodes - 1 + best_set[i]].n2;
                v_s.push_back(idx1, idx2);
                w_s.push_back(edges[n_nodes - 1 + best_set[i]].weight);
            }

            return best_sum;
        }

        /**
         * Calculates an upper bound for a maximum weight matching by selecting the best
         * subset of edges based on precomputed valid sets while skipping the first vertex
         * and its edges.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF get_upper_bound_skip_first_bforce_lookup(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            TypeSF best_sum = -std::numeric_limits<TypeSF>::max();
            size_t n_to_choose = n_endpoints / 2;

            best_set.resize(n_to_choose);

            const std::vector<size_t> &sets = valid_sets[n_nodes - 3][n_to_choose - 1];
            const size_t n_sets = sets.size() / n_to_choose;

            for (size_t i = 0; i < n_sets; ++i) {
                TypeSF sum = 0;
                for (size_t j = 0; j < n_to_choose; ++j) {
                    sum += edges[n_nodes - 1 + sets[i * n_to_choose + j]].weight;
                }

                if (sum > best_sum) {
                    best_sum = sum;
                    for (size_t j = 0; j < n_to_choose; ++j) {
                        best_set[j] = sets[i * n_to_choose + j];
                    }
                }
            }

            for (size_t i = 0; i < n_to_choose; ++i) {
                size_t idx1 = offset + edges[n_nodes - 1 + best_set[i]].n1;
                size_t idx2 = offset + edges[n_nodes - 1 + best_set[i]].n2;
                v_s.push_back(idx1, idx2);
                w_s.push_back(edges[n_nodes - 1 + best_set[i]].weight);
            }

            return best_sum;
        }

        /**
         * Initialize a dynamic array for caching edge weights and indices based on bitset representation.
         */
        void initialize_dynamic_array() {
            // Return if already initialized
            if (solution_ready[0] == 1) {
                return;
            }

            // Populate caches with edge data
            for (size_t i = 0; i < edges.size(); ++i) {
                size_t bitset = 0 | (1 << edges[i].n1) | (1 << edges[i].n2);
                cache_si[bitset] = edges[i].weight;
                cache_bitset_t1[bitset] = i;
                cache_bitset_t2[bitset] = i;
            }

            // Mark initialization as complete
            solution_ready[0] = 1;
        }

        /**
         * Update the dynamic array for caching edge weights and indices.
         *
         * @param idx The index to update.
         */
        void update_dynamic_array(size_t idx) {
            // Return if already initialized
            if (solution_ready[idx] == 1) {
                return;
            }

            // If idx is 0, initialize the dynamic array and return
            if (idx == 0) {
                initialize_dynamic_array();
                return;
            }

            // Recursively update the dynamic array for the previous index
            update_dynamic_array(idx - 1);

            size_t set_size = 2 + idx * 2;
            size_t last_bitset = (pow_2(set_size) - 1) << (n_nodes - set_size);
            size_t bitset0 = pow_2(set_size) - 1;

            // iterate over all node subsets of size set_size
            while (bitset0 <= last_bitset) {
                // determine bitset for the node set
                for (size_t i = 0; i < n_nodes; ++i) {
                    if (!((bitset0 >> i) & 1)) {
                        continue;
                    }
                    for (size_t j = i + 1; j < n_nodes; ++j) {
                        if (!((bitset0 >> j) & 1)) {
                            continue;
                        }
                        // split node subset into two bitsets, one for T[0] and one for T[idx-1]
                        size_t bitset1 = (1ULL << i) | (1LL << j); // set i and j to true
                        size_t bitset2 = bitset0 & ~(1ULL << i) & ~(1ULL << j); // set i and j to false

                        // check for the new minimum
                        if (cache_si[bitset1] + cache_si[bitset2] < cache_si[bitset0]) {
                            cache_si[bitset0] = cache_si[bitset1] + cache_si[bitset2];
                            cache_bitset_t1[bitset0] = bitset1;
                            cache_bitset_t2[bitset0] = bitset2;
                        }
                    }
                }
                bitset0 = next_subset(bitset0, (int) set_size);
            }
            solution_ready[idx] = 1;
        }

        /**
         * Get the upper bound for a maximum weight matching using dynamic programming.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF get_upper_bound_dynamic(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            size_t idx = (n_endpoints / 2) - 1;
            update_dynamic_array(idx);

            // iterate over all node subsets of size n_endpoints
            size_t set_size = 2 + idx * 2;
            size_t last_bitset = (pow_2(set_size) - 1) << (n_nodes - set_size);

            TypeSF best_sum = -std::numeric_limits<TypeSF>::max();
            size_t best_bitset = 0;
            size_t bitset0 = pow_2(set_size) - 1;

            while (bitset0 <= last_bitset) {
                // find maximum si
                if (cache_si[bitset0] > best_sum) {
                    best_sum = cache_si[bitset0];
                    best_bitset = bitset0;
                }
                bitset0 = next_subset(bitset0, (int) set_size);
            }

            // get the pairs from best_edge_idx
            size_t edge_idx;
            size_t bitset_temp = best_bitset;
            for (size_t i = 0; i < (n_endpoints / 2) - 2; ++i) {
                edge_idx = cache_bitset_t1[bitset_temp];
                bitset_temp = cache_bitset_t2[bitset_temp];

                edge_idx = cache_bitset_t1[edge_idx];
                v_s.push_back(edges[edge_idx].n1 + offset, edges[edge_idx].n2 + offset);
                w_s.push_back(edges[edge_idx].weight);
            }
            edge_idx = cache_bitset_t1[bitset_temp];
            bitset_temp = cache_bitset_t2[bitset_temp];

            edge_idx = cache_bitset_t1[edge_idx];
            v_s.push_back(edges[edge_idx].n1 + offset, edges[edge_idx].n2 + offset);
            w_s.push_back(edges[edge_idx].weight);

            edge_idx = cache_bitset_t2[bitset_temp];
            v_s.push_back(edges[edge_idx].n1 + offset, edges[edge_idx].n2 + offset);
            w_s.push_back(edges[edge_idx].weight);

            return best_sum;
        }

        /**
         * Get the upper bound for a maximum weight matching using dynamic programming,
         * but it will skip the node with index 0.
         *
         * @param v_s Vector containing vertices representing edges (two per edge).
         * @param w_s Vector storing weights for each edge.
         * @param n_endpoints Number of endpoints in the clique (selecting n_endpoints / 2 edges).
         * @param offset Offset for vertex conversion.
         * @return The sum of the weights of the edges in the subset that provides the upper bound.
         */
        TypeSF get_upper_bound_skip_first_dynamic(VectorOfVectors<size_t> &v_s, std::vector<TypeSF> &w_s, size_t n_endpoints, size_t offset) {
            ASSERT(v_s.size == 0);
            ASSERT(w_s.empty());

            size_t idx = (n_endpoints / 2) - 1;
            update_dynamic_array(idx);

            size_t set_size = 2 + idx * 2;
            size_t last_bitset = (pow_2(set_size) - 1) << (n_nodes - 1 - set_size);

            TypeSF best_sum = -std::numeric_limits<TypeSF>::max();
            size_t best_bitset = 0;
            size_t bitset0 = pow_2(set_size) - 1;

            while (bitset0 <= last_bitset) {
                // find maximum si
                size_t temp = bitset0 << 1;
                if (cache_si[temp] > best_sum) {
                    best_sum = cache_si[temp];
                    best_bitset = temp;
                }
                bitset0 = next_subset(bitset0, (int) set_size);
            }

            // get the pairs from best_edge_idx
            size_t edge_idx;
            size_t bitset_temp = best_bitset;
            for (size_t i = 0; i < (n_endpoints / 2) - 2; ++i) {
                edge_idx = cache_bitset_t1[bitset_temp];
                bitset_temp = cache_bitset_t2[bitset_temp];

                edge_idx = cache_bitset_t1[edge_idx];
                v_s.push_back(edges[edge_idx].n1 + offset, edges[edge_idx].n2 + offset);
                w_s.push_back(edges[edge_idx].weight);
            }
            edge_idx = cache_bitset_t1[bitset_temp];
            bitset_temp = cache_bitset_t2[bitset_temp];

            edge_idx = cache_bitset_t1[edge_idx];
            v_s.push_back(edges[edge_idx].n1 + offset, edges[edge_idx].n2 + offset);
            w_s.push_back(edges[edge_idx].weight);

            edge_idx = cache_bitset_t2[bitset_temp];
            v_s.push_back(edges[edge_idx].n1 + offset, edges[edge_idx].n2 + offset);
            w_s.push_back(edges[edge_idx].weight);

            return best_sum;
        }

        /**
         * Checks if the set contains no overlapping pairs.
         *
         * @param v The set.
         * @return True if it is a valid set.
         */
        bool valid_set(std::vector<size_t> &v, std::vector<size_t> &c, size_t offset) {
            // v is the set and c is the counter

            // reset count
            for (auto &i: c) {
                i = 0;
            }

            // count
            for (auto &i: v) {
                c[edges[offset + i].n1] += 1;
                c[edges[offset + i].n2] += 1;
            }

            // check if any has more than 2 counts
            bool valid = true;
            for (auto &i: c) {
                valid &= i <= 1;
            }

            return valid;
        }
    };

}

#endif //SUBSETOPTIMIZATION_UB2DALGORITHM_H
