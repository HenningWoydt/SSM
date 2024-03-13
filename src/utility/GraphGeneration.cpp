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

#include "GraphGeneration.h"

namespace SMSM {

    GraphNegativeGroupFarness<int> random_GraphNegativeGroupFarness(size_t n_nodes, size_t n_edges) {
        GraphNegativeGroupFarness<int> g(n_nodes, n_edges);

        size_t curr_n_edges = 0;
        std::vector<std::vector<uint32_t>> forrest;

        // forrest has n_nodes trees.
        for (uint32_t node = 0; node < n_nodes; ++node) {
            forrest.push_back({node});
        }

        // connect two trees (n_node-1) times
        for (uint32_t i = 0; i < n_nodes - 1; ++i) {
            size_t tree1_idx = random_int() % forrest.size();
            size_t tree2_idx = random_int() % forrest.size();

            while (tree1_idx == tree2_idx) {
                tree2_idx = random_int() % forrest.size();
            }

            uint32_t node1 = forrest[tree1_idx][random_int() % forrest[tree1_idx].size()];
            uint32_t node2 = forrest[tree2_idx][random_int() % forrest[tree2_idx].size()];
            g.add_edge(node1, node2);
            curr_n_edges += 1;

            forrest[tree1_idx].insert(forrest[tree1_idx].end(), forrest[tree2_idx].begin(), forrest[tree2_idx].end());
            forrest.erase(forrest.begin() + tree2_idx);
        }
        // one large tree

        std::vector<uint32_t> available_n1;
        std::vector<uint32_t> available_n2;
        for (uint32_t n1 = 0; n1 < n_nodes; ++n1) {
            for (uint32_t n2 = n1 + 1; n2 < n_nodes; ++n2) {
                if (!g.edge_exists(n1, n2)) {
                    available_n1.push_back(n1);
                    available_n2.push_back(n2);
                }
            }
        }

        // connect random graphs
        for (; curr_n_edges < n_edges; ++curr_n_edges) {
            size_t idx = random_int() % available_n1.size();
            g.add_edge(available_n1[idx], available_n2[idx]);

            std::swap(available_n1[idx], available_n1.back());
            std::swap(available_n2[idx], available_n2.back());
            available_n1.pop_back();
            available_n2.pop_back();
        }

        g.finalize();
        return g;
    }

    std::vector<GraphNegativeGroupFarness<int>> random_GraphNegativeGroupFarness_dataset(size_t n_nodes, size_t n_graphs) {
        std::vector<GraphNegativeGroupFarness<int>> vec;

        for (uint32_t n = 0; n < n_graphs; ++n) {
            uint32_t n_edges = (n_nodes - 1) + (random_int() % ((n_nodes * (n_nodes - 1) / 2) - (n_nodes - 1)));
            vec.emplace_back(random_GraphNegativeGroupFarness(n_nodes, n_edges));
        }

        return vec;
    }

    GraphPartialDominatingSet<int> random_GraphPartialDominatingSet(size_t n_nodes, size_t n_edges) {
        GraphPartialDominatingSet<int> g(n_nodes, n_edges);

        size_t curr_n_edges = 0;
        std::vector<std::vector<uint32_t>> forrest;

        // forrest has n_nodes trees.
        for (uint32_t node = 0; node < n_nodes; ++node) {
            forrest.push_back({node});
        }

        // connect two trees (n_node-1) times
        for (uint32_t i = 0; i < n_nodes - 1; ++i) {
            size_t tree1_idx = random_int() % forrest.size();
            size_t tree2_idx = random_int() % forrest.size();

            while (tree1_idx == tree2_idx) {
                tree2_idx = random_int() % forrest.size();
            }

            uint32_t node1 = forrest[tree1_idx][random_int() % forrest[tree1_idx].size()];
            uint32_t node2 = forrest[tree2_idx][random_int() % forrest[tree2_idx].size()];
            g.add_edge(node1, node2);
            curr_n_edges += 1;

            forrest[tree1_idx].insert(forrest[tree1_idx].end(), forrest[tree2_idx].begin(), forrest[tree2_idx].end());
            forrest.erase(forrest.begin() + tree2_idx);
        }
        // one large tree

        std::vector<uint32_t> available_n1;
        std::vector<uint32_t> available_n2;
        for (uint32_t n1 = 0; n1 < n_nodes; ++n1) {
            for (uint32_t n2 = n1 + 1; n2 < n_nodes; ++n2) {
                if (!g.edge_exists(n1, n2)) {
                    available_n1.push_back(n1);
                    available_n2.push_back(n2);
                }
            }
        }

        // connect random graphs
        for (; curr_n_edges < n_edges; ++curr_n_edges) {
            size_t idx = random_int() % available_n1.size();
            g.add_edge(available_n1[idx], available_n2[idx]);

            std::swap(available_n1[idx], available_n1.back());
            std::swap(available_n2[idx], available_n2.back());
            available_n1.pop_back();
            available_n2.pop_back();
        }

        g.finalize();
        return g;
    }

    std::vector<GraphPartialDominatingSet<int>> random_GraphPartialDominatingSet_dataset(size_t n_nodes, size_t n_graphs) {
        std::vector<GraphPartialDominatingSet<int>> vec;

        for (uint32_t n = 0; n < n_graphs; ++n) {
            uint32_t n_edges = (n_nodes - 1) + (random_int() % ((n_nodes * (n_nodes - 1) / 2) - (n_nodes - 1)));
            vec.emplace_back(random_GraphPartialDominatingSet(n_nodes, n_edges));
        }

        return vec;
    }

}
