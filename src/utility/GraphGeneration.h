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

#ifndef SUBSETOPTIMIZATION_GRAPHGENERATION_H
#define SUBSETOPTIMIZATION_GRAPHGENERATION_H

#include <cstdint>
#include <string>
#include <vector>
#include <filesystem>

#include "Util.h"
#include "../structures/Graph.h"
#include "../structures/GraphNegativeGroupFarness.h"
#include "../structures/GraphPartialDominatingSet.h"

namespace SMSM {

    /**
     * Generates a random graph.
     *
     * @param n_nodes Number of nodes.
     * @param n_edges Number of edges.
     * @return The graph.
     */
    GraphNegativeGroupFarness<int> random_GraphNegativeGroupFarness(size_t n_nodes, size_t n_edges);

    /**
     * Generates a random graph dataset.
     *
     * @param n_nodes Number of nodes.
     * @param n_graphs Number of graphs.
     * @return Vector holding the graphs.
     */
    std::vector<GraphNegativeGroupFarness<int>> random_GraphNegativeGroupFarness_dataset(size_t n_nodes, size_t n_graphs);

    /**
     * Generates a random graph.
     *
     * @param n_nodes Number of nodes.
     * @param n_edges Number of edges.
     * @return The graph.
     */
    GraphPartialDominatingSet<int> random_GraphPartialDominatingSet(size_t n_nodes, size_t n_edges);

    /**
     * Generates a random graph dataset.
     *
     * @param n_nodes Number of nodes.
     * @param n_graphs Number of graphs.
     * @return Vector holding the graphs.
     */
    std::vector<GraphPartialDominatingSet<int>> random_GraphPartialDominatingSet_dataset(size_t n_nodes, size_t n_graphs);

}

#endif //SUBSETOPTIMIZATION_GRAPHGENERATION_H
