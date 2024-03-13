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

#ifndef SUBSETOPTIMIZATION_GRAPH_H
#define SUBSETOPTIMIZATION_GRAPH_H

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <optional>
#include <vector>
#include <fstream>
#include <cstring>

#include "../utility/VectorOfVectors.h"
#include "../utility/Util.h"
#include "StructureInterface.h"

namespace SMSM {

/**
 * Class for defining a graph.
 * */
    template<typename TypeSF>
    class Graph : public StructureInterface<TypeSF> {
    public:
        // normal variables to model the graph
        size_t n_nodes = 0;
        size_t n_edges = 0;
        std::vector<std::vector<uint32_t>> adj_list;

    public:
        /**
         * Initializes a graph.
         *
         * @return A graph object.
         * */
        Graph() = default;

        /**
         * Initializes a graph.
         *
         * @return A graph object.
         */
        Graph(size_t n, size_t m) {
            reserve(n, m);
        };

        /**
        * Reads a graph from the specified file. File should have the following format:
        * idx1 idx2\n
        * idx3 idx4\n
        * ...
        *
        * Vertices have values from 0 to n-1.
        * Lines with '%' are ignored.
        *
        * @param file_path Path to the file.
        * @return The Graph.
        */
        explicit Graph(const std::string &file_path) {
            if (!file_exists(file_path)) {
                std::cout << "File " << file_path << " was not found!\n";
                std::cout << "Current working directory is " << std::filesystem::current_path() << "!" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::ifstream file(file_path);
            std::string line;
            std::vector<uint32_t> n1, n2;

            while (std::getline(file, line)) {
                if (line[0] != '%') {
                    std::istringstream iss(line);
                    uint32_t a, b;
                    iss >> a >> b;

                    n1.push_back(a);
                    n2.push_back(b);
                }
            }
            file.close();

            uint32_t n = std::max(max(n1), max(n2)) + 1;
            size_t m = n1.size();
            reserve(n, m);

            for (size_t i = 0; i < m; ++i) {
                add_edge(n1[i], n2[i]);
            }
        };

        /**
         * Reserves space, so the graph can hold n nodes and m edges.
         *
         * @param n Number of nodes.
         * @param m Number of edges.
         */
        inline void reserve(size_t n, size_t m) {
            n_nodes = n;
            n_edges = m;
            adj_list.resize(n, std::vector<uint32_t>());
        }

        /**
         * Returns the number of nodes.
         *
         * @return Number of nodes.
         */
        inline size_t get_n() override {
            return n_nodes;
        };

        /**
         * Adds an edge to the graph.
         *
         * @param n1 Node 1.
         * @param n2 Node 2.
         */
        inline void add_edge(uint32_t n1, uint32_t n2) {
            adj_list[n1].push_back(n2);
            adj_list[n2].push_back(n1);
        };

        /**
         * Determines if an edge exists between th two vertices.
         *
         * @param n1 Vertex a.
         * @param n2 Vertex b.
         * @return True if there is an edge, False else.
         */
        inline bool edge_exists(uint32_t n1, uint32_t n2) {
            return contains(adj_list[n1], n2) && contains(adj_list[n2], n1);
        }

        /**
         * Determines if the graph consists of only one component.
         *
         * @return True if one component, False else.
         */
        inline bool is_one_component() {
            std::vector<uint8_t> processed(n_nodes, 0);
            std::vector<uint32_t> stack1(n_nodes);
            std::vector<uint32_t> stack2(n_nodes);

            stack1[0] = 0;
            size_t stack1_size = 1;
            size_t stack2_size = 0;

            while (stack1_size != 0) {
                for (size_t i = 0; i < stack1_size; ++i) {
                    uint32_t node = stack1[i];
                    processed[node] = 1;

                    for (auto n: adj_list[node]) {
                        if (!processed[n]) {
                            processed[n] = 1;
                            stack2[stack2_size++] = n;
                        }
                    }
                }

                stack1.swap(stack2);
                std::swap(stack1_size, stack2_size);
                stack2_size = 0;
            }

            return std::all_of(processed.begin(), processed.end(), [](uint8_t x) { return x == 1; });

        }

        /**
        * Sorts the neighborhood of each node by node-id.
        * Also removes unnecessary allocated memory.
        */
        inline void sort_unique_neighbours() {
            for (size_t i = 0; i < n_nodes; ++i) {
                std::sort(adj_list[i].begin(), adj_list[i].end());
                adj_list[i].erase(std::unique(adj_list[i].begin(), adj_list[i].end()), adj_list[i].end());
                adj_list[i].shrink_to_fit();
            }
        };

        /**
         * Returns a string representation of the graph.
         *
         * @return A string.
         */
        inline std::string get_as_string() const {
            std::stringstream output;
            output << "%" << n_nodes << " " << n_edges << "\n";

            for (uint32_t node = 0; node < n_nodes; ++node) {
                for (auto n: adj_list[node]) {
                    if (node < n) {
                        output << node << " " << n << "\n";
                    }
                }
            }
            return output.str();
        };

        /**
         * Writes the graph to a file.
         *
         * @param file_path The path to the file.
         */
        inline void write_to_file(std::string &file_path) const {
            std::stringstream output;
            output << "%" << n_nodes << " " << n_edges << "\n";

            for (uint32_t node = 0; node < n_nodes; ++node) {
                print(adj_list[node], adj_list[node].size());
                for (auto n: adj_list[node]) {
                    if (node < n) {
                        output << node << " " << n << "\n";
                    }
                }
            }

            std::ofstream file;
            file.open(file_path);
            file << output.rdbuf();
            file.close();
        };

        inline void finalize() override { sort_unique_neighbours(); };
    };

}

#endif //SUBSETOPTIMIZATION_GRAPH_H
