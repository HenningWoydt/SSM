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

#include <vector>

#include <gtest/gtest.h>

#include "../src/structures/GraphPartialDominatingSet.h"
#include "../src/utility/GraphGeneration.h"
#include "../src/algorithms/BFSearch.h"
#include "../src/utility/AlgorithmConfiguration.h"
#include "../src/algorithms/TreeSearch.h"

namespace SMSM {

    void test_partial_dominating_set(size_t n_nodes, size_t n_graphs) {
        std::vector<GraphPartialDominatingSet<int>> graphs = random_GraphPartialDominatingSet_dataset(n_nodes, n_graphs);

        std::vector<AlgorithmConfiguration> acs = get_all_algorithm_configurations();

        auto bf_algorithm = BFSearch<GraphPartialDominatingSet<int>, int>(graphs[0], 1);

        for (auto &g: graphs) {
            for (uint32_t k = 1; k <= g.get_n(); ++k) {

                bf_algorithm.reinitialize(g, k);
                std::vector<uint32_t> bf_s = bf_algorithm.search();
                double bf_score = g.evaluate_general(bf_s, k);

                for (auto &ac: acs) {
                    ac.k = k;
                    ac.finalize();

                    g.initialize_helping_structures(k);
                    std::vector<uint32_t> s = TreeSearch<GraphPartialDominatingSet<int>, int>(g, ac).search();
                    double score = g.evaluate_general(s, k);

                    if (double_eq(bf_score, score, 0.0000001)) {
                        EXPECT_DOUBLE_EQ(bf_score, score);
                    } else {
                        EXPECT_DOUBLE_EQ(bf_score, score) <<
                                                          "\nScore Function: Partial-Dominating-Set" <<
                                                          "\nAlgorithm: " << ac.to_JSON() <<
                                                          "\nK: " << k <<
                                                          "\nBF Solution : " << to_string(bf_s) <<
                                                          "\nAlg Solution: " << to_string(s) <<
                                                          "\n" << g.get_as_string() << std::endl;
                    }
                }
            }
        }
    }

    TEST(GraphPartialDominatingSet, Manual_Hand1) {
        GraphPartialDominatingSet<int> g("../data/hand_drawn/hand1.edges");
        g.finalize();

        double score;
        auto s = std::vector<uint32_t>(5);

        score = g.evaluate_general(s, 0);
        EXPECT_DOUBLE_EQ(score, 0.0);

        s[0] = 0;
        score = g.evaluate_general(s, 1);
        EXPECT_DOUBLE_EQ(score, 2.0);

        s[0] = 1;
        score = g.evaluate_general(s, 1);
        EXPECT_DOUBLE_EQ(score, 4.0);

        s[0] = 2;
        score = g.evaluate_general(s, 1);
        EXPECT_DOUBLE_EQ(score, 3.0);

        s[0] = 3;
        score = g.evaluate_general(s, 1);
        EXPECT_DOUBLE_EQ(score, 3.0);

        s[0] = 4;
        score = g.evaluate_general(s, 1);
        EXPECT_DOUBLE_EQ(score, 3.0);

        s[0] = 0;
        s[1] = 1;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 4.0);

        s[0] = 0;
        s[1] = 2;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 4.0);

        s[0] = 0;
        s[1] = 3;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 4;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 4.0);

        s[0] = 1;
        s[1] = 2;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 1;
        s[1] = 3;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 1;
        s[1] = 4;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 2;
        s[1] = 3;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 4.0);

        s[0] = 2;
        s[1] = 4;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 4.0);

        s[0] = 3;
        s[1] = 4;
        score = g.evaluate_general(s, 2);
        EXPECT_DOUBLE_EQ(score, 4.0);

        s[0] = 0;
        s[1] = 1;
        s[2] = 2;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 1;
        s[2] = 3;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 1;
        s[2] = 4;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 2;
        s[2] = 3;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 2;
        s[2] = 4;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 3;
        s[2] = 4;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 1;
        s[1] = 2;
        s[2] = 3;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 1;
        s[1] = 2;
        s[2] = 4;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 1;
        s[1] = 3;
        s[2] = 4;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 2;
        s[1] = 3;
        s[2] = 4;
        score = g.evaluate_general(s, 3);
        EXPECT_DOUBLE_EQ(score, 4.0);

        s[0] = 0;
        s[1] = 1;
        s[2] = 2;
        s[3] = 3;
        score = g.evaluate_general(s, 4);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 1;
        s[2] = 2;
        s[3] = 4;
        score = g.evaluate_general(s, 4);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 1;
        s[2] = 3;
        s[3] = 4;
        score = g.evaluate_general(s, 4);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 2;
        s[2] = 3;
        s[3] = 4;
        score = g.evaluate_general(s, 4);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 1;
        s[1] = 2;
        s[2] = 3;
        s[3] = 4;
        score = g.evaluate_general(s, 4);
        EXPECT_DOUBLE_EQ(score, 5.0);

        s[0] = 0;
        s[1] = 1;
        s[2] = 2;
        s[3] = 3;
        s[4] = 4;
        score = g.evaluate_general(s, 5);
        EXPECT_DOUBLE_EQ(score, 5.0);
    }

    TEST(GraphPartialDominatingSet, Random4) {
        size_t n_nodes = 4;
        size_t n_graphs = 900;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random5) {
        size_t n_nodes = 5;
        size_t n_graphs = 850;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random6) {
        size_t n_nodes = 6;
        size_t n_graphs = 800;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random7) {
        size_t n_nodes = 7;
        size_t n_graphs = 750;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random8) {
        size_t n_nodes = 8;
        size_t n_graphs = 700;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random9) {
        size_t n_nodes = 9;
        size_t n_graphs = 650;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random10) {
        size_t n_nodes = 10;
        size_t n_graphs = 600;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random11) {
        size_t n_nodes = 11;
        size_t n_graphs = 550;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random12) {
        size_t n_nodes = 12;
        size_t n_graphs = 500;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random13) {
        size_t n_nodes = 13;
        size_t n_graphs = 450;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random14) {
        size_t n_nodes = 14;
        size_t n_graphs = 400;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random15) {
        size_t n_nodes = 15;
        size_t n_graphs = 350;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random16) {
        size_t n_nodes = 16;
        size_t n_graphs = 300;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random17) {
        size_t n_nodes = 17;
        size_t n_graphs = 250;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random18) {
        size_t n_nodes = 18;
        size_t n_graphs = 200;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random19) {
        size_t n_nodes = 19;
        size_t n_graphs = 150;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

    TEST(GraphPartialDominatingSet, Random20) {
        size_t n_nodes = 20;
        size_t n_graphs = 100;
        test_partial_dominating_set(n_nodes, n_graphs);
    }

}
