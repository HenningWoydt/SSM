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

#include "../src/structures/DataPointsEuclidianDistance.h"
#include "../src/utility/DataPointGeneration.h"
#include "../src/algorithms/BFSearch.h"
#include "../src/utility/AlgorithmConfiguration.h"
#include "../src/algorithms/TreeSearch.h"

namespace SMSM {

    void test_EuclidianDistance(size_t n, size_t d, size_t n_sets) {
        std::vector<DataPointsEuclidianDistance<double>> graphs = random_DataPointsEuclidianDistance_dataset(n, d, n_sets);

        std::vector<AlgorithmConfiguration> acs = get_all_algorithm_configurations();

        auto bf_algorithm = BFSearch<DataPointsEuclidianDistance<double>, double>(graphs[0], 1);

        for (auto &g: graphs) {
            for (uint32_t k = 1; k <= g.get_n(); ++k) {

                bf_algorithm.reinitialize(g, k);
                std::vector<uint32_t> bf_s = bf_algorithm.search();
                double bf_score = g.evaluate_general(bf_s, k);

                for (auto &ac: acs) {
                    ac.k = k;
                    ac.finalize();

                    g.initialize_helping_structures(k);
                    std::vector<uint32_t> s = TreeSearch<DataPointsEuclidianDistance<double>, double>(g, ac).search();
                    double score = g.evaluate_general(s, k);

                    if (double_eq(bf_score, score, 0.0000001)) {
                        EXPECT_DOUBLE_EQ(bf_score, score);
                    } else {
                        EXPECT_DOUBLE_EQ(bf_score, score) <<
                                                          "\nScore Function: Euclidian-Distance" <<
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

    TEST(DataPoints_EuclidianDistance, Random4) {
        size_t n = 4;
        size_t d = 5;
        size_t n_sets = 900;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random5) {
        size_t n = 5;
        size_t d = 5;
        size_t n_sets = 850;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random6) {
        size_t n = 6;
        size_t d = 5;
        size_t n_sets = 800;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random7) {
        size_t n = 7;
        size_t d = 5;
        size_t n_sets = 750;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random8) {
        size_t n = 8;
        size_t d = 5;
        size_t n_sets = 700;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random9) {
        size_t n = 9;
        size_t d = 5;
        size_t n_sets = 650;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random10) {
        size_t n = 10;
        size_t d = 5;
        size_t n_sets = 600;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random11) {
        size_t n = 11;
        size_t d = 5;
        size_t n_sets = 550;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random12) {
        size_t n = 12;
        size_t d = 5;
        size_t n_sets = 500;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random13) {
        size_t n = 13;
        size_t d = 5;
        size_t n_sets = 450;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random14) {
        size_t n = 14;
        size_t d = 5;
        size_t n_sets = 400;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random15) {
        size_t n = 15;
        size_t d = 5;
        size_t n_sets = 350;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random16) {
        size_t n = 16;
        size_t d = 5;
        size_t n_sets = 300;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random17) {
        size_t n = 17;
        size_t d = 5;
        size_t n_sets = 250;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random18) {
        size_t n = 18;
        size_t d = 5;
        size_t n_sets = 200;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random19) {
        size_t n = 19;
        size_t d = 5;
        size_t n_sets = 150;
        test_EuclidianDistance(n, d, n_sets);
    }

    TEST(DataPoints_EuclidianDistance, Random20) {
        size_t n = 20;
        size_t d = 5;
        size_t n_sets = 100;
        test_EuclidianDistance(n, d, n_sets);
    }

}
