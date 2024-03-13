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

#ifndef SUBSETOPTIMIZATION_DATAPOINTSEUCLIDIANDISTANCE_H
#define SUBSETOPTIMIZATION_DATAPOINTSEUCLIDIANDISTANCE_H

#include <cstddef>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>

#include "DataPoints.h"

namespace SMSM {

/**
 * Class to hold data for kMeans.
 */
    template<typename TypeSF>
    class DataPointsEuclidianDistance final : public DataPoints<TypeSF> {
    public:
        using DataPoints<TypeSF>::DataPoints;
        // distance matrix
        std::vector<std::vector<TypeSF>> dist_mtx;

        // structures to speed up score function evaluation
        size_t depth = 0;
        std::vector<std::vector<TypeSF>> min_dist;
        std::vector<TypeSF> temp_min;

        inline TypeSF evaluate_empty_set() override {
            TypeSF s = 0.0;
            for (size_t i = 0; i < DataPoints<TypeSF>::n_data_points; ++i) {
                s += sum(dist_mtx[i]);
            }
            return -s;
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
                return evaluate_empty_set();
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
            temp_min.resize(DataPoints<TypeSF>::n_data_points);
            initialize_dist_mtx();

            DataPoints<TypeSF>::max_reachable_score = 0;
        };

        inline void initialize_helping_structures(size_t k) override {
            min_dist.clear();
            min_dist.resize((k + 1), std::vector<TypeSF>(DataPoints<TypeSF>::n_data_points, std::numeric_limits<TypeSF>::max()));

            depth = 0;
        };

        inline void visit_new_depth(const std::vector<uint32_t> &s, size_t s_size) override {
            depth += 1;

            min(min_dist[depth], min_dist[depth - 1], dist_mtx[s[s_size - 1]]);
        };

        inline void return_from_last_depth() override {
            depth -= 1;
        };

        /**
        * Initializes the distance matrix.
        */
        inline void initialize_dist_mtx() {
            dist_mtx.resize(DataPoints<TypeSF>::n_data_points, std::vector<TypeSF>(DataPoints<TypeSF>::n_data_points, 0));
            for (size_t i = 0; i < DataPoints<TypeSF>::n_data_points; ++i) {
                for (size_t j = i + 1; j < DataPoints<TypeSF>::n_data_points; ++j) {

                    TypeSF distance = 0.0;
                    for (size_t d = 0; d < DataPoints<TypeSF>::dimensionality; ++d) {
                        distance += (DataPoints<TypeSF>::data_points[i][d] - DataPoints<TypeSF>::data_points[j][d]) * (DataPoints<TypeSF>::data_points[i][d] - DataPoints<TypeSF>::data_points[j][d]);
                    }
                    distance = sqrt(distance);

                    dist_mtx[i][j] = distance;
                    dist_mtx[j][i] = distance;
                }
            }
        };
    };

}

#endif //SUBSETOPTIMIZATION_DATAPOINTSEUCLIDIANDISTANCE_H
