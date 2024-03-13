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

#include "DataPointGeneration.h"

namespace SMSM {

    DataPointsEuclidianDistance<double> random_DataPointsEuclidianDistance(size_t n, size_t d) {
        DataPointsEuclidianDistance<double> data_points(n, d);

        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_real_distribution<> dist(0.0, 1.0);

        dist(rng);

        std::vector<double> point(d, 0.0);

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < d; ++j) {
                point[j] = dist(rng);
            }
            data_points.add_entry(point);
        }
        data_points.determine_n_and_d();

        data_points.finalize();
        return data_points;
    }

    std::vector<DataPointsEuclidianDistance<double>> random_DataPointsEuclidianDistance_dataset(size_t n, size_t d, size_t n_sets) {
        std::vector<DataPointsEuclidianDistance<double>> vec;

        for (uint32_t i = 0; i < n_sets; ++i) {
            vec.emplace_back(random_DataPointsEuclidianDistance(n, d));
        }

        return vec;
    }

}
