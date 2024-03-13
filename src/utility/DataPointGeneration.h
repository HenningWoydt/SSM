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

#ifndef SUBSETOPTIMIZATION_DATAPOINTGENERATION_H
#define SUBSETOPTIMIZATION_DATAPOINTGENERATION_H

#include <cstddef>
#include <random>

#include "Util.h"
#include "../structures/DataPoints.h"
#include "../structures/DataPointsEuclidianDistance.h"

namespace SMSM {

    /**
     * Generates random datapoints.
     *
     * @param n Number of datapoints.
     * @param d Dimensionality of each datapoint.
     * @return The Datapoint object.
     */
    DataPointsEuclidianDistance<double> random_DataPointsEuclidianDistance(size_t n, size_t d);

    /**
     * Generates a random datapoints dataset.
     *
     * @param n Number of datapoints.
     * @param d Dimensionality of each datapoint.
     * @param n_sets Number of sets.
     * @return Vector holding the Datapoint objects.
     */
    std::vector<DataPointsEuclidianDistance<double>> random_DataPointsEuclidianDistance_dataset(size_t n, size_t d, size_t n_sets);

}

#endif //SUBSETOPTIMIZATION_DATAPOINTGENERATION_H
