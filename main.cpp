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

#include "src/algorithms/TreeSearch.h"
#include "src/utility/AlgorithmConfiguration.h"
#include "src/utility/GraphGeneration.h"
#include "src/utility/DataPointGeneration.h"
#include "src/structures/FacilityLocation.h"

int main(int argc, char *argv[]) {
    SMSM::AlgorithmConfiguration ac = SMSM::parse_command_line(argc, argv, true);
    if (ac.invalid) {
        std::cout << ac.to_JSON() << std::endl;
        return -1;
    }

    if (ac.structure_type == "graph") {
        if (ac.score_function == "negative-group-farness") {
            SMSM::GraphNegativeGroupFarness<int> g(ac.input_file_path);
            g.finalize();
            g.initialize_helping_structures(ac.k);
            if (g.get_n() < ac.k) {
                std::cout << "n (" << g.get_n() << ") is smaller than k (" << ac.k << ")!" << std::endl;
                return -1;
            }
            auto ts = SMSM::TreeSearch<SMSM::GraphNegativeGroupFarness<int>, int>(g, ac);
            ts.search();
        } else if (ac.score_function == "partial-dominating-set") {
            SMSM::GraphPartialDominatingSet<int> g(ac.input_file_path);
            g.finalize();
            g.initialize_helping_structures(ac.k);
            if (g.get_n() < ac.k) {
                std::cout << "n (" << g.get_n() << ") is smaller than k (" << ac.k << ")!" << std::endl;
                return -1;
            }

            auto ts = SMSM::TreeSearch<SMSM::GraphPartialDominatingSet<int>, int>(g, ac);
            ts.search();
        }
    } else if (ac.structure_type == "k-medoid") {
        if (ac.score_function == "euclidian-distance") {
            SMSM::DataPointsEuclidianDistance<double> dp(ac.input_file_path);
            dp.finalize();
            dp.initialize_helping_structures(ac.k);
            if (dp.get_n() < ac.k) {
                std::cout << "n (" << dp.get_n() << ") is smaller than k (" << ac.k << ")!" << std::endl;
                return -1;
            }

            auto ts = SMSM::TreeSearch<SMSM::DataPointsEuclidianDistance<double>, double>(dp, ac);
            ts.search();
        } else {
            std::cout << "Score function '" << ac.score_function << "' not known for structure type '" << ac.structure_type << "'!" << std::endl;
            return -1;
        }
    } else if (ac.structure_type == "facility") {
        if (ac.score_function == "benefits") {
            SMSM::FacilityLocation<double> fl(ac.input_file_path);
            fl.finalize();
            fl.initialize_helping_structures(ac.k);
            if (fl.get_n() < ac.k) {
                std::cout << "n (" << fl.get_n() << ") is smaller than k (" << ac.k << ")!" << std::endl;
                return -1;
            }

            auto ts = SMSM::TreeSearch<SMSM::FacilityLocation<double>, double>(fl, ac);
            ts.search();
        } else {
            std::cout << "Score function '" << ac.score_function << "' not known for structure type '" << ac.structure_type << "'!" << std::endl;
            return -1;
        }
    } else {
        std::cout << "Structure type '" << ac.structure_type << "' not known!" << std::endl;
        return -1;
    }
}
