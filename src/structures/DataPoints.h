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

#ifndef SUBSETOPTIMIZATION_DATAPOINTS_H
#define SUBSETOPTIMIZATION_DATAPOINTS_H

#include <cstddef>
#include <filesystem>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "../utility/Util.h"
#include "StructureInterface.h"

namespace SMSM {

/**
 * Class to hold data for kMeans.
 */
    template<typename TypeSF>
    class DataPoints : public StructureInterface<TypeSF> {
    public:
        size_t n_data_points = 0;
        size_t dimensionality = 0;
        std::vector<std::vector<TypeSF>> data_points;

        DataPoints() = default;

        /**
         * Constructs Datapoints.
         *
         * @param n Number of points.
         * @param d Dimensionality of each point.
         */
        DataPoints(size_t n, size_t d) {
            n_data_points = n;
            dimensionality = d;
        }

        /**
         * Reads a file with datapoints.
         *
         * @param file_path Path to the file.
         */
        explicit DataPoints(const std::string &file_path) {
            if (!file_exists(file_path)) {
                std::cout << "File " << file_path << " was not found!\n";
                std::cout << "Current working directory is " << std::filesystem::current_path() << "!" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::ifstream file(file_path);
            std::string line;

            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '%') {
                    continue;
                }

                std::vector<std::string> str_entries = split(line, ' ');
                std::vector<double> entries;
                for (auto &s: str_entries) {
                    entries.push_back(std::stod(s));
                }
                add_entry(entries);
            }
            file.close();
            determine_n_and_d();
        }

        /**
         * Adds an entry to the datapoints.
         *
         * @param vec Vector holding the entry.
         */
        void add_entry(const std::vector<TypeSF> &vec) {
            data_points.push_back(vec);
        }

        /**
         * Determine the number of datapoints and the dimensionality.
         */
        void determine_n_and_d() {
            n_data_points = data_points.size();
            dimensionality = data_points[0].size();
        }

        /**
         * Returns the number of datapoints.
         *
         * @return Number of data points.
         */
        size_t get_n() override {
            return n_data_points;
        }

        /**
         * Returns a string representation of the datapoints.
         *
         * @return The string.
         */
        std::string get_as_string() const {
            std::stringstream output;
            output << "%" << n_data_points << " " << dimensionality << "\n";

            for (size_t i = 0; i < n_data_points; ++i) {
                for (size_t d = 0; d < dimensionality; ++d) {
                    output << data_points[i][d] << " ";
                }
                output << "\n";
            }
            return output.str();
        }

        /**
         * Writes the datapoints to a file.
         *
         * @param file_path The path to the file.
         */
        inline void write_to_file(std::string &file_path) const {
            std::stringstream output;
            output << get_as_string() << "\n";
            std::ofstream file;
            file.open(file_path);
            file << output.rdbuf();
            file.close();
            std::cout << "write to " << file_path << std::endl;
        };

        /**
         * Checks that each datapoint has the correct dimensionality.
         */
        inline void finalize() override {
            for (size_t i = 0; i < n_data_points; ++i) {
                auto &vec = data_points[i];
                if (vec.size() != dimensionality) {
                    std::cout << "Datapoint " << i << " does not have dimensionality " << dimensionality << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        };
    };

}

#endif //SUBSETOPTIMIZATION_DATAPOINTS_H
