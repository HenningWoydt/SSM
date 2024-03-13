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

#include "JSONUtil.h"

namespace SMSM {

    std::string to_JSON_key(const bool b) {
        return "\"" + std::to_string(b) + "\"";
    }

    std::string to_JSON_key(const double d) {
        return "\"" + std::to_string(d) + "\"";
    }

    std::string to_JSON_key(const uint32_t n) {
        return "\"" + std::to_string(n) + "\"";
    }

    std::string to_JSON_key(const size_t n) {
        return "\"" + std::to_string(n) + "\"";
    }

    std::string to_JSON_key(int n) {
        return "\"" + std::to_string(n) + "\"";
    }

    std::string to_JSON_key(const std::string &s) {
        return "\"" + s + "\"";
    }

    std::string to_JSON_value(bool b) {
        return std::to_string(b);
    }

    std::string to_JSON_value(double d) {
        return std::to_string(d);
    }

    std::string to_JSON_value(uint32_t n) {
        return std::to_string(n);
    }

    std::string to_JSON_value(size_t n) {
        return std::to_string(n);
    }

    std::string to_JSON_value(int n) {
        return std::to_string(n);
    }

    std::string to_JSON_value(const std::string &s) {
        return "\"" + s + "\"";
    }

    std::vector<uint32_t> read_initial_vector(std::string &file_path) {
        std::vector<uint32_t> s;

        // open and read file
        std::string line;
        std::string content;
        std::ifstream file(file_path);
        if (file.is_open()) {
            while (getline(file, line)) {
                content += line;
            }
            file.close();
        } else {
            std::cout << "Could not read JSON file '" << file_path << "' as the initial solution file!" << std::endl;
            exit(EXIT_FAILURE);
        }

        // remove all '\n', ' ', '{', '}', '[', ']' characters
        content.erase(std::remove(content.begin(), content.end(), '\n'), content.end());
        content.erase(std::remove(content.begin(), content.end(), ' '), content.end());
        content.erase(std::remove(content.begin(), content.end(), '{'), content.end());
        content.erase(std::remove(content.begin(), content.end(), '}'), content.end());
        content.erase(std::remove(content.begin(), content.end(), '['), content.end());
        content.erase(std::remove(content.begin(), content.end(), ']'), content.end());

        // string should have the format '"s":x_1, x_2, ..., x_n'
        std::string s_vec = split(content, ':')[1];

        // get all numbers as strings
        std::vector<std::string> elements = split(s_vec, ',');

        // convert all strings to integers.
        for (const std::string &s_e: elements) {
            s.push_back(std::stoi(s_e));
        }
        std::sort(s.begin(), s.end());

        return s;
    }

}
