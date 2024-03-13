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

#ifndef SUBSETOPTIMIZATION_VECTOROFVECTORS_H
#define SUBSETOPTIMIZATION_VECTOROFVECTORS_H

#include <vector>
#include <cstddef>

#include "Util.h"

namespace SMSM {

    template<class T>
    class VectorOfVectors {
    public:
        std::vector<std::vector<T>> vec;
        size_t size = 0;
        size_t capacity = 0;

        /**
         * Default constructor.
         */
        VectorOfVectors() = default;

        /**
         * Construct a VectorOfVectors with a given number of subvectors.
         *
         * @param n The initial number of subvectors.
         */
        explicit VectorOfVectors(size_t n) {
            vec.resize(n, std::vector<T>());
            size = 0;
            capacity = n;
        };

        /**
         * Copy subvectors from another VectorOfVectors.
         *
         * @param vec_of_vec The VectorOfVectors to copy from.
         */
        void copy_from(const VectorOfVectors<T> &vec_of_vec) {
            clear();

            for (size_t i = 0; i < vec_of_vec.size; ++i) {
                push_back(vec_of_vec.vec[i]);
            }
        };

        /**
         * Push a subvector to the end of the VectorOfVectors.
         *
         * @param v The subvector to push.
         */
        void push_back(const std::vector<T> &v) {
            if (size == capacity) {
                vec.push_back(v);
                size += 1;
                capacity += 1;
                return;
            }

            vec[size].resize(v.size());
            std::memcpy(vec[size].data(), v.data(), v.size() * sizeof(T));
            size += 1;
        };

        /**
         * Push a single element as a subvector to the end of the VectorOfVectors.
         *
         * @param t The element to push.
         */
        void push_back(const T &t) {
            if (size == capacity) {
                vec.push_back({t});
                size += 1;
                capacity += 1;
                return;
            }

            vec[size].resize(1);
            vec[size][0] = t;
            size += 1;
        }

        /**
         * Push two elements as a subvector to the end of the VectorOfVectors.
         *
         * @param t1 The first element.
         * @param t2 The second element.
         */
        void push_back(const T &t1, const T &t2) {
            if (size == capacity) {
                vec.push_back(std::vector<T>{t1, t2});
                size += 1;
                capacity += 1;
                return;
            }

            vec[size].resize(2);
            vec[size][0] = t1;
            vec[size][1] = t2;
            size += 1;
        }

        /**
         * Pushes the element into the last vector.
         *
         * @param t The element to push.
         */
        void push_back_to_last(const T t) {
            vec[size - 1].push_back(t);
        }

        /**
         * Reserve space for a specified number of subvectors.
         *
         * @param n The number of subvectors to reserve space for.
         */
        void reserve(size_t n) {
            vec.resize(n, std::vector<T>());
            capacity = n;
        };

        /**
         * Reserve space for a specified number of subvectors.
         *
         * @param n The number of subvectors to reserve space for.
         */
        void resize(size_t n) {
            vec.resize(n, std::vector<T>());
            capacity = n;
        };

        /**
         * Clear the VectorOfVectors by resetting the size.
         */
        void clear() {
            size = 0;
        };

        /**
         * Check if a specified element exists in any subvector.
         *
         * @param x The element to search for.
         * @return True if the element is found; false otherwise.
         */
        bool contains(const T &x) {
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = 0; j < vec[i].size(); ++j) {
                    if (vec[i][j] == x) {
                        return true;
                    }
                }
            }
            return false;
        }

        /**
         * Check if a specified element exists in any subvector.
         *
         * @param x The element to search for.
         * @return True if the element is found; false otherwise.
         */
        bool contains(const T &&x) {
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = 0; j < vec[i].size(); ++j) {
                    if (vec[i][j] == x) {
                        return true;
                    }
                }
            }
            return false;
        }
    };

}

#endif //SUBSETOPTIMIZATION_VECTOROFVECTORS_H
