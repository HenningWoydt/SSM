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

#ifndef SUBSETOPTIMIZATION_STRUCTUREINTERFACE_H
#define SUBSETOPTIMIZATION_STRUCTUREINTERFACE_H

#include <vector>
#include <cstddef>
#include <cstdint>

namespace SMSM {

/**
 * Interface for structures, so that all score functions can be called in the
 * same way.
 *
 * @tparam TypeSF Datatype of the score function.
 */
    template<typename TypeSF>
    class StructureInterface {
    public:
        TypeSF max_reachable_score = std::numeric_limits<TypeSF>::max(); // maximum reachable score

        /**
         * Returns n.
         *
         * @return n.
         */
        virtual inline size_t get_n() = 0;

        /**
         * Gives the score for the empty set.
         *
         * @return Score of empty set.
         */
        virtual inline TypeSF evaluate_empty_set() = 0;

        /**
         * Gives the score for the set that has only added one element.
         *
         * @return Score of the set.
         */
        virtual inline TypeSF evaluate_1D(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
         * Gives the score for the set that has only added two elements.
         *
         * @return Score of the set.
         */
        virtual inline TypeSF evaluate_2D(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
         * Evaluates the score function, for sets that have added more than two
         * elements.
         *
         * @param s The set S.
         * @param s_size The size of set S.
         */
        virtual inline TypeSF evaluate_XD(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
         * Evaluates the score function in the most general way. The function
         * should be able to accept any set of any size.
         *
         * @param s The set S.
         * @param s_size The size of set S.
         */
        virtual inline TypeSF evaluate_general(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
        * Processes the graph after all edges have been added.
        */
        virtual inline void finalize() = 0;

        /**
         * Will initialize helping structures, that can help with score function
         * evaluation.
         *
         * @param k Number of depths, the search algorithm can explore.
         */
        virtual inline void initialize_helping_structures(size_t k) = 0;

        /**
         * Signals to the structure, that a new depth will be explored. Use this
         * function to keep the helping structures up to date.
         *
         * @param s The set S, that is present at the new depth.
         * @param s_size Size of the set S.
         */
        virtual inline void visit_new_depth(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
         * Signals to the structure, that we returned from the last depth. Use this
         * function to keep the helping structures up to date.
         */
        virtual inline void return_from_last_depth() = 0;
    };

}

#endif //SUBSETOPTIMIZATION_STRUCTUREINTERFACE_H
