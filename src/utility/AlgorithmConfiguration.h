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

#ifndef SUBSETOPTIMIZATION_ALGORITHMCONFIGURATION_H
#define SUBSETOPTIMIZATION_ALGORITHMCONFIGURATION_H


#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <limits>

#include <boost/program_options.hpp>

#include "Util.h"
#include "JSONUtil.h"

namespace SMSM {
    /**
     * Class to configurate the algorithm.
     */
    class AlgorithmConfiguration {
    public:
        bool invalid = false; // if the ac is valid or not

        std::string structure_type; // 'graph', 'k-medoid'
        std::string input_file_path; // path to a file holding the data
        size_t k; // required size of the set
        std::string score_function; // string specifying the score function

        std::vector<uint32_t> partial_s; // partial solution to S, the algorithm will only fill in the last k' elements (optional)

        std::vector<uint32_t> init_candidates; // Initial candidate list (optional)
        std::vector<double> init_si; // Initial score improvements (optional)
        std::vector<uint8_t> init_acc; // Initial accurate state of the score improvements (optional)

        double score_threshold = std::numeric_limits<double>::max(); // if a set of size k with a greater score has been found, it will return

        bool write_output = false; // whether to write output to a file
        std::string output_file_path = "../out.JSON"; // path specifying the file

        bool measure_oracle_time = false; // whether to measure oracle time

        double time_limit = std::numeric_limits<double>::max(); // given time limit

        // Configuration of brute force
        size_t bf_threshold_n = 1;
        size_t bf_threshold_r = 1;

        bool greedy_local_search = false;

        // if to use the plain algorithm
        bool plain = false;

        // Configuration for SUB
        bool SUB_enabled = false;

        // Configuration for RPC
        bool RPC_enabled = false;

        // Configuration for Lazy Evaluation
        size_t LE_mode = 0; // 0 - disabled, 1 - score or rank, 2 - score and rank
        double LE_y_rank_value = 0; // the rank y value
        size_t LE_rank_var = 0; // 0 - none, 1 - n, 2 - k, 3 - r_n, 4 - r_k
        double LE_y_score_value = 0; // the score y value

        // Configuration for the Upper Bound 2D heuristic
        bool UB2D_enabled = false;
        size_t UB2D_l_func = 0; // 0 - invalid, 1 - sqrt()*y, 2 - id()*y, 3 - y
        size_t UB2D_l_var = 0; // 0 - no, 1 - n, 2 - k, 3 - r_n, 4 - r_k
        double UB2D_l_y = 0; // y value
        size_t UB2D_max_l = 0; // cap on how large the hyperparameter l can be
        uint32_t UB2D_alg_type = 0; // 0 - invalid, 1 - Greedy, 2 - Matching, 3 - Brute Force, 4 - Dynamic
        uint32_t UB2D_odd_type = 0; // 0 - invalid, 1 - A (1 + pairs), 2 - B (pairs + 1), 3 - min(A, B)
        bool UB2D_RPC_enabled = true; // If we use RPC to remove elements
        bool UB2D_safe_skip_enabled = false; // If we try to identify if we have tot start UB2D
        double UB2D_lazy_skip_start_value = 0.0; // start value for lazy skip heuristic
        double UB2D_lazy_skip_add_value = 1.0; // value to add to lazy skips
        size_t UB2D_low_depth = 0; // lowest depth to start UB2D
        size_t UB2D_high_depth = std::numeric_limits<size_t>::max(); // greatest depth to start UB2D
        std::vector<bool> UB2D_depth_enabled; // vector holding information whether to start UB2D on a given depth - automatically filled
        std::vector<bool> UB2D_manual_depth_enabled; // vector holding information to start UB2D on a given depth - manually filled (will be prioritized over the automatic vector)
        double UB2D_sub_bound_percentage = 0.0; // percentage specifying how close the sub bound must be to the best score, so that UB2D gets started
        std::vector<double> UB2D_sub_bound_percentage_vec; // vector holding a percentage for each depth - automatically filled
        std::vector<double> UB2D_manual_sub_bound_percentage_vec; // vector holding a percentage for each depth - manually filled (will be prioritized over the automatic vector)

        // Configuration for the Partial Brute Force heuristic
        bool PBF_enabled = false;
        size_t PBF_n_func = 0; // 0 - invalid, 1 - sqrt()*y, 2 - id()*y, 3 - y
        size_t PBF_n_var = 0;  // 0 - no, 1 - n, 2 - k, 3 - r_n, 4 - r_k
        double PBF_n_y = 0; // y value
        size_t PBF_max_n = 0; // cap on how large n can be
        size_t PBF_l_func = 0;  // 0 - invalid, 1 - sqrt()*y, 2 - id()*y, 3 - y
        size_t PBF_l_var = 0;  // 0 - no, 1 - n, 2 - k, 3 - r_n, 4 - r_k
        double PBF_l_y = 0; // y value
        size_t PBF_max_l = 0; // cap on how large l can be
        size_t PBF_alg_type = 0; // 0 - invalid, 1 - bforce, 2 - dynamic
        bool PBF_RPC_enabled = true; // if we use RPC to remove elements
        bool PBF_safe_skip_enabled = false; // // If we try to identify if we have to start PBF
        double PBF_lazy_skip_start_value = 0.0; // start value for lazy skip heuristic
        double PBF_lazy_skip_add_value = 1.0; // value to add to lazy skips
        size_t PBF_low_depth = 0; // lowest depth to start PBF
        size_t PBF_high_depth = std::numeric_limits<size_t>::max(); // greatest depth to start PBF
        std::vector<bool> PBF_depth_enabled; // vector holding information whether to start PBF on a given depth - automatically filled
        std::vector<bool> PBF_manual_depth_enabled; // vector holding information to start PBF on a given depth - manually filled (will be prioritized over the automatic vector)
        double PBF_sub_bound_percentage = 0.0; // percentage specifying how close the sub bound must be to the best score, so that PBF gets started
        std::vector<double> PBF_sub_bound_percentage_vec; // vector holding a percentage for each depth - automatically filled
        std::vector<double> PBF_manual_sub_bound_percentage_vec; // vector holding a percentage for each depth - manually filled (will be prioritized over the automatic vector)

        // Configuration for the Recursive heuristic
        bool REC_enabled = false;
        size_t REC_l_func = 0; // 0 - invalid, 1 - sqrt()*y, 2 - id()*y, 3 - y
        size_t REC_l_var = 0; // 0 - no, 1 - n, 2 - k, 3 - r_n, 4 - r_k
        double REC_l_y = 0; // y value
        size_t REC_max_l = 0; // cap on l
        bool REC_RPC_enabled = true; // if we use RPC to remove elements
        bool REC_safe_skip_enabled = false; // // If we try to identify if we have to start REC
        double REC_lazy_skip_start_value = 0.0; // start value for lazy skip heuristic
        double REC_lazy_skip_add_value = 1.0; // value to add to lazy skips
        size_t REC_low_depth = 0; // lowest depth to start REC
        size_t REC_high_depth = std::numeric_limits<size_t>::max(); // greatest depth to start REC
        std::vector<bool> REC_depth_enabled; // vector holding information whether to start REC on a given depth - automatically filled
        std::vector<bool> REC_manual_depth_enabled; // vector holding information to start REC on a given depth - manually filled (will be prioritized over the automatic vector)
        double REC_sub_bound_percentage = 0.0; // percentage specifying how close the sub bound must be to the best score, so that REC gets started
        std::vector<double> REC_sub_bound_percentage_vec; // vector holding a percentage for each depth - automatically filled
        std::vector<double> REC_manual_sub_bound_percentage_vec; // vector holding a percentage for each depth - manually filled (will be prioritized over the automatic vector)

        /**
         * Finalizes the configuartion.
         */
        void finalize() {
            determine_UB2D_depth_active();
            determine_UB2D_sub_bound_percentage_vec();
            determine_PBF_depth_active();
            determine_PBF_sub_bound_percentage_vec();
            determine_REC_depth_active();
            determine_REC_sub_bound_percentage_vec();
        }

        /**
         * Determines the Lazy Evaluation rank threshold.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @param t_r_n Number of remaining candidates.
         * @param t_r_k Remaining size of the solution.
         * @return The threshold.
         */
        size_t determine_LE_rank_threshold(size_t t_n, size_t t_k, size_t t_r_n, size_t t_r_k) const {
            size_t var_value = 1;
            if (LE_rank_var == 0) {
                var_value = 1;
            } else if (LE_rank_var == 1) {
                var_value = t_n;
            } else if (LE_rank_var == 2) {
                var_value = t_k;
            } else if (LE_rank_var == 3) {
                var_value = t_r_n;
            } else if (LE_rank_var == 4) {
                var_value = t_r_k;
            }

            auto threshold = (size_t) (((double) var_value) * LE_y_rank_value);

            return std::min(threshold, t_n);
        }

        /**
         * Determines the Lazy Evaluation score threshold.
         *
         * @param r_score_avg Remaining score per element.
         * @return The threshold.
         */
        double determine_LE_score_threshold(double r_score_avg) const {
            return r_score_avg * LE_y_score_value;
        }

        /**
         * Determines the maximum value for 'l' for Upper Bound 2D.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @return The maximum 'l'.
         */
        size_t determine_UB2D_max_l(size_t t_n, size_t t_k) const {
            size_t var_max = t_n;
            if (UB2D_l_var == 1 || UB2D_l_var == 3) {
                var_max = t_n;
            } else if (UB2D_l_var == 2 || UB2D_l_var == 4) {
                var_max = t_k;
            }

            auto res = (size_t) (sqrt((double) var_max) * UB2D_l_y);
            if (UB2D_l_func == 1) {
                res = (size_t) (sqrt((double) var_max) * UB2D_l_y);
            } else if (UB2D_l_func == 2) {
                res = (size_t) (((double) var_max) * UB2D_l_y);
            } else if (UB2D_l_func == 3) {
                res = (size_t) UB2D_l_y;
            }

            res = std::min(res, UB2D_max_l);
            return res;
        }

        /**
         * Determines the 'l' value for Upper Bound 2D.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @param t_r_n Number of remaining candidates.
         * @param t_r_k Remaining size of the solution.
         * @return The 'l' value.
         */
        size_t determine_UB2D_l(size_t t_n, size_t t_k, size_t t_r_n, size_t t_r_k) const {
            size_t var = t_n;
            if (UB2D_l_var == 1) {
                var = t_n;
            } else if (UB2D_l_var == 2) {
                var = t_k;
            } else if (UB2D_l_var == 3) {
                var = t_r_n;
            } else if (UB2D_l_var == 4) {
                var = t_r_k;
            }

            auto res = (size_t) (sqrt((double) var) * UB2D_l_y);
            if (UB2D_l_func == 1) {
                res = (size_t) (sqrt((double) var) * UB2D_l_y);
            } else if (UB2D_l_func == 2) {
                res = (size_t) (((double) var) * UB2D_l_y);
            } else if (UB2D_l_func == 3) {
                res = (size_t) UB2D_l_y;
            }

            res = std::min(res, UB2D_max_l);
            return res;
        }

        /**
         * Determines the maximum value for 'n' for Partial Brute Force.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @return The maximum 'n' value.
         */
        size_t determine_PBF_max_n(size_t t_n, size_t t_k) const {
            size_t var_max = t_n;
            if (PBF_n_var == 1 || PBF_n_var == 3) {
                var_max = t_n;
            } else if (PBF_n_var == 2 || PBF_n_var == 4) {
                var_max = t_k;
            }

            auto res = (size_t) (sqrt((double) var_max) * PBF_n_y);
            if (PBF_n_func == 1) {
                res = (size_t) (sqrt((double) var_max) * PBF_n_y);
            } else if (PBF_n_func == 2) {
                res = (size_t) (((double) var_max) * PBF_n_y);
            } else if (PBF_n_func == 3) {
                res = (size_t) PBF_n_y;
            }

            res = std::min(res, PBF_max_n);
            return res;
        }

        /**
         * Determines the maximum value for 'l' for Partial Brute Force.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @return The maximum 'l' value.
         */
        size_t determine_PBF_max_l(size_t t_n, size_t t_k) const {
            size_t var_max = t_n;
            if (PBF_l_var == 1 || PBF_l_var == 3) {
                var_max = t_n;
            } else if (PBF_l_var == 2 || PBF_l_var == 4) {
                var_max = t_k;
            }

            auto res = (size_t) (sqrt((double) var_max) * PBF_l_y);
            if (PBF_l_func == 1) {
                res = (size_t) (sqrt((double) var_max) * PBF_l_y);
            } else if (PBF_l_func == 2) {
                res = (size_t) (((double) var_max) * PBF_l_y);
            } else if (PBF_l_func == 3) {
                res = (size_t) PBF_l_y;
            }

            res = std::min(res, PBF_max_n);
            return res;
        }

        /**
         * Determines the 'n' value for Partial Brute Force.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @param t_r_n Number of remaining candidates.
         * @param t_r_k Remaining size of the solution.
         * @return The 'n' value.
         */
        size_t determine_PBF_n(size_t t_n, size_t t_k, size_t t_r_n, size_t t_r_k) const {
            size_t var = t_n;
            if (PBF_n_var == 1) {
                var = t_n;
            } else if (PBF_n_var == 2) {
                var = t_k;
            } else if (PBF_n_var == 3) {
                var = t_r_n;
            } else if (PBF_n_var == 4) {
                var = t_r_k;
            }

            auto res = (size_t) (sqrt((double) var) * PBF_n_y);
            if (PBF_n_func == 1) {
                res = (size_t) (sqrt((double) var) * PBF_n_y);
            } else if (PBF_n_func == 2) {
                res = (size_t) (((double) var) * PBF_n_y);
            } else if (PBF_n_func == 3) {
                res = (size_t) PBF_n_y;
            }

            res = std::min(res, PBF_max_n);
            return res;
        }

        /**
         * Determines the 'l' value for Partial Brute Force.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @param t_r_n Number of remaining candidates.
         * @param t_r_k Remaining size of the solution.
         * @return The 'l' value.
         */
        size_t determine_PBF_l(size_t t_n, size_t t_k, size_t t_r_n, size_t t_r_k) const {
            size_t var = t_n;
            if (PBF_l_var == 1) {
                var = t_n;
            } else if (PBF_l_var == 2) {
                var = t_k;
            } else if (PBF_l_var == 3) {
                var = t_r_n;
            } else if (PBF_l_var == 4) {
                var = t_r_k;
            }

            auto res = (size_t) (sqrt((double) var) * PBF_l_y);
            if (PBF_l_func == 1) {
                res = (size_t) (sqrt((double) var) * PBF_l_y);
            } else if (PBF_l_func == 2) {
                res = (size_t) (((double) var) * PBF_l_y);
            } else if (PBF_l_func == 3) {
                res = (size_t) PBF_l_y;
            }

            res = std::min(res, PBF_max_l);
            return res;
        }

        /**
         * Determines the maximum value for 'l' for Recursive.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @return The maximum 'l' value.
         */
        size_t determine_REC_max_l(size_t t_n, size_t t_k) const {
            size_t var_max = t_n;
            if (REC_l_var == 1 || REC_l_var == 3) {
                var_max = t_n;
            } else if (REC_l_var == 2 || REC_l_var == 4) {
                var_max = t_k;
            }

            auto res = (size_t) (sqrt((double) var_max) * REC_l_y);
            if (REC_l_func == 1) {
                res = (size_t) (sqrt((double) var_max) * REC_l_y);
            } else if (REC_l_func == 2) {
                res = (size_t) (((double) var_max) * REC_l_y);
            } else if (REC_l_func == 3) {
                res = (size_t) REC_l_y;
            }

            res = std::min(res, REC_max_l);
            return res;
        }

        /**
         * Determines the 'l' value for Recursive.
         *
         * @param t_n Number of candidates.
         * @param t_k Size of the solution.
         * @param t_r_n Number of remaining candidates.
         * @param t_r_k Remaining size of the solution.
         * @return The 'l' value.
         */
        size_t determine_REC_l(size_t t_n, size_t t_k, size_t t_r_n, size_t t_r_k) const {
            size_t var = t_n;
            if (REC_l_var == 1) {
                var = t_n;
            } else if (REC_l_var == 2) {
                var = t_k;
            } else if (REC_l_var == 3) {
                var = t_r_n;
            } else if (REC_l_var == 4) {
                var = t_r_k;
            }

            auto res = (size_t) (sqrt((double) var) * REC_l_y);
            if (REC_l_func == 1) {
                res = (size_t) (sqrt((double) var) * REC_l_y);
            } else if (REC_l_func == 2) {
                res = (size_t) (REC_l_y * (double) var);
            } else if (REC_l_func == 3) {
                res = (size_t) REC_l_y;
            }

            return res;
        }

        /**
         * Determines on which depth Upper Bound 2D is active.
         */
        void determine_UB2D_depth_active() {
            UB2D_depth_enabled.resize(k, false);

            if (UB2D_manual_depth_enabled.size() == k) {
                for (size_t i = 0; i < k; ++i) {
                    UB2D_depth_enabled[i] = UB2D_manual_depth_enabled[i];
                }
            } else {
                for (size_t i = 0; i < k; ++i) {
                    UB2D_depth_enabled[i] = false;
                    if (UB2D_low_depth <= i && i <= UB2D_high_depth) {
                        UB2D_depth_enabled[i] = true;
                    }
                }
            }
        }

        /**
         * Determines the Simple Upper Bound percentage for each depth for Upper
         * Bound 2D.
         */
        void determine_UB2D_sub_bound_percentage_vec() {
            UB2D_sub_bound_percentage_vec.resize(k, UB2D_sub_bound_percentage);

            if (UB2D_manual_sub_bound_percentage_vec.size() == k) {
                for (size_t i = 0; i < k; ++i) {
                    UB2D_sub_bound_percentage_vec[i] = UB2D_manual_sub_bound_percentage_vec[i];
                }
            } else {
                for (size_t i = 0; i < k; ++i) {
                    UB2D_sub_bound_percentage_vec[i] = UB2D_sub_bound_percentage;
                }
            }
        }

        /**
         * Determines on which depth Partial Brute Force is active.
         */
        void determine_PBF_depth_active() {
            PBF_depth_enabled.resize(k, false);

            if (PBF_manual_depth_enabled.size() == k) {
                for (size_t i = 0; i < k; ++i) {
                    PBF_depth_enabled[i] = PBF_manual_depth_enabled[i];
                }
            } else {
                for (size_t i = 0; i < k; ++i) {
                    PBF_depth_enabled[i] = false;
                    if (PBF_low_depth <= i && i <= PBF_high_depth) {
                        PBF_depth_enabled[i] = true;
                    }
                }
            }
        }

        /**
         * Determines the Simple Upper Bound percentage for each depth for
         * Partial Brute Force.
         */
        void determine_PBF_sub_bound_percentage_vec() {
            PBF_sub_bound_percentage_vec.resize(k, PBF_sub_bound_percentage);

            if (PBF_manual_sub_bound_percentage_vec.size() == k) {
                for (size_t i = 0; i < k; ++i) {
                    PBF_sub_bound_percentage_vec[i] = PBF_manual_sub_bound_percentage_vec[i];
                }
            } else {
                for (size_t i = 0; i < k; ++i) {
                    PBF_sub_bound_percentage_vec[i] = PBF_sub_bound_percentage;
                }
            }
        }

        /**
         * Determines on which depth Recursive is active.
         */
        void determine_REC_depth_active() {
            REC_depth_enabled.resize(k, false);
            if (REC_manual_depth_enabled.size() == k) {
                for (size_t i = 0; i < k; ++i) {
                    REC_depth_enabled[i] = REC_manual_depth_enabled[i];
                }
            } else {
                for (size_t i = 0; i < k; ++i) {
                    REC_depth_enabled[i] = false;
                    if (REC_low_depth <= i && i <= REC_high_depth) {
                        REC_depth_enabled[i] = true;
                    }
                }
            }
        }

        /**
         * Determines the Simple Upper Bound percentage for each depth for
         * Recursive.
         */
        void determine_REC_sub_bound_percentage_vec() {
            REC_sub_bound_percentage_vec.resize(k, REC_sub_bound_percentage);

            if (REC_manual_sub_bound_percentage_vec.size() == k) {
                for (size_t i = 0; i < k; ++i) {
                    REC_sub_bound_percentage_vec[i] = REC_manual_sub_bound_percentage_vec[i];
                }
            } else {
                for (size_t i = 0; i < k; ++i) {
                    REC_sub_bound_percentage_vec[i] = REC_sub_bound_percentage;
                }
            }
        }

        /**
         * Writes each variable in JSON format to one string.
         *
         * @return The string.
         */
        std::string to_JSON() {
            std::string content;

            content += "\"structure-type\" : " + to_JSON_value(structure_type) + ",\n";
            content += "\"input-file-path\" : " + to_JSON_value(input_file_path) + ",\n";
            content += "\"k\" : " + to_JSON_value(k) + ",\n";
            content += "\"score-function\" : " + to_JSON_value(score_function) + ",\n";

            content += "\"partial-solution\" : " + SMSM::to_JSON(partial_s) + ",\n";
            content += "\"initial-candidates\" : " + SMSM::to_JSON(init_candidates) + ",\n";
            content += "\"initial-score-improvements\" : " + SMSM::to_JSON(init_si) + ",\n";
            content += "\"initial-accurate\" : " + SMSM::to_JSON(init_acc) + ",\n";

            content += "\"score-threshold\" : " + SMSM::to_JSON_value(score_threshold) + ",\n";
            content += "\"write-output\" : " + to_JSON_value(write_output) + ",\n";
            content += "\"output-file-path\" : " + to_JSON_value(output_file_path) + ",\n";
            content += "\"measure-oracle-time\" : " + to_JSON_value(measure_oracle_time) + ",\n";
            content += "\"time-limit\" : " + to_JSON_value(time_limit) + ",\n";
            content += "\"bf-threshold-n\" : " + to_JSON_value(bf_threshold_n) + ",\n";
            content += "\"bf-threshold-r\" : " + to_JSON_value(bf_threshold_r) + ",\n";
            content += "\"SUB-enabled\" : " + to_JSON_value(SUB_enabled) + ",\n";
            content += "\"RPC-enabled\" : " + to_JSON_value(RPC_enabled) + ",\n";

            content += "\"LE-mode\" : " + to_JSON_value(LE_mode) + ",\n";
            content += "\"LE-y-rank-value\" : " + to_JSON_value(LE_y_rank_value) + ",\n";
            content += "\"LE-rank-var\" : " + to_JSON_value(LE_rank_var) + ",\n";
            content += "\"LE-y-score-value\" : " + to_JSON_value(LE_y_score_value) + ",\n";

            content += "\"UB2D-enabled\" : " + to_JSON_value(UB2D_enabled) + ",\n";
            content += "\"UB2D-l-func\" : " + to_JSON_value(UB2D_l_func) + ",\n";
            content += "\"UB2D-l-var\" : " + to_JSON_value(UB2D_l_var) + ",\n";
            content += "\"UB2D-l-y\" : " + to_JSON_value(UB2D_l_y) + ",\n";
            content += "\"UB2D-max-l\" : " + to_JSON_value(UB2D_max_l) + ",\n";
            content += "\"UB2D-alg-type\" : " + to_JSON_value(UB2D_alg_type) + ",\n";
            content += "\"UB2D-odd-type\" : " + to_JSON_value(UB2D_odd_type) + ",\n";
            content += "\"UB2D-RPC-enabled\" : " + to_JSON_value(UB2D_RPC_enabled) + ",\n";
            content += "\"UB2D-safe-skips-enabled\" : " + to_JSON_value(UB2D_safe_skip_enabled) + ",\n";
            content += "\"UB2D-lazy-skip-start-value\" : " + to_JSON_value(UB2D_lazy_skip_start_value) + ",\n";
            content += "\"UB2D-lazy-skip-skip-add-value\" : " + to_JSON_value(UB2D_lazy_skip_add_value) + ",\n";
            content += "\"UB2D-low-depth\" : " + to_JSON_value(UB2D_low_depth) + ",\n";
            content += "\"UB2D-high-depth\" : " + to_JSON_value(UB2D_high_depth) + ",\n";
            content += "\"UB2D-depth-enabled\" : " + to_JSON_value(UB2D_depth_enabled) + ",\n";
            content += "\"UB2D-manual-depth-enabled\" : " + to_JSON_value(UB2D_manual_depth_enabled) + ",\n";
            content += "\"UB2D-sub-bound-percentage\" : " + to_JSON_value(UB2D_sub_bound_percentage) + ",\n";
            content += "\"UB2D-sub-bound-percentage-vec\" : " + to_JSON_value(UB2D_sub_bound_percentage_vec) + ",\n";
            content += "\"UB2D-manual-sub-bound-percentage-vec\" : " + to_JSON_value(UB2D_manual_sub_bound_percentage_vec) + ",\n";

            content += "\"PBF-enabled\" : " + to_JSON_value(PBF_enabled) + ",\n";
            content += "\"PBF-n-func\" : " + to_JSON_value(PBF_n_func) + ",\n";
            content += "\"PBF-n-var\" : " + to_JSON_value(PBF_n_var) + ",\n";
            content += "\"PBF-n-y\" : " + to_JSON_value(PBF_n_y) + ",\n";
            content += "\"PBF-max-n\" : " + to_JSON_value(PBF_max_n) + ",\n";
            content += "\"PBF-l-func\" : " + to_JSON_value(PBF_l_func) + ",\n";
            content += "\"PBF-l-var\" : " + to_JSON_value(PBF_l_var) + ",\n";
            content += "\"PBF-l-y\" : " + to_JSON_value(PBF_l_y) + ",\n";
            content += "\"PBF-max-l\" : " + to_JSON_value(PBF_max_l) + ",\n";
            content += "\"PBF-algorithm-type\" : " + to_JSON_value(PBF_alg_type) + ",\n";
            content += "\"PBF-RPC-enabled\" : " + to_JSON_value(PBF_RPC_enabled) + ",\n";
            content += "\"PBF-safe-skip-enabled\" : " + to_JSON_value(PBF_safe_skip_enabled) + ",\n";
            content += "\"PBF-lazy-skip-start-value\" : " + to_JSON_value(PBF_lazy_skip_start_value) + ",\n";
            content += "\"PBF-lazy-skip-add-value\" : " + to_JSON_value(PBF_lazy_skip_add_value) + ",\n";
            content += "\"PBF-low-depth\" : " + to_JSON_value(PBF_low_depth) + ",\n";
            content += "\"PBF-high-depth\" : " + to_JSON_value(PBF_high_depth) + ",\n";
            content += "\"PBF-depth-enabled\" : " + to_JSON_value(PBF_depth_enabled) + ",\n";
            content += "\"PBF-manual-depth-enabled\" : " + to_JSON_value(PBF_manual_depth_enabled) + ",\n";
            content += "\"PBF-sub-bound-percentage\" : " + to_JSON_value(PBF_sub_bound_percentage) + ",\n";
            content += "\"PBF-sub-bound-percentage-vec\" : " + to_JSON_value(PBF_sub_bound_percentage_vec) + ",\n";
            content += "\"PBF-manual-sub-bound-percentage-vec\" : " + to_JSON_value(PBF_sub_bound_percentage_vec) + ",\n";

            content += "\"REC-enabled\" : " + to_JSON_value(REC_enabled) + ",\n";
            content += "\"REC-l-func\" : " + to_JSON_value(REC_l_func) + ",\n";
            content += "\"REC-l-var\" : " + to_JSON_value(REC_l_var) + ",\n";
            content += "\"REC-l-y\" : " + to_JSON_value(REC_l_y) + ",\n";
            content += "\"REC-max-l\" : " + to_JSON_value(REC_max_l) + ",\n";
            content += "\"REC-RPC-enabled\" : " + to_JSON_value(REC_RPC_enabled) + ",\n";
            content += "\"REC-safe-skip-enabled\" : " + to_JSON_value(REC_safe_skip_enabled) + ",\n";
            content += "\"REC-lazy-skip-start-value\" : " + to_JSON_value(REC_lazy_skip_start_value) + ",\n";
            content += "\"REC-lazy-skip-add-value\" : " + to_JSON_value(REC_lazy_skip_add_value) + ",\n";
            content += "\"REC-low-depth\" : " + to_JSON_value(REC_low_depth) + ",\n";
            content += "\"REC-high-depth\" : " + to_JSON_value(REC_high_depth) + ",\n";
            content += "\"REC-depth-enabled\" : " + to_JSON_value(REC_depth_enabled) + ",\n";
            content += "\"REC-manual-depth-enabled\" : " + to_JSON_value(REC_manual_depth_enabled) + ",\n";
            content += "\"REC-sub-bound-percentage\" : " + to_JSON_value(REC_sub_bound_percentage) + ",\n";
            content += "\"REC-sub-bound-percentage-vec\" : " + to_JSON_value(REC_sub_bound_percentage_vec) + ",\n";
            content += "\"REC-manual-sub-bound-percentage-vec\" : " + to_JSON_value(REC_manual_sub_bound_percentage_vec) + "\n";

            return content;

        }

        /**
         * Parses the Lazy Evaluation arguments.
         *
         * @param le_option Lazy Evaluation arguments.
         * @param verbose Whether to print an error message, if an error is thrown.
         */
        void parse_LE(std::vector<std::string> &le_option, bool verbose) {
            std::string err_msg = "--LE invalid input! Input is : >>";
            for (const auto &option: le_option) {
                err_msg += option + " ";
            }
            if (!le_option.empty()) {
                err_msg.pop_back();
            }
            err_msg += "<<\nEither use:\n\t--LE Disabled\n\t--LE {Avg*y_1} {or, and} {n*y_2, k*y_2, r_n*y_2, r_k*y_2} for any y_1, y_2 >= 0\n";

            if (le_option.size() != 1 && le_option.size() != 3) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            std::string mode_string = le_option[0];
            if (mode_string == "Disabled") {
                // Disabled case
                LE_mode = 0;
                return;
            }

            // check one more element is available
            if (le_option.size() != 3) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            std::string score_threshold_string = le_option[0];

            if (contains(score_threshold_string, '*')) {
                std::vector<std::string> s = split(score_threshold_string, '*');
                std::string avg_str = s[0];
                std::string y_str = s[1];

                if (avg_str == "Avg") {
                    LE_y_score_value = std::stod(y_str);
                } else {
                    if (verbose) { std::cout << err_msg << std::endl; }
                    invalid = true;
                    return;
                }
            } else {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            std::string comb_string = le_option[1];
            if (comb_string != "or" && comb_string != "and") {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }
            if (comb_string == "or") {
                LE_mode = 1;
            } else if (comb_string == "and") {
                LE_mode = 2;
            }

            std::string rank_threshold_string = le_option[2];

            if (contains(rank_threshold_string, '*')) {
                std::vector<std::string> sub_strings = split(rank_threshold_string, '*');
                std::string var_string = sub_strings[0];
                std::string y_string = sub_strings[1];

                LE_y_rank_value = std::stod(y_string);
                if (var_string == "n") {
                    LE_rank_var = 1;
                } else if (var_string == "k") {
                    LE_rank_var = 2;
                } else if (var_string == "r_n") {
                    LE_rank_var = 3;
                } else if (var_string == "r_k") {
                    LE_rank_var = 4;
                } else {
                    if (verbose) { std::cout << err_msg << std::endl; }
                    invalid = true;
                    return;
                }
            } else {
                LE_rank_var = 0;
                LE_y_score_value = std::stod(rank_threshold_string);
            }
        }

        /**
         * Parses the Upper Bound 2D arguments.
         *
         * @param ub2d_option Upper Bound 2D arguments.
         * @param verbose Whether to print an error message, if an error is thrown.
         */
        void parse_UB2D(std::vector<std::string> &ub2d_option, bool verbose) {
            std::string err_msg = "--UB2D invalid input! Input is : >>";
            for (const auto &option: ub2d_option) {
                err_msg += option + " ";
            }
            if (!ub2d_option.empty()) {
                err_msg.pop_back();
            }

            err_msg += "<<\nUse:\n\t--UB2D Disabled\n\t--UB2D {Sqrt[x]*y, x*y, y} {l_max} {Greedy, Matching, BForce, Dynamic} {A, B, AB} {0, 1} {0, 1} {lazyStart} {lazyAdd} {lowDepthPercentage} {highDepthPercentage} {subBoundPercentage} for x in {n, k, r_n, r_k},y >= 0.0 (float), l_max >= 0 (int), 0.0 <= lazyStart <= 1.0 (float), 0.0 <= lazyAdd <= 1.0 (float), 0.0 <= lowDepthPercentage <= highDepthPercentage <= 1.0 (float), 0.0 <= subBoundPercentage <= 1.0 (float).\n";

            if (ub2d_option.size() != 1 && ub2d_option.size() != 11) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            if (ub2d_option[0] == "Disabled") {
                // Disabled case
                UB2D_enabled = false;
                return;
            }
            UB2D_enabled = true;

            // check all arguments are available
            if (ub2d_option.size() != 11) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine value for l
            std::string l_string = ub2d_option[0];
            if (contains(l_string, '*')) {
                std::vector<std::string> sub_strings = split(l_string, '*');
                std::string var_string = sub_strings[0];
                std::string y_string = sub_strings[1];

                UB2D_l_y = std::stod(y_string);

                // something complex is given
                if (var_string.find("Sqrt") != std::string::npos) {
                    // sqrt is given
                    UB2D_l_func = 1;

                    var_string.pop_back();
                    var_string.erase(0, 5);

                    if (var_string == "n") {
                        UB2D_l_var = 1;
                    } else if (var_string == "k") {
                        UB2D_l_var = 2;
                    } else if (var_string == "r_n") {
                        UB2D_l_var = 3;
                    } else if (var_string == "r_k") {
                        UB2D_l_var = 4;
                    } else {
                        if (verbose) { std::cout << err_msg << std::endl; }
                        invalid = true;
                        return;
                    }
                } else {
                    // identity is given
                    UB2D_l_func = 2;

                    if (var_string == "n") {
                        UB2D_l_var = 1;
                    } else if (var_string == "k") {
                        UB2D_l_var = 2;
                    } else if (var_string == "r_n") {
                        UB2D_l_var = 3;
                    } else if (var_string == "r_k") {
                        UB2D_l_var = 4;
                    } else {
                        if (verbose) { std::cout << err_msg << std::endl; }
                        invalid = true;
                        return;
                    }
                }
            } else {
                // only y is given
                UB2D_l_func = 3;
                UB2D_l_var = 0;
                UB2D_l_y = std::stod(l_string);
            }

            // determine max value for l
            std::string max_l_string = ub2d_option[1];
            UB2D_max_l = std::stoi(max_l_string);

            // determine the algorithm
            std::string algorithm_string = ub2d_option[2];
            if (algorithm_string == "Greedy") {
                UB2D_alg_type = 1;
            } else if (algorithm_string == "Matching") {
                UB2D_alg_type = 2;
            } else if (algorithm_string == "BForce") {
                UB2D_alg_type = 3;
            } else if (algorithm_string == "Dynamic") {
                UB2D_alg_type = 4;
            } else {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine algorithm option
            std::string option_string = ub2d_option[3];
            if (option_string == "A") {
                UB2D_odd_type = 1;
            } else if (option_string == "B") {
                UB2D_odd_type = 2;
            } else if (option_string == "AB") {
                UB2D_odd_type = 3;
            } else {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine RPC flag
            std::string rpc_string = ub2d_option[4];
            UB2D_RPC_enabled = std::stoi(rpc_string);

            // determine safe skip flag
            std::string safe_skip_string = ub2d_option[5];
            UB2D_safe_skip_enabled = std::stoi(safe_skip_string);

            // determine lazy skip start
            std::string lazy_skip_start_value_string = ub2d_option[6];
            UB2D_lazy_skip_start_value = std::stod(lazy_skip_start_value_string);
            if (UB2D_lazy_skip_start_value > 1.0 || UB2D_lazy_skip_start_value < 0) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine lazy skip add
            std::string lazy_skip_add_value_string = ub2d_option[7];
            UB2D_lazy_skip_add_value = std::stod(lazy_skip_add_value_string);
            if (UB2D_lazy_skip_add_value > 1.0 || UB2D_lazy_skip_add_value < 0) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine depth enabled lower bound
            std::string depth_enabled_lower_bound = ub2d_option[8];
            UB2D_low_depth = (size_t) (std::stod(depth_enabled_lower_bound) * (double) k);

            // determine depth enabled upper bound
            std::string depth_enabled_upper_bound = ub2d_option[9];
            UB2D_high_depth = (size_t) (std::stod(depth_enabled_upper_bound) * (double) k);

            // determine sub bound percentage
            std::string sub_bound_percentage = ub2d_option[10];
            UB2D_sub_bound_percentage = std::stod(sub_bound_percentage);
        }

        /**
         * Parses the Partial Brute Force arguments.
         *
         * @param pbf_option Partial Brute Force arguments.
         * @param verbose Whether to print an error message, if an error is thrown.
         */
        void parse_PBF(std::vector<std::string> &pbf_option, bool verbose) {
            std::string err_msg = "--PBF invalid input! Input is : >>";
            for (const auto &option: pbf_option) {
                err_msg += option + " ";
            }
            if (!pbf_option.empty()) {
                err_msg.pop_back();
            }
            err_msg += "<<\nEither use:\n\t--PBF Disabled\n\t--PBF {Sqrt[x_1]*y_1, x_1*y_1, y_1} {n_max} {Sqrt[x_2]*y_2, x_2*y_2, y_2} {l_max} {BForce, Dynamic} {0, 1} {0, 1} {lazyStart} {lazyAdd} {lowDepthPercentage} {highDepthPercentage} {subBoundPercentage} for x_1, x_2 in {n, k, r_n, r_k}, y_1, y_2 >= 0.0 (float), n_max >= 0 (int), l_max >= 0 (int), 0.0 <= lazyStart <= 1.0 (float), 0.0 <= lazyAdd <= 1.0 (float), 0.0 <= lowDepthPercentage <= highDepthPercentage <= 1.0 (float), 0.0 <= subBoundPercentage <= 1.0 (float) .\n";

            if (pbf_option.size() != 1 && pbf_option.size() != 12) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            if (pbf_option[0] == "Disabled") {
                // Disabled case
                PBF_enabled = false;
                return;
            }

            // check one more element is available
            if (pbf_option.size() != 12) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }
            PBF_enabled = true;

            // determine n func
            std::string n_string = pbf_option[0];
            if (contains(n_string, '*')) {
                std::vector<std::string> sub_strings = split(n_string, '*');
                std::string var_string = sub_strings[0];
                std::string y_string = sub_strings[1];

                PBF_n_y = std::stod(y_string);

                // something complex is given
                if (var_string.find("Sqrt") != std::string::npos) {
                    // sqrt is given
                    PBF_n_func = 1;

                    var_string.pop_back();
                    var_string.erase(0, 5);

                    if (var_string == "n") {
                        PBF_n_var = 1;
                    } else if (var_string == "k") {
                        PBF_n_var = 2;
                    } else if (var_string == "r_n") {
                        PBF_n_var = 3;
                    } else if (var_string == "r_k") {
                        PBF_n_var = 4;
                    } else {
                        if (verbose) { std::cout << err_msg << std::endl; }
                        invalid = true;
                        return;
                    }
                } else {
                    // identity is given
                    PBF_n_func = 2;

                    if (var_string == "n") {
                        PBF_n_var = 1;
                    } else if (var_string == "k") {
                        PBF_n_var = 2;
                    } else if (var_string == "r_n") {
                        PBF_n_var = 3;
                    } else if (var_string == "r_k") {
                        PBF_n_var = 4;
                    } else {
                        if (verbose) { std::cout << err_msg << std::endl; }
                        invalid = true;
                        return;
                    }
                }
            } else {
                // only y is given
                PBF_n_func = 3;
                PBF_n_var = 0;
                PBF_n_y = std::stod(n_string);
            }

            // determine max n
            std::string n_max_string = pbf_option[1];
            PBF_max_n = std::stoi(n_max_string);

            // determine l func
            std::string l_string = pbf_option[2];
            if (contains(l_string, '*')) {
                std::vector<std::string> sub_strings = split(l_string, '*');
                std::string var_string = sub_strings[0];
                std::string y_string = sub_strings[1];

                PBF_l_y = std::stod(y_string);

                // something complex is given
                if (var_string.find("Sqrt") != std::string::npos) {
                    // sqrt is given
                    PBF_l_func = 1;

                    var_string.pop_back();
                    var_string.erase(0, 5);

                    if (var_string == "n") {
                        PBF_l_var = 1;
                    } else if (var_string == "k") {
                        PBF_l_var = 2;
                    } else if (var_string == "r_n") {
                        PBF_l_var = 3;
                    } else if (var_string == "r_k") {
                        PBF_l_var = 4;
                    } else {
                        if (verbose) { std::cout << err_msg << std::endl; }
                        invalid = true;
                        return;
                    }
                } else {
                    // identity is given
                    PBF_l_func = 2;

                    if (var_string == "n") {
                        PBF_l_var = 1;
                    } else if (var_string == "k") {
                        PBF_l_var = 2;
                    } else if (var_string == "r_n") {
                        PBF_l_var = 3;
                    } else if (var_string == "r_k") {
                        PBF_l_var = 4;
                    } else {
                        if (verbose) { std::cout << err_msg << std::endl; }
                        invalid = true;
                        return;
                    }
                }
            } else {
                // only y is given
                PBF_l_func = 3;
                PBF_l_var = 0;
                PBF_l_y = std::stod(l_string);
            }

            // determine max l
            std::string l_max_string = pbf_option[3];
            PBF_max_l = std::stoi(l_max_string);

            // determine algorithm
            std::string algo_string = pbf_option[4];
            if (algo_string == "BForce") {
                PBF_alg_type = 1;
            } else if (algo_string == "Dynamic") {
                PBF_alg_type = 2;
            } else {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine RPC flag
            std::string rpc_string = pbf_option[5];
            PBF_RPC_enabled = std::stoi(rpc_string);

            // determine safe skip flag
            std::string safe_skip_flag = pbf_option[6];
            PBF_safe_skip_enabled = std::stoi(safe_skip_flag);

            // determine lazy skip start
            std::string lazy_skip_start_value_string = pbf_option[7];
            PBF_lazy_skip_start_value = std::stod(lazy_skip_start_value_string);
            if (PBF_lazy_skip_start_value > 1.0 || PBF_lazy_skip_start_value < 0) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine lazy skip add
            std::string lazy_skip_add_value_string = pbf_option[8];
            PBF_lazy_skip_add_value = std::stod(lazy_skip_add_value_string);
            if (PBF_lazy_skip_add_value > 1.0 || PBF_lazy_skip_add_value < 0) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine depth enabled lower bound
            std::string depth_enabled_lower_bound = pbf_option[9];
            PBF_low_depth = (size_t) (std::stod(depth_enabled_lower_bound) * (double) k);

            // determine depth enabled upper bound
            std::string depth_enabled_upper_bound = pbf_option[10];
            PBF_high_depth = (size_t) (std::stod(depth_enabled_upper_bound) * (double) k);

            // determine sub bound percentage
            std::string sub_bound_percentage = pbf_option[11];
            PBF_sub_bound_percentage = std::stod(sub_bound_percentage);
        }

        /**
         * Parses the Recursive arguments.
         *
         * @param rec_option Recursive arguments.
         * @param verbose Whether to print an error message, if an error is thrown.
         */
        void parse_REC(std::vector<std::string> &rec_option, bool verbose) {
            std::string err_msg = "--DAC invalid input! Input is : >>";
            for (const auto &option: rec_option) {
                err_msg += option + " ";
            }
            if (!rec_option.empty()) {
                err_msg.pop_back();
            }
            err_msg += "<<\nEither use:\n\t--DAC Disabled\n\t--DAC {Sqrt[x]*y, x*y, y} {l_max} {0, 1} {0, 1} {lazyStart} {lazyAdd} {lowDepthPercentage} {highDepthPercentage} {subBoundPercentage} for x in {n, k, r_n, r_k}, y >= 0.0 (float), l_max >= 0 (int), 0.0 <= lazyStart <= 1.0 (float), 0.0 <= lazyAdd <= 1.0 (float), 0.0 <= lowDepthPercentage <= highDepthPercentage <= 1.0 (float), 0.0 <= subBoundPercentage <= 1.0 (float) .\n";

            if (rec_option.size() != 1 && rec_option.size() != 9) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            if (rec_option[0] == "Disabled") {
                // Disabled case
                REC_enabled = false;
                return;
            }

            // check one more element is available
            if (rec_option.size() != 9) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }
            REC_enabled = true;

            // determine l func
            std::string l_string = rec_option[0];
            if (contains(l_string, '*')) {
                std::vector<std::string> sub_strings = split(l_string, '*');
                std::string var_string = sub_strings[0];
                std::string y_string = sub_strings[1];

                REC_l_y = std::stod(y_string);

                // something complex is given
                if (var_string.find("Sqrt") != std::string::npos) {
                    // sqrt is given
                    REC_l_func = 1;

                    var_string.pop_back();
                    var_string.erase(0, 5);

                    if (var_string == "n") {
                        REC_l_var = 1;
                    } else if (var_string == "k") {
                        REC_l_var = 2;
                    } else if (var_string == "r_n") {
                        REC_l_var = 3;
                    } else if (var_string == "r_k") {
                        REC_l_var = 4;
                    } else {
                        if (verbose) { std::cout << err_msg << std::endl; }
                        invalid = true;
                        return;
                    }
                } else {
                    // identity is given
                    REC_l_func = 2;

                    if (var_string == "n") {
                        REC_l_var = 1;
                    } else if (var_string == "k") {
                        REC_l_var = 2;
                    } else if (var_string == "r_n") {
                        REC_l_var = 3;
                    } else if (var_string == "r_k") {
                        REC_l_var = 4;
                    } else {
                        if (verbose) { std::cout << err_msg << std::endl; }
                        invalid = true;
                        return;
                    }
                }
            } else {
                // only y is given
                REC_l_func = 3;
                REC_l_var = 0;
                REC_l_y = std::stod(l_string);
            }

            // determine max l
            std::string l_max_string = rec_option[1];
            REC_max_l = std::stoi(l_max_string);

            // determine RPC flag
            std::string rpc_string = rec_option[2];
            REC_RPC_enabled = std::stoi(rpc_string);

            // determine safe skip flag
            std::string safe_skip_flag = rec_option[3];
            REC_safe_skip_enabled = std::stoi(safe_skip_flag);

            // determine lazy skip start
            std::string lazy_skip_start_value_string = rec_option[4];
            REC_lazy_skip_start_value = std::stod(lazy_skip_start_value_string);
            if (REC_lazy_skip_start_value > 1.0 || REC_lazy_skip_start_value < 0) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine lazy skip add
            std::string lazy_skip_add_value_string = rec_option[5];
            REC_lazy_skip_add_value = std::stod(lazy_skip_add_value_string);
            if (REC_lazy_skip_add_value > 1.0 || REC_lazy_skip_add_value < 0) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            // determine depth enabled lower bound
            std::string depth_enabled_lower_bound = rec_option[6];
            REC_low_depth = (size_t) (std::stod(depth_enabled_lower_bound) * (double) k);

            // determine depth enabled upper bound
            std::string depth_enabled_upper_bound = rec_option[7];
            REC_high_depth = (size_t) (std::stod(depth_enabled_upper_bound) * (double) k);

            // determine sub bound percentage
            std::string sub_bound_percentage = rec_option[8];
            REC_sub_bound_percentage = std::stod(sub_bound_percentage);
        }

        /**
         * Parses the Brute Force Threshold arguments.
         *
         * @param bf_threshold_option Brute Force Threshold arguments.
         * @param verbose Whether to print an error message, if an error is thrown.
         */
        void parse_BFThreshold(std::vector<std::string> &bf_threshold_option, bool verbose) {
            std::string err_msg = "--BF Threshold invalid input! Input is : >>";
            for (const auto &option: bf_threshold_option) {
                err_msg += option + " ";
            }
            if (!bf_threshold_option.empty()) {
                err_msg.pop_back();
            }
            err_msg += "<<\nUse:\n\t--BFThreshold {n} {r} for 0 <= r <= n.\n";

            if (bf_threshold_option.size() != 2) {
                if (verbose) { std::cout << err_msg << std::endl; }
                invalid = true;
                return;
            }

            bf_threshold_n = std::stoi(bf_threshold_option[0]);
            bf_threshold_r = std::stoi(bf_threshold_option[1]);
        }
    };

    /**
     * Gets all valid configurations.
     *
     * @return Vector containing all configurations.
     */
    std::vector<AlgorithmConfiguration> get_all_algorithm_configurations();

    /**
     * Return the fastest (known) algorithm configuration.
     *
     * @return The algorithm configuration.
     */
    AlgorithmConfiguration get_fast_algorithm_configuration();

    /**
 * Parses the command line and return the options for the algorithm.
 *
 * @param argc Number of arguments.
 * @param argv Array to the arguments.
 * @param verbose If the program prints error messages.
 * @return Algorithm configuration
 */
    AlgorithmConfiguration parse_command_line(int argc, char *argv[], bool verbose = false);

}

#endif //SUBSETOPTIMIZATION_ALGORITHMCONFIGURATION_H
