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

#include <gtest/gtest.h>

#include "../src/utility/Util.h"
#include "../src/utility/AlgorithmConfiguration.h"
#include "../src/algorithms/UB2DAlgorithm.h"

namespace SMSM {

/**
 * Tests basic options for the Command Line parser.
 */
    TEST(Misc, CommandLineParsing_Basic) {
        AlgorithmConfiguration ac;
        std::string s;
        char **argv = nullptr;
        int argc = 0;

        s = "subsetoptimization";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        EXPECT_TRUE(ac.structure_type == "graph");
        free_argv(argc, argv);

        s = "subsetoptimization --type graph";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        EXPECT_TRUE(ac.structure_type == "graph");
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        EXPECT_TRUE(ac.input_file_path == "example.txt");
        free_argv(argc, argv);

        s = "subsetoptimization -t graph --input-file example.txt";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        EXPECT_TRUE(ac.input_file_path == "example.txt");
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        EXPECT_EQ(ac.k, 10);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_TRUE(ac.score_function == "example-function");
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 --score-function example-function";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_TRUE(ac.score_function == "example-function");
        EXPECT_EQ(ac.write_output, false);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function -o out.JSON";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.write_output, true);
        EXPECT_TRUE(ac.output_file_path == "out.JSON");
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --output-file out.JSON";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.write_output, true);
        EXPECT_TRUE(ac.output_file_path == "out.JSON");
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --time-limit 10";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_DOUBLE_EQ(ac.time_limit, 10.0);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --time-limit 0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_DOUBLE_EQ(ac.time_limit, std::numeric_limits<double>::max());
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --SUB 1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.SUB_enabled, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --CR 1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.RPC_enabled, true);
        free_argv(argc, argv);
    }


/**
 * Tests LE options for the Command Line parser.
 */
    TEST(Misc, CommandLineParsing_LE) {
        AlgorithmConfiguration ac;
        std::string s;
        char **argv = nullptr;
        int argc = 0;

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE 1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE Disabled";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.LE_mode, 0);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE Score";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE Rank";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE INVALID";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE Avg";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE Avg*1 or n*2";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.LE_mode, 1);
        EXPECT_DOUBLE_EQ(ac.LE_y_score_value, 1.0);
        EXPECT_EQ(ac.LE_rank_var, 1);
        EXPECT_DOUBLE_EQ(ac.LE_y_rank_value, 2.0);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE Avg*0.5 and k*0.8";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.LE_mode, 2);
        EXPECT_DOUBLE_EQ(ac.LE_y_score_value, 0.5);
        EXPECT_EQ(ac.LE_rank_var, 2);
        EXPECT_DOUBLE_EQ(ac.LE_y_rank_value, 0.8);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE Avg*2.5 and r_n*0.6";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.LE_mode, 2);
        EXPECT_DOUBLE_EQ(ac.LE_y_score_value, 2.5);
        EXPECT_EQ(ac.LE_rank_var, 3);
        EXPECT_DOUBLE_EQ(ac.LE_y_rank_value, 0.6);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --LE Avg*1.5 or r_k*1.8";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.LE_mode, 1);
        EXPECT_DOUBLE_EQ(ac.LE_y_score_value, 1.5);
        EXPECT_EQ(ac.LE_rank_var, 4);
        EXPECT_DOUBLE_EQ(ac.LE_y_rank_value, 1.8);
        free_argv(argc, argv);
    }


/**
 * Tests UB2D options for the Command Line parser.
 */
    TEST(Misc, CommandLineParsing_UB2D) {
        AlgorithmConfiguration ac;
        std::string s;
        char **argv = nullptr;
        int argc = 0;

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D 1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Disabled";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, false);
        free_argv(argc, argv);

        // Default configuration
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_l_func, 1);
        EXPECT_EQ(ac.UB2D_l_var, 1);
        EXPECT_DOUBLE_EQ(ac.UB2D_l_y, 1.0);
        EXPECT_EQ(ac.UB2D_max_l, 20);
        EXPECT_EQ(ac.UB2D_alg_type, 1);
        EXPECT_EQ(ac.UB2D_odd_type, 1);
        EXPECT_EQ(ac.UB2D_RPC_enabled, 1);
        EXPECT_EQ(ac.UB2D_safe_skip_enabled, 1);
        EXPECT_DOUBLE_EQ(ac.UB2D_lazy_skip_start_value, 1.0);
        EXPECT_DOUBLE_EQ(ac.UB2D_lazy_skip_add_value, 1.0);
        EXPECT_EQ(ac.UB2D_low_depth, 2);
        EXPECT_EQ(ac.UB2D_high_depth, 8);
        EXPECT_DOUBLE_EQ(ac.UB2D_sub_bound_percentage, 1.0);
        free_argv(argc, argv);

        // try different l_var
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[k]*1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_l_var, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[r_n]*1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_l_var, 3);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[r_k]*1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_l_var, 4);
        free_argv(argc, argv);

        // try different l_func
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D n*1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_l_func, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D 1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_l_func, 3);
        free_argv(argc, argv);

        // try different l_y
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*2 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_l_y, 2.0);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D r_k*0.5 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_l_y, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D 10.1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_l_y, 10.1);
        free_argv(argc, argv);

        // try different max l
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 5 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_max_l, 5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 100 Greedy A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_max_l, 100);
        free_argv(argc, argv);

        // try different algorithm
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Matching A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_alg_type, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 BForce A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_alg_type, 3);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Dynamic A 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_alg_type, 4);
        free_argv(argc, argv);

        // try different odd types enabled
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy B 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_odd_type, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy AB 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_odd_type, 3);
        free_argv(argc, argv);

        // try different RPC enabled
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 0 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_RPC_enabled, 0);
        free_argv(argc, argv);

        // try different safe skip enabled
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 0 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_safe_skip_enabled, 0);
        free_argv(argc, argv);

        // try different lazy skip start values
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 0.5 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_lazy_skip_start_value, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 0.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_lazy_skip_start_value, 0.0);
        free_argv(argc, argv);

        // try different lazy skip add values
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 0.5 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_lazy_skip_add_value, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 0.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_lazy_skip_add_value, 0.0);
        free_argv(argc, argv);

        // try different low depths
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 1.0 0.4 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_low_depth, 4);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 1.0 0.1 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_low_depth, 1);
        free_argv(argc, argv);

        // try different high depths
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 1.0 0.2 0.3 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_high_depth, 3);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 1.0 0.2 0.6 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_EQ(ac.UB2D_high_depth, 6);
        free_argv(argc, argv);

        // try different sub bound percentages
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 0.5";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_sub_bound_percentage, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --UB2D Sqrt[n]*1 20 Greedy A 1 1 1.0 1.0 0.2 0.8 0.1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.UB2D_enabled, true);
        EXPECT_DOUBLE_EQ(ac.UB2D_sub_bound_percentage, 0.1);
        free_argv(argc, argv);
    }


/**
 * Tests PBF options for the Command Line parser.
 */
    TEST(Misc, CommandLineParsing_PBF) {
        AlgorithmConfiguration ac;
        std::string s;
        char **argv = nullptr;
        int argc = 0;

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF 1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Disabled";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, false);
        free_argv(argc, argv);

        // Default configuration
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_n_func, 1);
        EXPECT_EQ(ac.PBF_n_var, 1);
        EXPECT_DOUBLE_EQ(ac.PBF_n_y, 1.0);
        EXPECT_EQ(ac.PBF_max_n, 20);
        EXPECT_EQ(ac.PBF_l_func, 2);
        EXPECT_EQ(ac.PBF_l_var, 2);
        EXPECT_DOUBLE_EQ(ac.PBF_l_y, 2.0);
        EXPECT_EQ(ac.PBF_max_l, 10);
        EXPECT_EQ(ac.PBF_alg_type, 1);
        EXPECT_EQ(ac.PBF_RPC_enabled, 1);
        EXPECT_EQ(ac.PBF_safe_skip_enabled, 1);
        EXPECT_DOUBLE_EQ(ac.PBF_lazy_skip_start_value, 1.0);
        EXPECT_DOUBLE_EQ(ac.PBF_lazy_skip_add_value, 1.0);
        EXPECT_EQ(ac.PBF_low_depth, 2);
        EXPECT_EQ(ac.PBF_high_depth, 8);
        EXPECT_DOUBLE_EQ(ac.PBF_sub_bound_percentage, 1.0);
        free_argv(argc, argv);

        // try different n_funcs
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF n*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_n_func, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF 1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_n_func, 3);
        free_argv(argc, argv);

        // try different n_vars
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[k]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_n_var, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[r_n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_n_var, 3);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[r_k]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_n_var, 4);
        free_argv(argc, argv);

        // try different n_ys
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*0.5 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_n_y, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF 2.5 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_n_y, 2.5);
        free_argv(argc, argv);

        // try different max_n
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 10 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_max_n, 10);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 1000 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_max_n, 1000);
        free_argv(argc, argv);

        // try different l_funcs
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 Sqrt[k]*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_l_func, 1);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_l_func, 3);
        free_argv(argc, argv);

        // try different l_vars
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 n*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_l_var, 1);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 r_n*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_l_var, 3);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 r_k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_l_var, 4);
        free_argv(argc, argv);

        // try different l_ys
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*3.1 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_l_y, 3.1);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 3.7 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_l_y, 3.7);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*0.4 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_l_y, 0.4);
        free_argv(argc, argv);

        // try different max_ls
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 1 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_max_l, 1);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 25 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_max_l, 25);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 1000 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_max_l, 1000);
        free_argv(argc, argv);

        // try different algos
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 Dynamic 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_alg_type, 2);
        free_argv(argc, argv);

        // try different RPC enabled
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 0 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_RPC_enabled, 0);
        free_argv(argc, argv);

        // try different safe skip enabled
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 0 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_EQ(ac.PBF_safe_skip_enabled, 0);
        free_argv(argc, argv);

        // try different lazy skip start values
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 0.5 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_lazy_skip_start_value, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 0.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_lazy_skip_start_value, 0.0);
        free_argv(argc, argv);

        // try different lazy skip add values
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 0.5 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_lazy_skip_add_value, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 0.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_lazy_skip_add_value, 0.0);
        free_argv(argc, argv);

        // try different low depths
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_low_depth, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.0 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_low_depth, 0);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.7 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_low_depth, 7);
        free_argv(argc, argv);

        // try different high depths
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.3 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_high_depth, 3);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.5 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_high_depth, 5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.9 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_high_depth, 9);
        free_argv(argc, argv);

        // try different sub percentages
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 0.5";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_sub_bound_percentage, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --PBF Sqrt[n]*1 20 k*2 10 BForce 1 1 1.0 1.0 0.2 0.8 0.2";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.PBF_enabled, true);
        EXPECT_DOUBLE_EQ(ac.PBF_sub_bound_percentage, 0.2);
        free_argv(argc, argv);
    }

/**
 * Tests REC options for the Command Line parser.
 */
    TEST(Misc, CommandLineParsing_REC) {
        AlgorithmConfiguration ac;
        std::string s;
        char **argv = nullptr;
        int argc = 0;

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC 1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, true);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC Disabled";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, false);
        free_argv(argc, argv);

        // Default configuration
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC Sqrt[n]*1 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_l_func, 1);
        EXPECT_EQ(ac.REC_l_var, 1);
        EXPECT_DOUBLE_EQ(ac.REC_l_y, 1.0);
        EXPECT_EQ(ac.REC_max_l, 20);
        EXPECT_EQ(ac.REC_RPC_enabled, 1);
        EXPECT_EQ(ac.REC_safe_skip_enabled, 1);
        EXPECT_DOUBLE_EQ(ac.REC_lazy_skip_start_value, 1.0);
        EXPECT_DOUBLE_EQ(ac.REC_lazy_skip_add_value, 1.0);
        EXPECT_EQ(ac.REC_low_depth, 2);
        EXPECT_EQ(ac.REC_high_depth, 8);
        EXPECT_DOUBLE_EQ(ac.REC_sub_bound_percentage, 1.0);
        free_argv(argc, argv);

        // try different l_funcs
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC n*1 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_l_func, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC 1 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_l_func, 3);
        free_argv(argc, argv);

        // try different l_vars
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*1 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_l_var, 2);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC r_n*1 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_l_var, 3);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC r_k*1 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_l_var, 4);
        free_argv(argc, argv);

        // try different l_ys
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_l_y, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*2.5 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_l_y, 2.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*10.1 20 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_l_y, 10.1);
        free_argv(argc, argv);

        // try different max_ls
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 5 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_max_l, 5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_max_l, 500);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 5000 1 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_max_l, 5000);
        free_argv(argc, argv);

        // try different RPC enabled
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 0 1 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_RPC_enabled, 0);
        free_argv(argc, argv);

        // try different safe skip enabled
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 0 1.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_safe_skip_enabled, 0);
        free_argv(argc, argv);

        // try different lazy skip start value
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 0.5 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_lazy_skip_start_value, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 0.0 1.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_lazy_skip_start_value, 0.0);
        free_argv(argc, argv);

        // try different lazy skip add value
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 0.5 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_lazy_skip_add_value, 0.5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 0.0 0.2 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_lazy_skip_add_value, 0.0);
        free_argv(argc, argv);

        // try different low depths
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.3 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_low_depth, 3);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.5 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_low_depth, 5);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.0 0.8 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_low_depth, 0);
        free_argv(argc, argv);

        // try different high depths
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.2 1.0 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_high_depth, 10);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.2 0.7 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_high_depth, 7);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.2 0.3 1.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_high_depth, 3);
        free_argv(argc, argv);

        // try different sub percentages
        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.2 0.8 0.8";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_sub_bound_percentage, 0.8);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.2 0.8 0.1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_sub_bound_percentage, 0.1);
        free_argv(argc, argv);

        s = "subsetoptimization -t graph -i example.txt -k 10 -s example-function --DAC k*0.5 500 1 1 1.0 1.0 0.2 0.8 0.0";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_sub_bound_percentage, 0.0);
        free_argv(argc, argv);
    }

    TEST(Misc, CommandLineParsing_MISC) {
        AlgorithmConfiguration ac;
        std::string s;
        char **argv = nullptr;
        int argc = 0;

        s = "subsetoptimization -t k-medoid -i data/clustering/cs.uef.fi_sipu_datasets_dim064.txt -s manhattan-distance --SUB 1 --CR 1 --DAC r_n*0.5 100000 1 1 0.8 0.4 0.1 0.9 0.95 -k 2 -o results/ara/1793_2_profile.JSON --time-limit 1";
        argv = string_to_argv(s, argc);
        ac = parse_command_line(argc, argv);
        EXPECT_EQ(ac.invalid, false);
        EXPECT_EQ(ac.structure_type, "k-medoid");
        EXPECT_EQ(ac.input_file_path, "data/clustering/cs.uef.fi_sipu_datasets_dim064.txt");
        EXPECT_EQ(ac.score_function, "manhattan-distance");
        EXPECT_EQ(ac.k, 2);
        EXPECT_EQ(ac.SUB_enabled, true);
        EXPECT_EQ(ac.RPC_enabled, true);
        EXPECT_EQ(ac.REC_enabled, true);
        EXPECT_EQ(ac.REC_l_func, 2);
        EXPECT_EQ(ac.REC_l_var, 3);
        EXPECT_DOUBLE_EQ(ac.REC_l_y, 0.5);
        EXPECT_EQ(ac.REC_max_l, 100000);
        EXPECT_EQ(ac.REC_RPC_enabled, true);
        EXPECT_EQ(ac.REC_safe_skip_enabled, true);
        EXPECT_DOUBLE_EQ(ac.REC_lazy_skip_start_value, 0.8);
        EXPECT_DOUBLE_EQ(ac.REC_lazy_skip_add_value, 0.4);
        EXPECT_DOUBLE_EQ(ac.REC_sub_bound_percentage, 0.95);
        std::cout << ac.to_JSON() << std::endl;
        free_argv(argc, argv);
    }

}
