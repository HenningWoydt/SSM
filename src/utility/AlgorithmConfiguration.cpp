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

#include "AlgorithmConfiguration.h"

namespace SMSM {

    namespace boost_po = boost::program_options;


    AlgorithmConfiguration parse_command_line(int argc, char *argv[], bool verbose) {
        AlgorithmConfiguration ac = AlgorithmConfiguration();

        boost_po::options_description description("Allowed options");
        description.add_options()
                ("help,h", "Produces a help message")
                ("type,t", boost_po::value<std::string>(&ac.structure_type), "Which structure is given (important for reading the input file)")
                ("input-file,i", boost_po::value<std::string>(&ac.input_file_path), "Path to the file holding the structure")
                ("k,k", boost_po::value<size_t>(&ac.k), "Size of the set.")
                ("score-function,s", boost_po::value<std::string>(&ac.score_function), "The function to maximize")
                ("output-file,o", boost_po::value<std::string>(&ac.output_file_path), "Path to the output file")
                ("time-limit", boost_po::value<double>(&ac.time_limit), "Time-limit in seconds (0 == infinite)")
                ("SUB", boost_po::value<bool>(&ac.SUB_enabled)->default_value(false), "Enables Simple-Upper-Bound heuristic")
                ("CR", boost_po::value<bool>(&ac.RPC_enabled)->default_value(false), "Enables Candidate-Reduction heuristic")
                ("Greedy-Local-Search", boost_po::value<bool>(&ac.greedy_local_search)->default_value(false), "Whether to enable Greedy Local Search")
                ("LE", boost_po::value<std::vector<std::string>>()->multitoken(), "Enables Lazy-Evaluation")
                ("UB2D", boost_po::value<std::vector<std::string>>()->multitoken(), "Enables Upper-Bound-2D-Dynamic heuristic")
                ("PBF", boost_po::value<std::vector<std::string>>()->multitoken(), "Enables Partial-Brute-Force heuristic")
                ("DAC", boost_po::value<std::vector<std::string>>()->multitoken(), "Enables Divide-And-Conquer heuristic")
                ("BFThreshold", boost_po::value<std::vector<std::string>>()->multitoken(), "Sets the Brute Force Threshold")
                ("Plain", boost_po::value<bool>(&ac.plain)->default_value(false), "Whether to check for the plain algorithm")
                ("measure-oracle", boost_po::value<bool>(&ac.measure_oracle_time)->default_value(false), "Whether to measure oracle time");

        boost_po::variables_map vm;
        boost_po::store(boost_po::parse_command_line(argc, argv, description), vm);
        boost_po::notify(vm);

        if (vm.count("help")) {
            std::cout << description << "\n";
            exit(EXIT_SUCCESS);
        }

        // check required arguments
        if (!vm.count("type") && !vm.count("t")) {
            if (verbose) {
                std::cout << "-t [  --type  ] not specified\n";
            }
            ac.invalid = true;
        }
        if (!vm.count("input-file") && !vm.count("i")) {
            if (verbose) {
                std::cout << "-i [  --input-file  ] not specified\n";
            }
            ac.invalid = true;
        }
        if (!vm.count("k")) {
            if (verbose) {
                std::cout << "-k [ --k  ] not specified\n";
            }
            ac.invalid = true;
        }
        if (!vm.count("score-function") && !vm.count("s")) {
            if (verbose) {
                std::cout << "-s [ --score-function  ] not specified\n";
            }
            ac.invalid = true;
        }

        // check optional arguments
        if (vm.count("o") || vm.count("output-file")) {
            ac.write_output = true;
        }

        if (vm.count("time-limit")) {
            if (vm["time-limit"].as<double>() == 0.0) {
                ac.time_limit = std::numeric_limits<double>::max();
            }
        } else {
            ac.time_limit = std::numeric_limits<double>::max();
        }

        if (ac.plain) {
            ac.bf_threshold_n = 1;
            ac.bf_threshold_r = 1;
            ac.SUB_enabled = false;
            ac.RPC_enabled = false;
            ac.LE_mode = 0;
            ac.UB2D_enabled = false;
            ac.PBF_enabled = false;
            ac.REC_enabled = false;

            ac.finalize();
            return ac;
        }

        // check complex heuristics
        if (vm.count("LE")) {
            std::vector<std::string> le_option = vm["LE"].as<std::vector<std::string>>();
            ac.parse_LE(le_option, verbose);
        }

        if (vm.count("UB2D")) {
            std::vector<std::string> ub2d_option = vm["UB2D"].as<std::vector<std::string>>();
            ac.parse_UB2D(ub2d_option, verbose);
        }

        if (vm.count("PBF")) {
            std::vector<std::string> pbf_option = vm["PBF"].as<std::vector<std::string>>();
            ac.parse_PBF(pbf_option, verbose);
        }

        if (vm.count("DAC")) {
            std::vector<std::string> rec_option = vm["DAC"].as<std::vector<std::string>>();
            ac.parse_REC(rec_option, verbose);
        }

        if (vm.count("BFThreshold")) {
            std::vector<std::string> bf_threshold_option = vm["BFThreshold"].as<std::vector<std::string>>();
            ac.parse_BFThreshold(bf_threshold_option, verbose);
        }

        ac.finalize();

        return ac;
    }

    std::vector<AlgorithmConfiguration> get_all_algorithm_configurations() {
        std::vector<AlgorithmConfiguration> configs;

        // SUB_enabled: test enabled
        std::vector<bool> SUB_enabled = {true};
        for (auto SUB: SUB_enabled) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.SUB_enabled = SUB;
            ac.finalize();
            configs.push_back(ac);
        }

        // RPC_enabled: test enabled
        std::vector<bool> RPC_enabled = {true};
        for (auto RPC: RPC_enabled) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.RPC_enabled = RPC;
            ac.finalize();
            configs.push_back(ac);
        }

        // UB2D: test different values of l (when using a function)
        std::vector<size_t> UB2D_l_func = {1, 2};
        std::vector<size_t> UB2D_l_var = {1, 2, 3, 4};
        std::vector<double> UB2D_l_y = {0.5, 1.0, 2.0};
        for (auto l_func: UB2D_l_func) {
            for (auto l_var: UB2D_l_var) {
                for (auto l_y: UB2D_l_y) {
                    AlgorithmConfiguration ac = AlgorithmConfiguration();
                    ac.UB2D_enabled = true;
                    ac.UB2D_l_func = l_func;
                    ac.UB2D_l_var = l_var;
                    ac.UB2D_l_y = l_y;
                    ac.UB2D_max_l = 10;
                    ac.UB2D_alg_type = 2;
                    ac.UB2D_odd_type = 1;
                    ac.UB2D_RPC_enabled = false;
                    ac.UB2D_safe_skip_enabled = false;
                    ac.UB2D_lazy_skip_start_value = 1.0;
                    ac.UB2D_lazy_skip_add_value = 1.0;
                    ac.UB2D_low_depth = 0;
                    ac.UB2D_high_depth = std::numeric_limits<size_t>::max();
                    ac.UB2D_sub_bound_percentage = 0.0;
                    ac.finalize();
                    configs.push_back(ac);
                }
            }
        }

        // UB2D: test different values of l, when not function specific
        std::vector<double> UB2D_l_y_2 = {2.0, 5.0, 10.0};
        for (auto l_y: UB2D_l_y_2) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.UB2D_enabled = true;
            ac.UB2D_l_func = 3;
            ac.UB2D_l_var = 0;
            ac.UB2D_l_y = l_y;
            ac.UB2D_max_l = 10;
            ac.UB2D_alg_type = 2;
            ac.UB2D_odd_type = 1;
            ac.UB2D_RPC_enabled = false;
            ac.UB2D_safe_skip_enabled = false;
            ac.UB2D_lazy_skip_start_value = 1.0;
            ac.UB2D_lazy_skip_add_value = 1.0;
            ac.UB2D_low_depth = 0;
            ac.UB2D_high_depth = std::numeric_limits<size_t>::max();
            ac.UB2D_sub_bound_percentage = 0.0;
            ac.finalize();
            configs.push_back(ac);
        }

        // UB2D: test different algorithm options
        std::vector<size_t> UB2D_alg_types = {1, 2, 3, 4};
        std::vector<size_t> UB2D_odd_types = {1, 2, 3};
        for (auto alg_type: UB2D_alg_types) {
            for (auto odd_type: UB2D_odd_types) {
                AlgorithmConfiguration ac = AlgorithmConfiguration();
                ac.UB2D_enabled = true;
                ac.UB2D_l_func = 3;
                ac.UB2D_l_var = 0;
                ac.UB2D_l_y = 10;
                ac.UB2D_max_l = 10;
                ac.UB2D_alg_type = alg_type;
                ac.UB2D_odd_type = odd_type;
                ac.UB2D_RPC_enabled = false;
                ac.UB2D_safe_skip_enabled = false;
                ac.UB2D_lazy_skip_start_value = 1.0;
                ac.UB2D_lazy_skip_add_value = 1.0;
                ac.UB2D_low_depth = 0;
                ac.UB2D_high_depth = std::numeric_limits<size_t>::max();
                ac.UB2D_sub_bound_percentage = 0.0;
                ac.finalize();
                configs.push_back(ac);
            }
        }

        // UB2D: test different values for rpc
        std::vector<bool> UB2D_rpc_enabled = {true};
        for (auto rpc_enabled: UB2D_rpc_enabled) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.UB2D_enabled = true;
            ac.UB2D_l_func = 3;
            ac.UB2D_l_var = 0;
            ac.UB2D_l_y = 10;
            ac.UB2D_max_l = 10;
            ac.UB2D_alg_type = 2;
            ac.UB2D_odd_type = 1;
            ac.UB2D_RPC_enabled = rpc_enabled;
            ac.UB2D_safe_skip_enabled = false;
            ac.UB2D_lazy_skip_start_value = 1.0;
            ac.UB2D_lazy_skip_add_value = 1.0;
            ac.UB2D_low_depth = 0;
            ac.UB2D_high_depth = std::numeric_limits<size_t>::max();
            ac.UB2D_sub_bound_percentage = 0.0;
            ac.finalize();
            configs.push_back(ac);
        }

        // UB2D: test different values for safe skip
        std::vector<bool> UB2D_safe_skip_enabled = {true};
        for (auto safe_skip: UB2D_safe_skip_enabled) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.UB2D_enabled = true;
            ac.UB2D_l_func = 3;
            ac.UB2D_l_var = 0;
            ac.UB2D_l_y = 10;
            ac.UB2D_max_l = 10;
            ac.UB2D_alg_type = 2;
            ac.UB2D_odd_type = 1;
            ac.UB2D_RPC_enabled = false;
            ac.UB2D_safe_skip_enabled = safe_skip;
            ac.UB2D_lazy_skip_start_value = 1.0;
            ac.UB2D_lazy_skip_add_value = 1.0;
            ac.UB2D_low_depth = 0;
            ac.UB2D_high_depth = std::numeric_limits<size_t>::max();
            ac.UB2D_sub_bound_percentage = 0.0;
            ac.finalize();
            configs.push_back(ac);
        }

        // UB2D: test different values for lazy skip
        std::vector<double> UB2D_lazy_skip_start_values = {0.0, 0.25, 0.5, 0.75, 1.0};
        std::vector<double> UB2D_lazy_skip_add_values = {0.1, 0.3, 0.5, 0.7, 0.9};
        for (auto add_val: UB2D_lazy_skip_add_values) {
            for (auto start_val: UB2D_lazy_skip_start_values) {
                AlgorithmConfiguration ac = AlgorithmConfiguration();
                ac.UB2D_enabled = true;
                ac.UB2D_l_func = 3;
                ac.UB2D_l_var = 0;
                ac.UB2D_l_y = 10;
                ac.UB2D_max_l = 10;
                ac.UB2D_alg_type = 2;
                ac.UB2D_odd_type = 1;
                ac.UB2D_RPC_enabled = false;
                ac.UB2D_safe_skip_enabled = false;
                ac.UB2D_lazy_skip_start_value = start_val;
                ac.UB2D_lazy_skip_add_value = add_val;
                ac.UB2D_low_depth = 0;
                ac.UB2D_high_depth = std::numeric_limits<size_t>::max();
                ac.UB2D_sub_bound_percentage = 0.0;
                ac.finalize();
                configs.push_back(ac);
            }
        }

        // UB2D: test different values for depths
        std::vector<size_t> UB2D_low_depths = {0, 1, 2, 3, 4};
        std::vector<size_t> UB2D_high_depths = {0, 4, 8, 12, 16, 20};
        for (auto low_d: UB2D_low_depths) {
            for (auto high_d: UB2D_high_depths) {
                if (low_d <= high_d) {
                    AlgorithmConfiguration ac = AlgorithmConfiguration();
                    ac.UB2D_enabled = true;
                    ac.UB2D_l_func = 3;
                    ac.UB2D_l_var = 0;
                    ac.UB2D_l_y = 10;
                    ac.UB2D_max_l = 10;
                    ac.UB2D_alg_type = 2;
                    ac.UB2D_odd_type = 1;
                    ac.UB2D_RPC_enabled = false;
                    ac.UB2D_safe_skip_enabled = false;
                    ac.UB2D_lazy_skip_start_value = 1.0;
                    ac.UB2D_lazy_skip_add_value = 1.0;
                    ac.UB2D_low_depth = low_d;
                    ac.UB2D_high_depth = high_d;
                    ac.UB2D_sub_bound_percentage = 0.0;
                    ac.finalize();
                    configs.push_back(ac);
                }
            }
        }

        // UB2D: test different values for sub percentage
        std::vector<double> UB2D_sub_percentages = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
        for (auto sub_p: UB2D_sub_percentages) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.SUB_enabled = true;
            ac.UB2D_enabled = true;
            ac.UB2D_l_func = 3;
            ac.UB2D_l_var = 0;
            ac.UB2D_l_y = 10;
            ac.UB2D_max_l = 10;
            ac.UB2D_alg_type = 2;
            ac.UB2D_odd_type = 1;
            ac.UB2D_RPC_enabled = false;
            ac.UB2D_safe_skip_enabled = false;
            ac.UB2D_lazy_skip_start_value = 1.0;
            ac.UB2D_lazy_skip_add_value = 1.0;
            ac.UB2D_low_depth = 0;
            ac.UB2D_high_depth = std::numeric_limits<size_t>::max();
            ac.UB2D_sub_bound_percentage = sub_p;
            ac.finalize();
            configs.push_back(ac);
        }

        // Partial Brute Force: test different values for n_b, n_l
        std::vector<double> PBF_n_y = {1, 2, 3};
        std::vector<double> PBF_l_y = {2, 3, 4, 5};
        std::vector<size_t> PBF_algo_type = {1, 2};
        for (auto n_y: PBF_n_y) {
            for (auto l_y: PBF_l_y) {
                for (auto alg_type: PBF_algo_type) {
                    AlgorithmConfiguration ac = AlgorithmConfiguration();
                    ac.PBF_enabled = true;
                    ac.PBF_n_func = 3;
                    ac.PBF_n_var = 0;
                    ac.PBF_n_y = n_y;
                    ac.PBF_l_func = 5;
                    ac.PBF_l_var = 0;
                    ac.PBF_l_y = l_y;
                    ac.PBF_alg_type = alg_type;
                    ac.PBF_RPC_enabled = false;
                    ac.PBF_safe_skip_enabled = false;
                    ac.PBF_lazy_skip_start_value = 1.0;
                    ac.PBF_lazy_skip_add_value = 1.0;
                    ac.PBF_low_depth = 0;
                    ac.PBF_high_depth = std::numeric_limits<size_t>::max();
                    ac.PBF_sub_bound_percentage = 0.0;
                    ac.finalize();
                    configs.push_back(ac);
                }
            }
        }

        // PBF: test rpc enabled
        std::vector<bool> PBF_rpc_enabled = {true};
        for (auto rpc_enabled: PBF_rpc_enabled) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.PBF_enabled = true;
            ac.PBF_n_func = 3;
            ac.PBF_n_var = 0;
            ac.PBF_n_y = 2;
            ac.PBF_l_func = 5;
            ac.PBF_l_var = 0;
            ac.PBF_l_y = 3;
            ac.PBF_alg_type = 2;
            ac.PBF_RPC_enabled = rpc_enabled;
            ac.PBF_safe_skip_enabled = false;
            ac.PBF_lazy_skip_start_value = 1.0;
            ac.PBF_lazy_skip_add_value = 1.0;
            ac.PBF_low_depth = 0;
            ac.PBF_high_depth = std::numeric_limits<size_t>::max();
            ac.PBF_sub_bound_percentage = 0.0;
            ac.finalize();
            configs.push_back(ac);
        }

        // PBF: test safe skip enabled
        std::vector<bool> PBF_safe_skip_enabled = {true};
        for (auto safe_skip_enabled: PBF_safe_skip_enabled) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.PBF_enabled = true;
            ac.PBF_n_func = 3;
            ac.PBF_n_var = 0;
            ac.PBF_n_y = 2;
            ac.PBF_l_func = 5;
            ac.PBF_l_var = 0;
            ac.PBF_l_y = 3;
            ac.PBF_alg_type = 2;
            ac.PBF_RPC_enabled = false;
            ac.PBF_safe_skip_enabled = safe_skip_enabled;
            ac.PBF_lazy_skip_start_value = 1.0;
            ac.PBF_lazy_skip_add_value = 1.0;
            ac.PBF_low_depth = 0;
            ac.PBF_high_depth = std::numeric_limits<size_t>::max();
            ac.PBF_sub_bound_percentage = 0.0;
            ac.finalize();
            configs.push_back(ac);
        }

        // Partial Brute Force: test different values for n_b, n_l
        std::vector<double> PBF_lazy_skip_start_values = {0.0, 0.25, 0.5, 0.75, 1.0};
        std::vector<double> PBF_lazy_skip_add_values = {0.1, 0.3, 0.5, 0.7, 0.9};
        for (auto add_value: PBF_lazy_skip_add_values) {
            for (auto start_value: PBF_lazy_skip_start_values) {
                AlgorithmConfiguration ac = AlgorithmConfiguration();
                ac.PBF_enabled = true;
                ac.PBF_n_func = 3;
                ac.PBF_n_var = 0;
                ac.PBF_n_y = 2;
                ac.PBF_l_func = 5;
                ac.PBF_l_var = 0;
                ac.PBF_l_y = 3;
                ac.PBF_alg_type = 2;
                ac.PBF_RPC_enabled = false;
                ac.PBF_safe_skip_enabled = false;
                ac.PBF_lazy_skip_start_value = start_value;
                ac.PBF_lazy_skip_add_value = add_value;
                ac.PBF_low_depth = 0;
                ac.PBF_high_depth = std::numeric_limits<size_t>::max();
                ac.PBF_sub_bound_percentage = 0.0;
                ac.finalize();
                configs.push_back(ac);
            }
        }

        // PBF: test different values for depths
        std::vector<size_t> PBF_low_depths = {0, 1, 2, 3, 4};
        std::vector<size_t> PBF_high_depths = {0, 4, 8, 12, 16, 20};
        for (auto low_d: PBF_low_depths) {
            for (auto high_d: PBF_high_depths) {
                if (low_d <= high_d) {
                    AlgorithmConfiguration ac = AlgorithmConfiguration();
                    ac.PBF_enabled = true;
                    ac.PBF_n_func = 3;
                    ac.PBF_n_var = 0;
                    ac.PBF_n_y = 2;
                    ac.PBF_l_func = 5;
                    ac.PBF_l_var = 0;
                    ac.PBF_l_y = 3;
                    ac.PBF_alg_type = 2;
                    ac.PBF_RPC_enabled = false;
                    ac.PBF_safe_skip_enabled = false;
                    ac.PBF_lazy_skip_start_value = 1.0;
                    ac.PBF_lazy_skip_add_value = 1.0;
                    ac.PBF_low_depth = low_d;
                    ac.PBF_high_depth = high_d;
                    ac.PBF_sub_bound_percentage = 0.0;
                    ac.finalize();
                    configs.push_back(ac);
                }
            }
        }

        // PBF: test different values for sub percentage
        std::vector<double> PBF_sub_percentages = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
        for (auto sub_p: PBF_sub_percentages) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.SUB_enabled = true;
            ac.PBF_enabled = true;
            ac.PBF_n_func = 3;
            ac.PBF_n_var = 0;
            ac.PBF_n_y = 2;
            ac.PBF_l_func = 5;
            ac.PBF_l_var = 0;
            ac.PBF_l_y = 3;
            ac.PBF_alg_type = 2;
            ac.PBF_RPC_enabled = false;
            ac.PBF_safe_skip_enabled = false;
            ac.PBF_lazy_skip_start_value = 1.0;
            ac.PBF_lazy_skip_add_value = 1.0;
            ac.PBF_low_depth = 0;
            ac.PBF_high_depth = std::numeric_limits<size_t>::max();
            ac.PBF_sub_bound_percentage = sub_p;
            ac.finalize();
            configs.push_back(ac);
        }

        // Lazy Evaluation: test mode 1
        std::vector<double> LE_score_y_values = {0.5, 1.0, 2.0};
        std::vector<size_t> LE_vars = {1, 2, 3, 4};
        std::vector<double> LE_rank_y_values = {0.25, 0.5, 0.75, 1.0, 2.0, 5.0};
        std::vector<size_t> LE_modes = {1, 2};
        for (auto score_y: LE_score_y_values) {
            for (auto rank_var: LE_vars) {
                for (auto rank_y: LE_rank_y_values) {
                    for (auto mode: LE_modes) {
                        AlgorithmConfiguration ac = AlgorithmConfiguration();
                        ac.LE_mode = mode;
                        ac.LE_y_score_value = score_y;
                        ac.LE_rank_var = rank_var;
                        ac.LE_y_rank_value = rank_y;
                        ac.finalize();
                        configs.push_back(ac);
                    }
                }
            }
        }

        // Recursive heuristic
        std::vector<size_t> REC_l_func = {1, 2};
        std::vector<size_t> REC_l_var = {1, 2, 3, 4};
        std::vector<double> REC_l_y = {0.5, 1.0, 2.0};
        for (auto l_func: REC_l_func) {
            for (auto l_var: REC_l_var) {
                for (auto l_y: REC_l_y) {
                    AlgorithmConfiguration ac = AlgorithmConfiguration();
                    ac.REC_enabled = true;
                    ac.REC_l_func = l_func;
                    ac.REC_l_var = l_var;
                    ac.REC_l_y = l_y;
                    ac.finalize();
                    configs.push_back(ac);
                }
            }
        }

        // Recursive heuristic
        std::vector<double> REC_l = {3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        for (auto l: REC_l) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.REC_enabled = true;
            ac.REC_l_func = 3;
            ac.REC_l_var = 0;
            ac.REC_l_y = l;
            ac.finalize();
            configs.push_back(ac);
        }

        // REC: test different values for rpc
        std::vector<bool> REC_rpc_enabled = {true};
        for (auto rpc_enabled: REC_rpc_enabled) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.REC_enabled = true;
            ac.REC_l_func = 3;
            ac.REC_l_var = 0;
            ac.REC_l_y = 10;
            ac.REC_max_l = 10;
            ac.REC_RPC_enabled = rpc_enabled;
            ac.REC_safe_skip_enabled = false;
            ac.REC_lazy_skip_start_value = 1.0;
            ac.REC_lazy_skip_add_value = 1.0;
            ac.REC_low_depth = 0;
            ac.REC_high_depth = std::numeric_limits<size_t>::max();
            ac.REC_sub_bound_percentage = 0.0;
            ac.finalize();
            configs.push_back(ac);
        }

        // REC: test different values for safe skip
        std::vector<bool> REC_safe_skip_enabled = {true};
        for (auto safe_skip: REC_safe_skip_enabled) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.REC_enabled = true;
            ac.REC_l_func = 3;
            ac.REC_l_var = 0;
            ac.REC_l_y = 10;
            ac.REC_max_l = 10;
            ac.REC_RPC_enabled = false;
            ac.REC_safe_skip_enabled = safe_skip;
            ac.REC_lazy_skip_start_value = 1.0;
            ac.REC_lazy_skip_add_value = 1.0;
            ac.REC_low_depth = 0;
            ac.REC_high_depth = std::numeric_limits<size_t>::max();
            ac.REC_sub_bound_percentage = 0.0;
            ac.finalize();
            configs.push_back(ac);
        }

        // REC: test different values for lazy skip
        std::vector<double> REC_lazy_skip_start_values = {0.0, 0.25, 0.5, 0.75, 1.0};
        std::vector<double> REC_lazy_skip_add_values = {0.1, 0.3, 0.5, 0.7, 0.9};
        for (auto add_val: REC_lazy_skip_add_values) {
            for (auto start_val: REC_lazy_skip_start_values) {
                AlgorithmConfiguration ac = AlgorithmConfiguration();
                ac.REC_enabled = true;
                ac.REC_l_func = 3;
                ac.REC_l_var = 0;
                ac.REC_l_y = 10;
                ac.REC_max_l = 10;
                ac.REC_RPC_enabled = false;
                ac.REC_safe_skip_enabled = false;
                ac.REC_lazy_skip_start_value = start_val;
                ac.REC_lazy_skip_add_value = add_val;
                ac.REC_low_depth = 0;
                ac.REC_high_depth = std::numeric_limits<size_t>::max();
                ac.REC_sub_bound_percentage = 0.0;
                ac.finalize();
                configs.push_back(ac);
            }
        }

        // REC: test different values for depths
        std::vector<size_t> REC_low_depths = {0, 1, 2, 3, 4};
        std::vector<size_t> REC_high_depths = {0, 4, 8, 12, 16, 20};
        for (auto low_d: REC_low_depths) {
            for (auto high_d: REC_high_depths) {
                if (low_d <= high_d) {
                    AlgorithmConfiguration ac = AlgorithmConfiguration();
                    ac.REC_enabled = true;
                    ac.REC_l_func = 3;
                    ac.REC_l_var = 0;
                    ac.REC_l_y = 10;
                    ac.REC_max_l = 10;
                    ac.REC_RPC_enabled = false;
                    ac.REC_safe_skip_enabled = false;
                    ac.REC_lazy_skip_start_value = 1.0;
                    ac.REC_lazy_skip_add_value = 1.0;
                    ac.REC_low_depth = low_d;
                    ac.REC_high_depth = high_d;
                    ac.REC_sub_bound_percentage = 0.0;
                    ac.finalize();
                    configs.push_back(ac);
                }
            }
        }

        // REC: test different values for sub percentage
        std::vector<double> REC_sub_percentages = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
        for (auto sub_p: REC_sub_percentages) {
            AlgorithmConfiguration ac = AlgorithmConfiguration();
            ac.SUB_enabled = true;
            ac.REC_enabled = true;
            ac.REC_l_func = 3;
            ac.REC_l_var = 0;
            ac.REC_l_y = 10;
            ac.REC_max_l = 10;
            ac.REC_RPC_enabled = false;
            ac.REC_safe_skip_enabled = false;
            ac.REC_lazy_skip_start_value = 1.0;
            ac.REC_lazy_skip_add_value = 1.0;
            ac.REC_low_depth = 0;
            ac.REC_high_depth = std::numeric_limits<size_t>::max();
            ac.REC_sub_bound_percentage = sub_p;
            ac.finalize();
            configs.push_back(ac);
        }

        return configs;
    }

    AlgorithmConfiguration get_fast_algorithm_configuration() {
        AlgorithmConfiguration ac = AlgorithmConfiguration();
        ac.SUB_enabled = true;
        ac.RPC_enabled = true;
        ac.LE_mode = 1;
        ac.LE_y_score_value = 1;
        ac.finalize();
        return ac;
    }

}
