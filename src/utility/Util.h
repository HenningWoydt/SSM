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

#ifndef SUBSETOPTIMIZATION_UTIL_H
#define SUBSETOPTIMIZATION_UTIL_H

#include <sys/stat.h>
#include <cstddef>
#include <vector>
#include <iostream>
#include <map>
#include <bitset>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>
#include <cmath>
#include <random>
#include <filesystem>
#include <bit>
#include <bitset>
#include <cstdint>
#include <cstring>

namespace SMSM {

#define ASSERT_ENABLED 0
#if ASSERT_ENABLED
#define ASSERT(condition) if(!(condition)){std::cout << "Error in file " __FILE__ << " at line " << __LINE__ << std::endl;abort();} /*if(!(condition)) {__builtin_unreachable();}*/ ((void)0)
#else
#define ASSERT(condition) /*if(!(condition)) {__builtin_unreachable();}*/ ((void)0)
#endif

    /**
     * Pads each string in the vector to the same size.
     *
     * @param v The vector holding the strings.
     * @param front Whether to pad at the font or the back.
     */
    void pad(std::vector<std::string> &v,
             bool front = true);

    /**
     * Print a portion of the elements in a vector.
     *
     * This function prints the first `size` elements of the input vector `v` to the standard
     * output. The elements are enclosed in square brackets and separated by commas. If the vector
     * is empty or the specified size is zero, "[]" is printed. A newline character is also added
     * at the end of the output.
     *
     * @tparam T Type of the elements.
     * @param v The vector whose elements are to be printed.
     * @param size The number of elements from the beginning of the vector to be printed.
     */
    template<typename T>
    void print(const std::vector<T> &v,
               size_t size) {
        if (size == 0) {
            std::cout << "[]\n";
            return;
        }

        std::cout << "[";
        for (size_t i = 0; i < size - 1; ++i) {
            std::cout << +v[i] << ", ";
        }
        std::cout << +v[size - 1] << "] " << "\n";
    }

    /**
     * Print the elements of an array.
     *
     * This function prints the elements of the input array `arr` to the standard output.
     * The elements are enclosed in square brackets and separated by commas. If the array
     * is empty, "[]" is printed. A newline character is also added at the end of the output.
     *
     * @tparam T Type of the elements.
     * @param arr The array whose elements are to be printed.
     * @param size The number of elements in the array to be printed.
     */
    template<typename T>
    void print(const T *arr,
               const size_t size) {
        if (size == 0) {
            std::cout << "[]\n";
            return;
        }

        std::cout << "[";
        for (size_t i = 0; i < size - 1; ++i) {
            std::cout << arr[i] << ", ";
        }
        std::cout << arr[size - 1] << "] " << "\n";
    }

    /**
     * Print the elements of a std::map as a JSON-like representation.
     *
     * This function prints the elements of the input std::map `m` as a JSON-like representation
     * to the standard output. It constructs a string with key-value pairs enclosed in curly
     * braces, where keys and values are converted to strings using the `to_JSON` function.
     * Each key-value pair is separated by a colon, and pairs are separated by commas and newline
     * characters. If the map is empty, an empty curly braces "{}" is printed.
     *
     * @tparam T1 Type of the elements.
     * @tparam T2 Type of the elements.
     * @param m The std::map to be printed as a JSON-like representation.
     */
    template<typename T1, typename T2>
    void print(const std::map<T1, T2> &m) {
        std::string content = "{";

        if (m.empty()) {
            content += "}";
            std::cout << content << std::endl;
            return;
        }

        for (auto entry: m) {
            content += to_JSON(entry.first) + " : " + to_JSON(entry.second) + ",\n";
        }
        content.pop_back();
        content.back() = '}';

        std::cout << content << std::endl;
    }

    /**
     * Print the elements of an array as binary representations.
     *
     * This function prints the binary representations of the elements in the input array `arr`
     * to the standard output. Each binary representation is enclosed in square brackets and
     * separated by commas. If the array is empty, "[]" is printed. A newline character is also
     * added at the end of the output.
     *
     * @tparam T Type of the elements.
     * @param arr The array whose elements are to be printed as binary representations.
     * @param size The number of elements in the array to be printed.
     */
    template<typename T>
    void print_as_bits(const T *arr,
                       const size_t size) {
        if (size == 0) {
            std::cout << "[]\n";
            return;
        }

        std::cout << "[";
        for (size_t i = 0; i < size - 1; ++i) {
            std::bitset<sizeof(T) * 8> bits(arr[i]);
            std::string str(bits.to_string());
            std::cout << str << ", ";
        }
        std::bitset<sizeof(T) * 8> bits(arr[size - 1]);
        std::string str(bits.to_string());
        reverse(str.begin(), str.end());
        std::cout << str << "] " << "\n";
    }

    /**
     * Overload the << operator to stream a vector to an ostream.
     *
     * This overload allows a vector `v` to be streamed into an ostream `os`.
     * The elements of the vector are enclosed in square brackets and separated by commas.
     * If the vector is empty, "[]" is streamed. A newline character is also added at the end.
     *
     * @tparam T Type of the elements.
     * @param os The ostream to which the vector is streamed.
     * @param v The vector to be streamed.
     * @return A reference to the ostream `os` after streaming.
     */
    template<typename T>
    std::ostream &operator<<(std::ostream &os,
                             const std::vector<T> &v) {
        if (v.empty()) {
            os << "[]";
            return os;
        }

        os << "[";
        for (size_t i = 0; i < v.size() - 1; ++i) {
            os << v[i] << ", ";
        }
        os << v.back() << "]";
        return os;
    }

    /**
     * Overload the << operator to stream a vector to a basic_stringstream.
     *
     * This overload allows a vector `v` to be streamed into a basic_stringstream `os`.
     * The elements of the vector are enclosed in square brackets and separated by commas.
     * If the vector is empty, "[]" is streamed. A newline character is also added at the end.
     *
     * @tparam T Type of the elements.
     * @param os The basic_stringstream to which the vector is streamed.
     * @param v The vector to be streamed.
     * @return A reference to the basic_stringstream `os` after streaming.
     */
    template<typename T>
    std::basic_stringstream<char> &operator<<(std::basic_stringstream<char> &os,
                                              std::vector<T> &v) {
        if (v.empty()) {
            os << "[]\n";
            return os;
        }

        os << "[";
        for (size_t i = 0; i < v.size() - 1; ++i) {
            os << v[i] << ", ";
        }
        os << v.back() << "]";
        return os;
    }

    /**
     * Convert a vector to a string representation.
     *
     * This function converts the elements of the input vector `v` into a string representation,
     * where elements are enclosed in square brackets and separated by commas. The resulting
     * string is returned.
     *
     * @tparam T Type of the elements.
     * @param v The vector to be converted to a string representation.
     * @return A string representation of the vector, e.g., "[1, 2, 3]".
     */
    template<typename T>
    std::string to_string(const std::vector<T> &v) {
        if (v.empty()) {
            return std::string{"[]"};
        }

        std::string s = "[";
        for (size_t i = 0; i < v.size() - 1; ++i) {
            s += std::to_string(v[i]) + ", ";
        }
        s += std::to_string(v.back()) + "]";
        return s;
    }


    /**
     * Converts a double into a decimal string representation.
     *
     * @param d The value.
     * @param precision The precision.
     * @return The value as a string.
     */
    std::string double_to_string(double d,
                                 int precision);

    /**
     * Check if a vector contains a specified element.
     *
     * This function uses the standard library's `std::find` algorithm to check if the input
     * vector `vec` contains the element `x`. It returns `true` if the element is found and
     * `false` otherwise.
     *
     * @tparam T Type of the elements.
     * @param vec The vector to be checked for the presence of the element.
     * @param x The element to search for within the vector.
     * @return `true` if the element is found in the vector, `false` otherwise.
     */
    template<typename T>
    bool contains(std::vector<T> &vec,
                  T x) {
        return std::find(vec.begin(), vec.end(), x) != vec.end();
    }

    /**
     * Check if a vector contains a specified element.
     *
     * This function uses the standard library's `std::find` algorithm to check if the input
     * vector `vec` contains the element `x`. It returns `true` if the element is found and
     * `false` otherwise.
     *
     * @tparam T Type of the elements.
     * @param vec The vector to be checked for the presence of the element.
     * @param x The element to search for within the vector.
     * @return `true` if the element is found in the vector, `false` otherwise.
     */
    template<typename T>
    bool contains(const std::vector<T> &vec,
                  T x) {
        return std::find(vec.begin(), vec.end(), x) != vec.end();
    }

    /**
     * Checks if a given char is in a string.
     *
     * @param s The string.
     * @param c The char.
     * @return True if the element is present in the string, false else.
     */
    bool contains(std::string &s,
                  char c);

    /**
     * Check if a vector contains a specified element within a specified size.
     *
     * This function checks if the input vector `vec` contains the element `x` within the
     * first `size` elements. It iterates through the vector to determine if the element
     * is present within the specified range.
     *
     * @tparam T Type of the elements.
     * @param vec The vector to be checked for the presence of the element.
     * @param size The number of elements within the vector to consider for the search.
     * @param x The element to search for within the vector.
     * @return `true` if the element is found in the vector, `false` otherwise.
     */
    template<typename T>
    bool contains(std::vector<T> &vec,
                  size_t size,
                  T x) {
        for (size_t i = 0; i < size; ++i) {
            if (vec[i] == x) {
                return true;
            }
        }
        return false;
    }

    /**
     * Checks if a given element in in a set.
     *
     * @param vec The vector.
     * @param x The element.
     * @return True if the element is present in the vector, false else.
     */
    bool contains(std::set<uint32_t> &set,
                  uint32_t x);

    /**
     * Check if an array contains a specified element.
     *
     * This function checks if the input array `arr` of the given `size` contains the
     * element `x`. It iterates through the array to determine if the element is present.
     *
     * @tparam T Type of the elements.
     * @param arr The array to be checked for the presence of the element.
     * @param size The number of elements in the array.
     * @param x The element to search for within the array.
     * @return `true` if the element is found in the array, `false` otherwise.
     */
    template<typename T>
    bool contains(const T *arr,
                  const size_t size,
                  const T x) {
        for (size_t i = 0; i < size; ++i) {
            if (arr[i] == x) {
                return true;
            }
        }
        return false;
    }

    /**
     * Check if a vector contains duplicate elements within a specified size.
     *
     * This function checks if the input vector `vec` contains any duplicate elements
     * within the first `size` elements. It compares each element with all subsequent
     * elements to determine if any duplicates are present.
     *
     * @tparam T Type of the elements.
     * @param vec The vector to be checked for duplicate elements.
     * @param size The number of elements within the vector to consider for duplicates.
     * @return `true` if duplicates are found, `false` otherwise.
     */
    template<typename T>
    bool contains_duplicate(std::vector<T> &vec,
                            size_t size) {
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = i + 1; j < size; ++j) {
                if (vec[i] == vec[j]) {
                    return true;
                }
            }
        }
        return false;
    }


    /**
     * Splits a string at the specified char.
     *
     * @param s The string.
     * @param c The char.
     * @return Vector holding all substrings.
     */
    std::vector<std::string> split(std::string &s,
                                   char c);

    /**
     * Converts a string into an argv, argc representation.
     *
     * @param s The string.
     * @param argc The resulting number of arguments.
     * @return Pointer to argv.
     */
    char **string_to_argv(std::string &s, int &argc);

    /**
     * Frees all memory allocated by argv.
     *
     * @param argc Number of arguments.
     * @param argv The arguments.
     */
    void free_argv(int argc, char **argv);

    /**
     * Prints argv onto the standard output.
     *
     * @param argc Number of arguments.
     * @param argv The arguments.
     */
    void print_argv(int argc, char **argv);

    /**
     * Find the maximum element in a vector.
     *
     * This function finds and returns the maximum element in the input vector `vec` using
     * the standard library's `std::max_element` algorithm.
     *
     * @tparam T Type of the elements.
     * @param vec The vector in which to find the maximum element.
     * @return The maximum element in the vector.
     */
    template<typename T>
    T max(std::vector<T> &vec) {
        return *std::max_element(vec.begin(), vec.end());
    }

    /**
     * Find the minimum element in a vector.
     *
     * This function finds and returns the minimum element in the input vector `vec` using
     * the standard library's `std::min_element` algorithm.
     *
     * @tparam T Type of the elements.
     * @param vec The vector in which to find the minimum element.
     * @return The minimum element in the vector.
     */
    template<typename T>
    T min(const std::vector<T> &vec) {
        return *std::min_element(vec.begin(), vec.end());
    }

    /**
     * Copy the contents of one vector into another.
     *
     * This function takes a source vector `vec` and copies its elements into the destination
     * vector `res`. The destination vector `res` is resized to match the size of the source
     * vector `vec` before copying the elements. After execution, both vectors will contain
     * identical elements.
     *
     * @tparam T Type of the elements.
     * @param res Destination vector to store the copied elements.
     * @param vec Source vector from which elements are copied.
     */
    template<typename T>
    void copy(std::vector<T> &res, const std::vector<T> &vec) {
        res.resize(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            res[i] = vec[i];
        }
    }

    /**
     * Copy the contents of one vector into another.
     *
     * @tparam T Type of the elements.
     * @param res Destination vector.
     * @param vec Source vector.
     * @param size Number of elements to copy.
     */
    template<typename T>
    void copy(std::vector<T> &res, const std::vector<T> &vec, const size_t size) {
        res.resize(size);

        // Use memmove to copy the memory block
        std::memcpy(res.data(), vec.data(), size * sizeof(T));
    }

    /**
     * Copy the contents of one vector into another.
     *
     * @tparam T Type of the elements.
     * @param src Source vector.
     * @param src_offset Offset into the source vector.
     * @param dest Destination vector.
     * @param dest_offset Offset into the destination vector.
     * @param n Number of elements to copy.
     */
    template<typename T>
    void copy(const std::vector<T> &src, const size_t src_offset, std::vector<T> &dest, const size_t dest_offset, const size_t n) {
        // Use memmove to copy the memory block
        std::memcpy(dest.data() + dest_offset, src.data() + src_offset, n * sizeof(T));
    }

    /**
     * Calculate the sum of elements in a vector up to a specified size.
     *
     * This function computes the sum of the elements in the input vector `vec` up to a given
     * size `size`. It iterates through the vector and accumulates the values, returning the
     * result as the summation.
     *
     * @tparam T Type of the elements.
     * @param vec The vector containing elements to be summed.
     * @param size The number of elements to include in the summation.
     * @return The sum of the first `size` elements in the vector.
     */
    template<typename T>
    T sum(const std::vector<T> &vec, size_t size) {
        T s = 0;
        for (size_t i = 0; i < size; ++i) {
            s += vec[i];
        }
        return s;
    }

    /**
     * Calculate the sum of elements in an array up to a specified size.
     *
     * This function computes the sum of the elements in the input array `arr` up to a given
     * size `size`. It iterates through the array and accumulates the values, returning the
     * result as the summation.
     *
     * @tparam T Type of the elements.
     * @param arr The array containing elements to be summed.
     * @param size The number of elements to include in the summation.
     * @return The sum of the first `size` elements in the array.
     */
    template<typename T>
    T sum(const T *arr, size_t size) {
        T s = 0;
        for (size_t i = 0; i < size; ++i) {
            s += arr[i];
        }
        return s;
    }

    /**
     * Returns the factorial of n.
     *
     * @param n The non negative integer.
     * @return n!
     */
    size_t fact(size_t n);

    /**
     * Return the number of combinations when choosing k elements from a pool of n
     * elements with no duplicates.
     *
     * @param n Non negative integer.
     * @param k Non negative integer.
     * @return n over k.
     */
    size_t nCk(size_t n,
               size_t k);

    /**
     * Calculate the binomial coefficient "n choose k."
     *
     *
     * @param n The total number of items in the set.
     * @param k The number of items to choose from the set.
     * @return The binomial coefficient "n choose k."
     */
    size_t n_choose_k(size_t n, size_t k);

    /**
     * Calculate the power of 2 raised to the given exponent.
     *
     * @param n The exponent for the power of 2.
     * @return The result of 2 raised to the power of 'n'.
     */
    size_t pow_2(size_t n);

    /**
     * Sort the vector and additionally applies the moves on the other vectors.
     *
     * @tparam T1 Type of the vector to sort.
     * @tparam T2 Type of the second vector
     * @tparam T3 Type of the third vector.
     * @param vec_to_sort The vector, that will be sorted.
     * @param vec2 The second vector.
     * @param vec3 The third vector.
     * @param size The size of the vectors.
     * @param ascending If it should be sorted ascending.
     */
    template<class T1, class T2, class T3>
    void co_sort(std::vector<T1> &vec_to_sort,
                 std::vector<T2> &vec2,
                 std::vector<T3> &vec3,
                 size_t size,
                 bool ascending) {
        bool sorted = false;
        for (size_t i = 0; !sorted && i < size; ++i) {
            sorted = true;
            for (size_t j = 0; j < size - i - 1; ++j) {
                bool descending_unsorted = !ascending && (vec_to_sort[j] < vec_to_sort[j + 1]);
                bool ascending_unsorted = ascending && (vec_to_sort[j] > vec_to_sort[j + 1]);

                if (descending_unsorted || ascending_unsorted) {
                    std::swap(vec_to_sort[j], vec_to_sort[j + 1]);
                    std::swap(vec2[j], vec2[j + 1]);
                    std::swap(vec3[j], vec3[j + 1]);
                    sorted = false;
                }
            }
        }
    }

    /**
     * Sort the vector and additionally applies the moves on the second vector.
     * The second vector holds twice as much elements as the first vector and two
     * adjacent entries are considered one pair.
     *
     * @tparam T1 Type of the vector to sort.
     * @tparam T2 Type of the second vector
     * @param vec_to_sort The vector, that will be sorted.
     * @param vec2 The second vector.
     * @param size The size of the first vector.
     * @param ascending If it should be sorted ascending.
     */
    template<class T1, class T2>
    void co_sort(std::vector<T1> &vec_to_sort,
                 std::vector<T2> &vec2,
                 size_t size,
                 bool ascending) {
        bool sorted = false;
        for (size_t i = 0; !sorted && i < size; ++i) {
            sorted = true;
            for (size_t j = 0; j < size - i - 1; ++j) {
                bool descending_unsorted = !ascending && (vec_to_sort[j] < vec_to_sort[j + 1]);
                bool ascending_unsorted = ascending && (vec_to_sort[j] > vec_to_sort[j + 1]);

                if (descending_unsorted || ascending_unsorted) {
                    std::swap(vec_to_sort[j], vec_to_sort[j + 1]);
                    std::swap(vec2[j * 2], vec2[(j + 1) * 2]);
                    std::swap(vec2[j * 2 + 1], vec2[(j + 1) * 2 + 1]);
                    sorted = false;
                }
            }
        }
    }

    /**
     * Sort the vector and additionally applies the moves on the other vectors.
     *
     * @param vec_to_sort The vector, that will be sorted.
     * @param vec2 The second vector.
     * @param vec3 The third vector.
     * @param size The size of the vectors.
     * @param ascending If it should be sorted ascending.
     */
    template<>
    void co_sort<double, uint32_t, bool>(std::vector<double> &vec_to_sort,
                                         std::vector<uint32_t> &vec2,
                                         std::vector<bool> &vec3,
                                         size_t size,
                                         bool ascending);

    /**
     * Round up a value to the nearest multiple of another value.
     *
     * This function rounds up the value `n` to the nearest multiple of the value `m`. If `n`
     * i s already a multiple of `m`, it returns `n` unchanged. Otherwise, it returns the smallest
     * multiple of `m` that is greater than or equal to `n`.
     *
     * https://stackoverflow.com/questions/3407012/rounding-up-to-the-nearest-multiple-of-a-number
     *
     * @tparam T Type of the elements.
     * @param n The value to be rounded up.
     * @param m The value to which `n` is rounded up to become a multiple.
     * @return The rounded-up value.
     */
    template<typename T>
    T round_up(T n,
               T m) {
        T r = n % m;
        if (r == 0) {
            return n;
        }

        return n + m - r;
    }

    size_t ceil(size_t x, size_t y);

    /**
     * Checks if all the values are close to each other.
     *
     * @param vec Vector of double values.
     * @return True if all are close to each other, false else.
     */
    bool all_close(std::vector<double> &vec);

    /**
     * Get a random integer.
     *
     * @return Random integer.
     */
    size_t random_int();

    /**
     * Gets all files in a specific directory.
     *
     * @param directory_path Path to the directory.
     * @return Vector holding all paths.
     */
    std::vector<std::string> get_directory_files(std::string &directory_path);

    /**
     * Compute a compile-time hash value for a constant character string.
     *
     * This function calculates a hash value for a constant character string `input` at compile time.
     * The hash value is computed using a simple algorithm that iterates over the characters in the
     * string, applying a hash formula. The resulting hash value can be used as a compile-time constant.
     *
     * @param input The constant character string for which to calculate the hash.
     * @return The computed hash value for the input string.
     */
    size_t constexpr const_hash(char const *input) {
        size_t i = 0;
        size_t hash = 0;
        while (input[i] != '\0') {
            hash += static_cast<unsigned int>(input[i]) + 33 * hash;
            i++;
        }
        hash += 5381;

        return hash;
    }

    /**
     * Compute the cumulative sum of elements in a subvector.
     *
     * This function calculates the cumulative sum of elements in a subvector of the input
     * vector `vec`. It computes the cumulative sum starting from the `offset` position and
     * includes `n` elements in the result vector `res`. The result vector `res` is resized
     * to have `n + 1` elements, where the first element is always set to zero.
     *
     * @tparam T Type of the elements.
     * @param vec The input vector containing elements to compute the cumulative sum from.
     * @param res The resulting vector to store the cumulative sums.
     * @param offset The starting index within the input vector `vec`.
     * @param n The number of elements from the starting index to include in the cumulative sum.
     */
    template<typename T>
    void cumulative_sum(const std::vector<T> &vec,
                        std::vector<T> &res,
                        size_t offset,
                        size_t n) {
        res.resize(n + 1);

        T sum = 0.0;
        res[0] = sum;
        for (size_t i = 0; i < n; ++i) {
            sum += vec[offset + i];
            res[i + 1] = sum;
        }
    }


    /**
     * Calculate the sum of elements from a cumulative sum array.
     *
     * This function calculates the sum of elements in a range of a cumulative sum array `csum`.
     * It subtracts the cumulative sum at the `start_idx` position from the cumulative sum at
     * the position `start_idx + n`. The provided cumulative sum array must have a sufficient
     * size to access the specified indices.
     *
     * @tparam T Type of the elements.
     * @param csum The cumulative sum array from which to calculate the sum.
     * @param start_idx The starting index within the cumulative sum array.
     * @param n The number of elements to include in the sum.
     * @return The sum of elements within the specified range of the cumulative sum array.
     */
    template<typename T>
    inline T get_sum_from_csum(const std::vector<T> &csum,
                               size_t start_idx,
                               size_t n) {
        ASSERT(start_idx + n < csum.size());
        return csum[start_idx + n] - csum[start_idx];
    }

    /**
     * Checks if a file exists.
     *
     * @param file_path Path to the file.
     * @return True if the file exists, false else.
     */
    bool file_exists(const std::string &file_path);

    /**
     * Checks if two numbers are almost equal.
     *
     * @param a Number a.
     * @param b Number b.
     * @param epsilon The error tolerance.
     * @return True if they are almost equal, false else.
     */
    bool double_eq(double a, double b, double epsilon);

    /**
     * Returns the current time point.
     *
     * @return The time point.
     */
    std::chrono::steady_clock::time_point get_time_point();

    /**
     * Get the elapsed seconds between the both time points.
     *
     * @param sp Start point.
     * @param ep End point.
     * @return Elapsed time in seconds
     */
    double get_elapsed_seconds(std::chrono::steady_clock::time_point sp,
                               std::chrono::steady_clock::time_point ep);

    /**
     * Iterates over all possible subsets of {0, ..., n - 1}.
     *
     * @param set The set.
     * @param set_size Size of the set.
     * @param n Greatest possible element.
     * @return True if the subset is valid.
     */
    bool next_subset(std::vector<size_t> &set, size_t set_size, size_t n);

    /**
     * Iterates over all possible subsets of {0, ..., 63} with size k.
     *
     * @param bits The set.
     * @param k Size of the set.
     * @return True if the subset is valid.
     */
    size_t next_subset(size_t bits, int k);

    /**
     * Remove all occurrences of a specific element from the end of a vector.
     *
     * This function removes all occurrences of the specified element `x` from the end of
     * the input vector `vec`. It repeatedly removes elements from the end until a non-matching
     * element is encountered or the vector becomes empty.
     *
     * @tparam T Type of the elements.
     * @param vec The vector from which to remove elements.
     * @param x The element to be removed from the end of the vector.
     */
    template<typename T>
    void pop_all_x(std::vector<T> &vec, T x) {
        while (!vec.empty() && vec.back() == x) {
            vec.pop_back();
        }
    }

    /**
     * Division operation.
     *
     * @param res Result Vector.
     * @param v1 Numerator.
     * @param v2 Denominator.
     */
    void div(std::vector<double> &res, const std::vector<double> &v1, const std::vector<double> &v2);

    /**
     * Division operation.
     *
     * @param res Result Vector.
     * @param v1 Numerator.
     * @param v2 Denominator.
     */
    void div(std::vector<double> &res, const std::vector<size_t> &v1, const std::vector<size_t> &v2);

    /**
     * Division operation.
     *
     * @param res Result Vector.
     * @param v1 Numerator.
     * @param v2 Denominator.
     */
    void div(std::vector<double> &res, const std::vector<double> &v1, const std::vector<size_t> &v2);

    /**
     * Division operation.
     *
     * @param res Result Vector.
     * @param v1 Numerator.
     * @param v2 Denominator.
     */
    void div(std::vector<double> &res, const std::vector<size_t> &v1, const std::vector<double> &v2);

    /***** Array Functions *****/

    /**
     * Elementwise minimum.
     *
     * @tparam T Type of the elements.
     * @param res Resulting vector.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @param n Number of elements.
     */
    template<typename T>
    void min(T *__restrict__ res, const T *__restrict__ v1, const T *__restrict__ v2, const size_t n) {
        for (size_t i = 0; i < n; ++i) {
            res[i] = std::min(v1[i], v2[i]);
        }
    }

    /**
     * Elementwise minimum in-place.
     *
     * @tparam T Type of the elements.
     * @param res Resulting vector and vector 1.
     * @param v2 Vector 2.
     * @param n Number of elements.
     */
    template<typename T>
    void min_in_place(T *__restrict__ res, const T *__restrict__ v2, size_t n) {
        for (size_t i = 0; i < n; ++i) {
            res[i] = std::min(res[i], v2[i]);
        }
    }

    /**
     * Sum of the elementwise minimum.
     *
     * @tparam T Type of the elements.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @param n Number of elements.
     * @return The sum.
     */
    template<typename T>
    T sum_of_min(const T *__restrict__ v1, const T *__restrict__ v2, size_t n) {
        T sum = 0;
        for (size_t i = 0; i < n; ++i) {
            sum += std::min(v1[i], v2[i]);
        }
        return sum;
    }

    /**
     * Sum of the elementwise minimum.
     *
     * @tparam T Type of the elements.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @param v3 Vector 3.
     * @param n Number of elements.
     * @return The sum.
     */
    template<typename T>
    T sum_of_min(const T *__restrict__ v1, const T *__restrict__ v2, const T *__restrict__ v3, size_t n) {
        T sum = 0;
        for (size_t i = 0; i < n; ++i) {
            sum += std::min(std::min(v1[i], v2[i]), v3[i]);
        }
        return sum;
    }

    /**
     * Overwrites the first vector with content of the second vector.
     *
     * @tparam T Type of the elements.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @param n Number of elements.
     */
    template<typename T>
    void overwrite(T *__restrict__ v1, const T *__restrict__ v2, size_t n) {
        for (size_t i = 0; i < n; ++i) {
            v1[i] = v2[i];
        }
    }

    /***** Vector Functions *****/

    /**
     * Sums the elements in the vector.
     *
     * @tparam T Type of the elements.
     * @param vec The vector.
     * @return Sum of its elements.
     */
    template<typename T>
    T sum(const std::vector<T> &vec) {
        T s = 0;
        for (size_t i = 0; i < vec.size(); ++i) {
            s += vec[i];
        }
        return s;
    }

    /**
     * Elementwise minimum.
     *
     * @tparam T Type of the elements.
     * @param res Resulting vector.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     */
    template<typename T>
    void min(std::vector<T> &res, const std::vector<T> &v1, const std::vector<T> &v2) {
        min(res.data(), v1.data(), v2.data(), res.size());
    }

    /**
     * Elementwise minimum in-place.
     *
     * @tparam T Type of the elements.
     * @param res Resulting vector and vector 1.
     * @param v2 Vector 2.
     */
    template<typename T>
    void min_in_place(std::vector<T> &res, const std::vector<T> &v2) {
        min_in_place(res.data(), v2.data(), res.size());
    }

    /**
     * Sum of elementwise minimum.
     *
     * @tparam T Type of the elements.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @return The sum.
     */
    template<typename T>
    T sum_of_min(const std::vector<T> &v1, const std::vector<T> &v2) {
        return sum_of_min(v1.data(), v2.data(), v1.size());
    }

    /**
     * Sum of elementwise minimum.
     *
     * @tparam T Type of the elements.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @param v3 Vector 3.
     * @return The sum.
     */
    template<typename T>
    T sum_of_min(const std::vector<T> &v1, const std::vector<T> &v2, const std::vector<T> &v3) {
        return sum_of_min(v1.data(), v2.data(), v3.data(), v1.size());
    }

    /**
     * Overwrites vector 1 with vector 2.
     *
     * @tparam T Type of the elements.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     */
    template<typename T>
    void overwrite(std::vector<T> &v1, const std::vector<T> &v2) {
        overwrite(v1.data(), v2.data(), v1.size());
    }

}

#endif //SUBSETOPTIMIZATION_UTIL_H
