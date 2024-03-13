cmake_minimum_required(VERSION 3.24)
project(SSM)

set(CMAKE_CXX_STANDARD 20)
set(Boost_NO_WARN_NEW_VERSIONS 1)


find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(CMAKE_CXX_FLAGS_RELEASE "-std=c++20 -O3 -DNDEBUG -fopenmp -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -fopenmp -std=c++20 -Wall -Wextra -pedantic")


enable_testing()

include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

file(GLOB Dir1_Sources "3rd_party_tools/blossom5-v2.05.src/*.cpp")
file(GLOB Dir2_Sources "3rd_party_tools/blossom5-v2.05.src/*.h")
file(GLOB Dir3_Sources "3rd_party_tools/blossom5-v2.05.src/GEOM/*.cpp")
file(GLOB Dir4_Sources "3rd_party_tools/blossom5-v2.05.src/GEOM/*.h")
file(GLOB Dir5_Sources "3rd_party_tools/blossom5-v2.05.src/MinCost/*.cpp")
file(GLOB Dir6_Sources "3rd_party_tools/blossom5-v2.05.src/MinCost/*.h")
set(Blossom5_files ${Dir1_Sources} ${Dir2_Sources} ${Dir3_Sources} ${Dir4_Sources} ${Dir5_Sources} ${Dir6_Sources})
list(FILTER Blossom5_files EXCLUDE REGEX "example.cpp")


add_executable(SSM
        main.cpp
        ${Blossom5_files}
        src/algorithms/BFSearch.h
        src/algorithms/CandidateManager.h
        src/algorithms/PBFAlgorithm.h
        src/algorithms/SICache.h
        src/algorithms/TreeSearch.h
        src/algorithms/UB2DAlgorithm.h
        src/algorithms/UB2DStorage.cpp
        src/algorithms/UB2DAlgorithm.h
        src/algorithms/UpperBoundManager.h
        src/structures/DataPoints.h
        src/structures/DataPointsEuclidianDistance.h
        src/structures/Graph.h
        src/structures/GraphNegativeGroupFarness.h
        src/structures/GraphPartialDominatingSet.h
        src/structures/StructureInterface.h
        src/utility/AlgorithmConfiguration.cpp
        src/utility/AlgorithmConfiguration.h
        src/utility/DataPointGeneration.cpp
        src/utility/DataPointGeneration.h
        src/utility/GraphGeneration.cpp
        src/utility/GraphGeneration.h
        src/utility/JSONUtil.cpp
        src/utility/JSONUtil.h
        src/utility/Util.cpp
        src/utility/Util.h
        src/utility/VectorOfVectors.h
        src/structures/FacilityLocation.h)

add_executable(SSM_test
        ${Blossom5_files}
        src/algorithms/BFSearch.h
        src/algorithms/CandidateManager.h
        src/algorithms/PBFAlgorithm.h
        src/algorithms/SICache.h
        src/algorithms/TreeSearch.h
        src/algorithms/UB2DAlgorithm.h
        src/algorithms/UB2DStorage.cpp
        src/algorithms/UB2DAlgorithm.h
        src/algorithms/UpperBoundManager.h
        src/structures/DataPoints.h
        src/structures/DataPointsEuclidianDistance.h
        src/structures/Graph.h
        src/structures/GraphNegativeGroupFarness.h
        src/structures/GraphPartialDominatingSet.h
        src/structures/StructureInterface.h
        src/utility/AlgorithmConfiguration.cpp
        src/utility/AlgorithmConfiguration.h
        src/utility/DataPointGeneration.cpp
        src/utility/DataPointGeneration.h
        src/utility/GraphGeneration.cpp
        src/utility/GraphGeneration.h
        src/utility/JSONUtil.cpp
        src/utility/JSONUtil.h
        src/utility/Util.cpp
        src/utility/Util.h
        src/utility/VectorOfVectors.h
        tests/DataPointsEuclidianDistanceTest.cpp
        tests/GraphNegativeGroupFarnessTest.cpp
        tests/GraphPartialDominatingSetTest.cpp
        tests/MiscTest.cpp
        src/structures/FacilityLocation.h)

if (supported)
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET SSM PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET SSM_test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else ()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif ()

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(SSM OpenMP::OpenMP_CXX)
    target_link_libraries(SSM_test OpenMP::OpenMP_CXX)
endif ()

target_link_libraries(SSM ${Boost_LIBRARIES})
target_link_libraries(SSM_test ${Boost_LIBRARIES})

target_link_libraries(SSM_test GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(SSM_test)
