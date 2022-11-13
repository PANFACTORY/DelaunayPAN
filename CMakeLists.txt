cmake_minimum_required(VERSION 3.24)

project(delaunaypan_test CXX)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/release-1.12.0.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()

add_executable(DelaunayPANTest test/node_test.cpp)
target_link_libraries(DelaunayPANTest gtest_main)

include(GoogleTest)
gtest_discover_tests(DelaunayPANTest)
