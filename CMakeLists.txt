cmake_minimum_required(VERSION 2.8)

include(cmake/base.cmake)
include(cmake/eigen.cmake)
include(cmake/boost.cmake)

set(PROJECT_NAME gram_savitzky_golay)
set(PROJECT_DESCRIPTION "Implementation of Savitzky-Golay Filtering using Peter A. Gory Gram Polynomials method")
set(PROJECT_URL "")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

setup_project()
search_for_eigen()
enable_testing()

find_package(Boost 1.54 COMPONENTS unit_test_framework REQUIRED)

# Enable C++11
add_definitions(-std=c++11 -O3)

include_directories(include)
add_subdirectory(src)

setup_project_finalize()

