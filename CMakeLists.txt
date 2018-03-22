# Copyright 2017-2018 CNRS-UM LIRMM
# Copyright 2017-2018 Arnaud TANGUY <arnaud.tanguy@lirmm.fr>
#
# This file is part of robcalib.
#
# robcalib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# robcalib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with robcalib.  If not, see <http://www.gnu.org/licenses/>.

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

