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

cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME gram_savitzky_golay)
set(PROJECT_DESCRIPTION "Implementation of Savitzky-Golay Filtering using Peter A. Gory Gram Polynomials method")
set(PROJECT_URL "https://github.com/arntanguy/gram_savitzky_golay")
set(PROJECT_VERSION 1.0)
set(PROJECT_DEBUG_POSTFIX "_d")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(PROJECT_USE_KEYWORD_LINK_LIBRARIES TRUE)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(DOXYGEN_USE_MATHJAX "YES")
set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

option(BUILD_TESTS "Build unit tests" ON)

include(cmake/base.cmake)
include(cmake/msvc-specific.cmake)

project(${PROJECT_NAME} LANGUAGES CXX)

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif()


add_project_dependency(Eigen3 REQUIRED PKG_CONFIG_REQUIRES eigen3)
add_project_dependency(Boost REQUIRED)
add_subdirectory(src)

if(BUILD_TESTS)
  find_package(Boost COMPONENTS unit_test_framework REQUIRED)
  enable_testing()
  add_subdirectory(tests)
endif()
