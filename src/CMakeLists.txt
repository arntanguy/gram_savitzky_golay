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

add_library(gram_savitzky_golay SHARED
  gram_savitzky_golay.cpp
  spatial_filters.cpp
)

if(TARGET Eigen3::Eigen)
  target_link_libraries(gram_savitzky_golay Eigen3::Eigen)
else()
  target_include_directories(gram_savitzky_golay
    PUBLIC
      ${EIGEN3_INCLUDE_DIRS})
endif()

install(TARGETS gram_savitzky_golay
  EXPORT gram_savitzky_golayTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(EXPORT gram_savitzky_golayTargets
  FILE gram_savitzky_golayTargets.cmake
  NAMESPACE gram_savitzky_golay::
  DESTINATION lib/cmake/gram_savitzky_golay
)

install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
