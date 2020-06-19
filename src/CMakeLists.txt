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


set(SRC
gram_savitzky_golay.cpp
spatial_filters.cpp
)

set(HEADERS
../include/gram_savitzky_golay/gram_savitzky_golay.h
../include/gram_savitzky_golay/spatial_filters.h
)

add_library(gram_savitzky_golay SHARED ${SRC} ${HEADERS})
target_include_directories(gram_savitzky_golay PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
if(TARGET Eigen3::Eigen)
  target_link_libraries(gram_savitzky_golay PUBLIC Eigen3::Eigen)
else()
  target_include_directories(gram_savitzky_golay PUBLIC ${EIGEN3_INCLUDE_DIRS})
endif()
set_target_properties(gram_savitzky_golay PROPERTIES SOVERSION ${PROJECT_VERSION})

install(TARGETS gram_savitzky_golay
  EXPORT "${TARGETS_EXPORT_NAME}"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
