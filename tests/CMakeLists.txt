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


macro(add_test_simple NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} PUBLIC Boost::unit_test_framework Boost::disable_autolinking gram_savitzky_golay) 
  target_compile_definitions(${NAME} PRIVATE -DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)
  add_test(${NAME} ${NAME})
  # Adding a project configuration file (for MSVC only)
  generate_msvc_dot_user_file(${NAME})
ENDMACRO()

add_test_simple(test_gram_savitzky_golay)
add_test_simple(test_spatial_filters)
