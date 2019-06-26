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


#Add compile target
add_executable(test_gram_savitzky_golay test_gram_savitzky_golay.cpp)
#link to Boost libraries AND your targets and dependencies
target_link_libraries(test_gram_savitzky_golay Boost::unit_test_framework gram_savitzky_golay)

#Finally add it to test execution - 
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME test_gram_savitzky_golay
              WORKING_DIRECTORY ${CMAKE_CURRENT_DIRECTORY}
              COMMAND ${CMAKE_CURRENT_DIRECTORY}/test_gram_savitzky_golay )



#Add compile target
add_executable(test_spatial_filters test_spatial_filters.cpp)
#link to Boost libraries AND your targets and dependencies
target_link_libraries(test_spatial_filters Boost::unit_test_framework gram_savitzky_golay)

#Finally add it to test execution - 
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME test_spatial_filters
              WORKING_DIRECTORY ${CMAKE_CURRENT_DIRECTORY}
              COMMAND ${CMAKE_CURRENT_DIRECTORY}/test_spatial_filters )

