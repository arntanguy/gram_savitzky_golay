add_library(gram_savitzky_golay SHARED
  gram_savitzky_golay.cpp
  )

install(TARGETS gram_savitzky_golay
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)


#Add compile target
add_executable(test_gram_savitzky_golay test/test_gram_savitzky_golay.cpp)
target_link_libraries(test_gram_savitzky_golay gram_savitzky_golay)

#link to Boost libraries AND your targets and dependencies
target_link_libraries(test_gram_savitzky_golay ${Boost_LIBRARIES} gram_savitzky_golay)

#Finally add it to test execution - 
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME test_gram_savitzky_golay
              WORKING_DIRECTORY ${CMAKE_CURRENT_DIRECTORY}
              COMMAND ${CMAKE_CURRENT_DIRECTORY}/test_gram_savitzky_golay )
