# Copyright 2017-2018 CNRS-UM LIRMM
# Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL

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
