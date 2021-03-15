# Copyright 2017-2018 CNRS-UM LIRMM
# Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL

macro(add_test_simple NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} PUBLIC Boost::unit_test_framework Boost::disable_autolinking gram_savitzky_golay) 
  target_compile_definitions(${NAME} PRIVATE -DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)
  add_test(${NAME} ${NAME})
  # Adding a project configuration file (for MSVC only)
  generate_msvc_dot_user_file(${NAME})
endmacro()

add_test_simple(test_gram_savitzky_golay)
add_test_simple(test_spatial_filters)

if(TARGET gram_savitzky_golay_sva)
  message(STATUS "Building tests for SpaceVecAlg support")
  macro(add_sva_test NAME)
    add_test_simple(${NAME} ${ARGS})
    target_link_libraries(${NAME} PUBLIC gram_savitzky_golay_sva)
  endmacro()

  add_sva_test(test_sva_filters)
endif()
