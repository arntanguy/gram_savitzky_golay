# Copyright 2017-2018 CNRS-UM LIRMM
# Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL

set(SRC
gram_savitzky_golay.cpp
filters/rotation.cpp
)

set(HEADERS
../include/gram_savitzky_golay/gram_savitzky_golay.h
../include/gram_savitzky_golay/spatial_filters.h
)

add_library(gram_savitzky_golay SHARED ${SRC} ${HEADERS})
set_target_properties(gram_savitzky_golay PROPERTIES COMPILE_FLAGS "-DGRAM_SAVITZKY_GOLAY_EXPORTS")
target_include_directories(gram_savitzky_golay PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
if(TARGET Eigen3::Eigen)
  target_link_libraries(gram_savitzky_golay PUBLIC Eigen3::Eigen)
else()
  target_include_directories(gram_savitzky_golay PUBLIC ${EIGEN3_INCLUDE_DIRS})
endif()
if(TARGET Boost::boost)
  target_link_libraries(gram_savitzky_golay PUBLIC Boost::boost)
else()
  target_include_directories(gram_savitzky_golay PUBLIC ${BOOST_INCLUDE_DIRS})
  target_link_libraries(gram_savitzky_golay PUBLIC ${BOOST_LIBRARIES})
endif()

# pkg_config_append_libs(gram_savitzky_golay)
set_target_properties(gram_savitzky_golay PROPERTIES SOVERSION ${PROJECT_VERSION})
install(TARGETS gram_savitzky_golay
  EXPORT "${TARGETS_EXPORT_NAME}"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)

if(TARGET SpaceVecAlg::SpaceVecAlg)
  add_library(gram_savitzky_golay_sva INTERFACE)
  target_link_libraries(gram_savitzky_golay_sva INTERFACE gram_savitzky_golay SpaceVecAlg::SpaceVecAlg)

  install(TARGETS gram_savitzky_golay_sva
    EXPORT "${TARGETS_EXPORT_NAME}")
endif()

