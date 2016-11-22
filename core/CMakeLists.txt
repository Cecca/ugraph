## Main build

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/git_info.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp" @ONLY)

set( HEADER_FILES
  git_info.hpp
  prelude.hpp
  termcolor.hpp
  logging.hpp
  rand.hpp
  experiment_reporter.hpp
  io.hpp
  sampler.hpp
  sequential_clustering.hpp
)

set( LIB_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp"
  io.cpp
  logging.cpp
  rand.cpp
  sampler.cpp
  sequential_clustering.cpp
)
add_library(ugraph STATIC ${HEADER_FILES} ${LIB_FILES} )

add_executable(scratch scratch.cpp)
target_link_libraries(scratch
  ugraph
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(ugraph-sc sequential.cpp)
target_link_libraries(ugraph-sc
  ugraph
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY})
