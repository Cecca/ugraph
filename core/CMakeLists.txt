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
  connected_components.hpp
  cc_sampler.hpp
  bfs_sampler.hpp
  cluster_vertex.hpp
  guesser.hpp
  mcpc.hpp
  acpc.hpp
  counts_cache.hpp
)

set( LIB_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp"
  io.cpp
  logging.cpp
  rand.cpp
  connected_components.cpp
  cc_sampler.cpp
  bfs_sampler.cpp
  scores.cpp
)
add_library(ugraph STATIC ${HEADER_FILES} ${LIB_FILES} )

add_executable(ugraph-mcpc mcpc-run.cpp)
target_link_libraries(ugraph-mcpc
  ugraph
  ${Boost_LIBRARIES}
  ${BZIP2_LIBRARIES})

add_executable(ugraph-acpc acpc-run.cpp)
target_link_libraries(ugraph-acpc
        ugraph
        ${Boost_LIBRARIES}
        ${BZIP2_LIBRARIES})
      
add_executable(ugraph-gmm gmm.cpp)
target_link_libraries(ugraph-gmm
        ugraph
        ${Boost_LIBRARIES}
        ${BZIP2_LIBRARIES})

add_executable(ugraph-scores score_computer.cpp)
target_link_libraries(ugraph-scores
  ugraph
  ${Boost_LIBRARIES}
  ${BZIP2_LIBRARIES})

add_executable(ugraph-avg-prob avg-prob.cpp)
target_link_libraries(ugraph-avg-prob
  ugraph
  ${Boost_LIBRARIES}
  ${BZIP2_LIBRARIES})
