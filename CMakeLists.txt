cmake_minimum_required(VERSION 3.0)
project(reliable_clustering)

set(CMAKE_BUILD_TYPE Release)

# Git info handling
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIBE)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wall -std=c++11")

find_package(BZip2 REQUIRED)

find_package(Boost REQUIRED
  COMPONENTS graph system filesystem date_time iostreams
             program_options unit_test_framework)

include_directories(${Boost_INCLUDE_DIRS})

## Add the `include` directory to the include path, in order to be able to
## include header-only libraries stored there and managed by Peru
include_directories(include)

add_subdirectory(core ${CMAKE_BINARY_DIR}/core)
