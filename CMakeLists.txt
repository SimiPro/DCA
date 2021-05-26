cmake_minimum_required(VERSION 3.16)

include(FetchContent)
include(CMakeDependentOption)

macro(fetch what)
    FetchContent_GetProperties("${what}")
    if(NOT ${${what}_POPULATED})
        message(STATUS "fetching ${what} ...")
        FetchContent_Populate(${what})
    endif()
    mark_as_advanced(${${what}_SOURCE_DIR})
endmacro()

project("DCA" VERSION 0.1
             DESCRIPTION "Differentiable Collision Avoidance"
             HOMEPAGE_URL "https://github.com/SimiPro/DCA")

set(CMAKE_CXX_STANDARD 20)
file(GLOB DCA_SOURCES "src/*cpp")
file(GLOB_RECURSE DCA_HEADERS "include/DCA/*.h")

option(DCA_BUILD_DEV_TOOLS "Build DEV tools (like FD and example app)" OFF)

# These options are only available if we build dev tools (finite difference)
cmake_dependent_option(DCA_BUILD_EXAMPLES "Build examples" ON DCA_BUILD_DEV_TOOLS OFF)

add_library(${PROJECT_NAME} ${DCA_SOURCES} ${DCA_HEADERS})
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

if(DCA_BUILD_DEV_TOOLS)
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DDCA_DEV_TOOLS=1)
else(DCA_BUILD_DEV_TOOLS)
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DDCA_DEV_TOOLS=0)
endif(DCA_BUILD_DEV_TOOLS)

target_include_directories(
  ${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

find_package (Eigen3 3.3 NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  # Remove entry that has not been found
  unset(Eigen3_DIR CACHE)

  # Now download eigen3
  # download as zip, so it's a bit smaller than downloading the whole repo
  FetchContent_Declare(
      eigen #
      URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip #
      URL_HASH MD5=888aab45512cc0c734b3e8f60280daba #
  )
  fetch(eigen)
  add_library(eigen INTERFACE)
  add_library (Eigen3::Eigen ALIAS eigen)
  target_include_directories(eigen INTERFACE ${eigen_SOURCE_DIR})
endif(NOT TARGET Eigen3::Eigen)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

if(${DCA_BUILD_EXAMPLES})
  add_subdirectory(examples)
endif(${DCA_BUILD_EXAMPLES})