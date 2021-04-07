cmake_minimum_required(VERSION 3.16)

# FetchContent command is available with cmake >= 3.11
include(FetchContent)

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

find_package(Eigen3 REQUIRED NO_MODULE)
set(CMAKE_CXX_STANDARD 20)

cmake_policy(SET CMP0077 OLD) # let us overwrite options in subprojects

option(BUILD_EXAMPLES "Build examples" ON)
option(USE_COMPACT_N_SEARCH "Build a Pair generator with the help of CompactNSearch" ON)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

if(${USE_COMPACT_N_SEARCH})

  target_compile_definitions(${PROJECT_NAME} INTERFACE -DBUILD_COMPACT_N_SEARCH=1)
  
  FetchContent_Declare(
    compactnsearch #
    GIT_REPOSITORY https://github.com/digitalillusions/CompactNSearch.git #
    GIT_TAG 011eb171c2cbfa536cd47b6a6ffdf90c9b0d40fd #
  )
  set(USE_DOUBLE_PRECISION ON) # Tell CompactNSearch to use doubles instead of floats
  fetch(compactnsearch)
  add_subdirectory(${compactnsearch_SOURCE_DIR} CompactNSearch)
  
  target_link_libraries(${PROJECT_NAME} INTERFACE CompactNSearch)

else(${USE_COMPACT_N_SEARCH})

  target_compile_definitions(${PROJECT_NAME} INTERFACE -DBUILD_COMPACT_N_SEARCH=0)

endif(${USE_COMPACT_N_SEARCH})

if(${BUILD_EXAMPLES})

  add_subdirectory(examples)

  endif(${BUILD_EXAMPLES})