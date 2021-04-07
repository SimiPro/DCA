cmake_minimum_required(VERSION 3.16)

project("DCA" VERSION 0.1
             DESCRIPTION "Differentiable Collision Avoidance"
             HOMEPAGE_URL "https://github.com/SimiPro/DCA")

find_package(Eigen3 REQUIRED NO_MODULE)
set(CMAKE_CXX_STANDARD 20)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

option(BUILD_EXAMPLES "Build examples" ON)

if(${BUILD_EXAMPLES})
  add_subdirectory(examples)
endif(${BUILD_EXAMPLES})