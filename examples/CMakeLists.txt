cmake_minimum_required(VERSION 3.16)

project(example)

add_executable(fd fd.cpp)
target_link_libraries(fd PUBLIC DCA)

add_executable(example example.cpp)
target_link_libraries(example PUBLIC DCA)
