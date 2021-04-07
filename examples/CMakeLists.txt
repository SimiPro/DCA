cmake_minimum_required(VERSION 3.16)

project(examples)

add_executable(example example.cpp)

# Link to DCA and make the header files accessible only to this executable.
target_link_libraries(example PRIVATE DCA)
