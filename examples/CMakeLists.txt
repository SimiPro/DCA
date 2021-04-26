cmake_minimum_required(VERSION 3.16)

project(example)

add_executable(fd fd.cpp)
target_link_libraries(fd PRIVATE DCA)

add_executable(example example.cpp)

# Link to DCA and make the header files accessible only to this executable.
target_link_libraries(example PRIVATE DCA)

if(${DCA_USE_COMPACT_N_SEARCH})
  target_compile_definitions(example PRIVATE -DBUILD_COMPACT_N_SEARCH=1)
else(${DCA_USE_COMPACT_N_SEARCH})
  target_compile_definitions(example PRIVATE -DBUILD_COMPACT_N_SEARCH=0)
endif(${DCA_USE_COMPACT_N_SEARCH})