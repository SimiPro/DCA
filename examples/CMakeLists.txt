cmake_minimum_required(VERSION 3.16)

project(example)

add_executable(${PROJECT_NAME} fd.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE DCA)

# add_executable(${PROJECT_NAME} example.cpp)

# # Link to DCA and make the header files accessible only to this executable.
# target_link_libraries(${PROJECT_NAME} PRIVATE DCA)

# if(${USE_COMPACT_N_SEARCH})

#   target_compile_definitions(${PROJECT_NAME} PRIVATE -DBUILD_COMPACT_N_SEARCH=1)

# else(${USE_COMPACT_N_SEARCH})

#   target_compile_definitions(${PROJECT_NAME} PRIVATE -DBUILD_COMPACT_N_SEARCH=0)

# endif(${USE_COMPACT_N_SEARCH})