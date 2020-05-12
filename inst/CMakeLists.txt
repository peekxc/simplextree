cmake_minimum_required(VERSION 3.0.0)
project(st_test VERSION 0.1.0)

# set(cmake.cmakePath "/usr/local/Cellar/cmake/3.16.1")
SET(CMAKE_CXX_STANDARD 17)
SET(CXX_COMPILE_FLAGS "-Wall -pedantic")

## If using armadillo 
# SET(ARMA_CXX_FLAGS "-DARMA_DONT_USE_WRAPPER -framework Accelerate")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CXX_COMPILE_FLAGS}") # ${ARMA_CXX_FLAGS}

# include(CTest)
# enable_testing()

add_executable(st_test st_test.cpp)
target_include_directories(st_test PUBLIC include)

# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)
