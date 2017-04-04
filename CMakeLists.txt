cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(pydistmesh)

enable_language(CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_library(distance_functions SHARED distmesh/src/distance_functions.cpp)
