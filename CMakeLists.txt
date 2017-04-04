cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(pydistmesh)

enable_language(C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

add_library(distance_functions SHARED distmesh/src/distance_functions.c)
