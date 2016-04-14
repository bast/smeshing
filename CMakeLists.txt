cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(pydistmesh)

enable_language(C)

add_library(distance_functions SHARED distmesh/src/distance_functions.c)

# FIXME hardcoded, use find_library instead
target_link_libraries(distance_functions /usr/lib/liblapack.so)
