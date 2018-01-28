cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(custom_functions CXX)

add_library(
    custom_functions
    SHARED
    custom_functions.cpp
    )
