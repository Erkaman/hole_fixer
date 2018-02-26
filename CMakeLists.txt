cmake_minimum_required(VERSION 3.1)
project(hole_fixer)

include_directories(libigl/ libigl/eigen/)
add_executable(hole_fixer main.cpp)

