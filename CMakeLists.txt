cmake_minimum_required(VERSION 2.6)
project(simple_kalman)

add_executable(simple_kalman main.cpp simple_kalman.cpp)

install(TARGETS simple_kalman RUNTIME DESTINATION bin)

target_link_libraries (simple_kalman armadillo)	