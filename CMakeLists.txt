cmake_minimum_required(VERSION 2.6)

project(simple_kalman)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Enable this to run the first example
add_executable(simple_kalman main1.cpp simple_kalman.cpp ekf.cpp)

# Enable this to run the second example
#add_executable(simple_kalman main2.cpp simple_kalman.cpp) 

# Enable this to run the third example
#add_executable(simple_kalman main3.cpp simple_kalman.cpp)

# Enable this to run the third example
#add_executable(simple_kalman main4.cpp simple_kalman.cpp ekf.cpp)


install(TARGETS simple_kalman RUNTIME DESTINATION bin)

target_link_libraries (simple_kalman armadillo)	