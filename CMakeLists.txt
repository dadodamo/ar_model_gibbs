cmake_minimum_required(VERSION 3.24)
project(ar_gibbs)

set(CMAKE_CXX_STANDARD 17)

## External libraries
set(EIGEN_DIR "/users/daniel/Desktop/libraries/Eigen3")
set(BOOST_DIR "/users/daniel/Desktop/libraries/boost")
set(H_FILES "/users/daniel/Desktop/libraries/h_files")
if(NOT EXISTS ${EIGEN_DIR})
    message(FATAL_ERROR "Please check that the set Eigen directory is valid!")
endif()
if(NOT EXISTS ${BOOST_DIR})
    message(FATAL_ERROR "Please check that the set boost directory is valid!")
endif()
if(NOT EXISTS ${H_FILES})
    message(FATAL_ERROR "Please check that the set sampler directory is valid!")
endif()
include_directories(${EIGEN_DIR})
include_directories(${BOOST_DIR})
include_directories(${H_FILES})

add_executable(ar_gibbs
        main.cpp
        ${H_FILES}/eigenmvn.h
        matern.h
        coordinates.h dist_sample.cpp
        gibbs_sampler.h
        gibbs_sampler.cpp)
