cmake_minimum_required(VERSION 3.12)
project(loss_rate)
enable_language(Fortran)

find_package(PFUNIT REQUIRED)
find_package(MPI REQUIRED)
enable_testing()

# handle real size
# build with, e.g., -DREAL=128
if(REAL)
  message(STATUS "Configuring build for ${REAL}-bit reals")
  add_definitions(-DREAL${REAL})
else()
  message(STATUS "Configuring build for 32-bit reals")
  add_definitions(-DREAL32)
endif()

add_subdirectory(src)
# TODO tests
# add_subdirectory(test)
