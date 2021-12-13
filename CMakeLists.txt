cmake_minimum_required(VERSION 3.12)
project(loss_rate)
enable_language(Fortran)

# include JSON-Fortran (my fork thereof)
# https://github.com/jphaupt/json-fortran
# include(ExternalProject)

# set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

# add_executable(FirstTest main.cpp)
# add_dependencies(FirstTest googletest)
# target_link_libraries(FirstTest gtest gtest_main pthread)

# TODO make PFUNIT optional (but recommended of course)
find_package(PFUNIT REQUIRED)
find_package(MPI REQUIRED)
find_package ( jsonfortran-${CMAKE_Fortran_COMPILER_ID} 8.2.5 REQUIRED )
message(STUFF "stuff")
message(STUFF "${jsonfortran_INCLUDE_DIRS}")
include_directories ( "${jsonfortran_INCLUDE_DIRS}" )
enable_testing()

# handle real size
# build with, e.g., -DREAL=128
if(REAL)
  message(STATUS "Configuring build for ${REAL}-bit reals")
  add_definitions(-DREAL${REAL})
#   ExternalProject_Add(jsonfortran
#       GIT_REPOSITORY https://github.com/jphaupt/json-fortran
#       CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DSKIP_DOC_GEN=FALSE -DJSON_REAL_KIND=REAL${REAL}
#   )
else()
  message(STATUS "Configuring build for 32-bit reals")
  add_definitions(-DREAL32)
#   ExternalProject_Add(jsonfortran
#   GIT_REPOSITORY https://github.com/jphaupt/json-fortran
#   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DSKIP_DOC_GEN=FALSE -DJSON_REAL_KIND=REAL32
#   )
endif()


# include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
# link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_subdirectory(src)
# TODO tests
# add_subdirectory(test)
