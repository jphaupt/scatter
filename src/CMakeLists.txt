add_library(
    sut
    precision.F90
    constants.f90
    special_functions.f90
    potential_type.f90
    integrate.f90
)

# not totally sure what these two do, but they are necessary for testing
set_target_properties (sut PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(sut PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_executable(main main.f90)

target_include_directories(sut PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(sut PUBLIC ${MPI_Fortran_LIBRARIES})
target_link_libraries(main sut jsonfortran)
# add_dependencies(main jsonfortran)
