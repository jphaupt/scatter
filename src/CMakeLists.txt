add_library(
    sut
    integrate.f90
    # TODO
)

add_executable(main main.f90)

target_link_libraries(main sut)
add_dependencies(main jsonfortran)
