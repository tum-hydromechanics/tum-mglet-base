add_library(particle STATIC
    particlecore_mod.F90
    particle_list_mod.F90
)

target_include_directories(particle
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(particle
    PRIVATE core
    PRIVATE ${MPI_Fortran_LIBRARIES}
)
