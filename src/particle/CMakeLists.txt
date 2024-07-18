add_library(particle STATIC
    particle_config_mod.F90
    particlecore_mod.F90
    particle_list_mod.F90
    particle_dict_mod.F90
    particle_snapshot_mod.F90
    particle_connect_mod.F90
    particle_timeintegration_mod.F90
    particle_mod.F90
)

target_include_directories(particle
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(particle
    PRIVATE core
    PRIVATE flow
    PRIVATE ${MPI_Fortran_LIBRARIES}
)
