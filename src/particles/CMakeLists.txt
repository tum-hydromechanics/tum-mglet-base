add_library(particles STATIC
    particle_config_mod.F90
    particle_core_mod.F90
    particle_list_mod.F90
    particle_dict_mod.F90
    particle_snapshot_mod.F90
    particle_exchange_mod.F90
    particle_interpolation_mod.F90
    particle_timeintegration_mod.F90
    particle_mod.F90
    particle_interpolation_mod.F90
)

target_include_directories(particles
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(particles
    PRIVATE core
    PRIVATE flow
    PRIVATE ${MPI_Fortran_LIBRARIES}
)
