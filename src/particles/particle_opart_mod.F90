MODULE particle_opart_mod
    USE precision_mod, ONLY: intk, realk
    USE err_mod, ONLY: errr

    USE particle_config_mod
    USE particle_basetype_mod
    USE particle_list_mod

    IMPLICIT NONE

    INTEGER(intk), ALLOCATABLE :: particle_ipart_arr(:)
    INTEGER(intk), ALLOCATABLE :: particle_igrid_arr(:)
    INTEGER(intk), ALLOCATABLE :: particle_icell_arr(:)
    INTEGER(intk), ALLOCATABLE :: particle_jcell_arr(:)
    INTEGER(intk), ALLOCATABLE :: particle_kcell_arr(:)
    INTEGER(intk), ALLOCATABLE :: particle_seed_arr(:)

    REAL(realk), ALLOCATABLE :: particle_x_arr(:)
    REAL(realk), ALLOCATABLE :: particle_y_arr(:)
    REAL(realk), ALLOCATABLE :: particle_z_arr(:)

    !$omp declare target(particle_ipart_arr, particle_igrid_arr)
    !$omp declare target(particle_icell_arr, particle_jcell_arr, particle_kcell_arr)
    !$omp declare target(particle_seed_arr)
    !$omp declare target(particle_x_arr, particle_y_arr, particle_z_arr)

CONTAINS

    SUBROUTINE init_particle_arrays()

        INTEGER(intk) :: array_length

        array_length = my_particle_list%max_np

        ALLOCATE(particle_ipart_arr(array_length))
        ALLOCATE(particle_igrid_arr(array_length))
        ALLOCATE(particle_icell_arr(array_length))
        ALLOCATE(particle_jcell_arr(array_length))
        ALLOCATE(particle_kcell_arr(array_length))
        ALLOCATE(particle_seed_arr(array_length))
        ALLOCATE(particle_x_arr(array_length))
        ALLOCATE(particle_y_arr(array_length))
        ALLOCATE(particle_z_arr(array_length))

        !$omp target enter data map(alloc: particle_ipart_arr, particle_igrid_arr)
        !$omp target enter data map(alloc: particle_icell_arr, particle_jcell_arr, particle_kcell_arr)
        !$omp target enter data map(alloc: particle_seed_arr)
        !$omp target enter data map(alloc: particle_x_arr, particle_y_arr, particle_z_arr)

    END SUBROUTINE init_particle_arrays    


    SUBROUTINE copy_particle_data_to()

        INTEGER(intk) :: i

        DO i = 1, my_particle_list%ifinal

            particle_ipart_arr(i) = my_particle_list%particles(i)%ipart
            particle_igrid_arr(i) = my_particle_list%particles(i)%igrid
            particle_icell_arr(i) = my_particle_list%particles(i)%ijkcell(1)
            particle_jcell_arr(i) = my_particle_list%particles(i)%ijkcell(2)
            particle_kcell_arr(i) = my_particle_list%particles(i)%ijkcell(3)
            particle_seed_arr(i) = my_particle_list%particles(i)%seed
            particle_x_arr(i) = my_particle_list%particles(i)%x
            particle_y_arr(i) = my_particle_list%particles(i)%y
            particle_z_arr(i) = my_particle_list%particles(i)%z

        END DO

        !$omp target update to(particle_ipart_arr, particle_igrid_arr)
        !$omp target update to(particle_icell_arr, particle_jcell_arr, particle_kcell_arr)
        !$omp target update to(particle_seed_arr)
        !$omp target update to(particle_x_arr, particle_y_arr, particle_z_arr)

    END SUBROUTINE copy_particle_data_to

    SUBROUTINE copy_particle_data_from()

        INTEGER(intk) :: i

        !$omp target update from(particle_ipart_arr, particle_igrid_arr)
        !$omp target update from(particle_icell_arr, particle_jcell_arr, particle_kcell_arr)
        !$omp target update from(particle_seed_arr)
        !$omp target update from(particle_x_arr, particle_y_arr, particle_z_arr)

        DO i = 1, my_particle_list%ifinal

            my_particle_list%particles(i)%ipart = particle_ipart_arr(i)
            my_particle_list%particles(i)%igrid = particle_igrid_arr(i)
            my_particle_list%particles(i)%ijkcell(1) = particle_icell_arr(i)
            my_particle_list%particles(i)%ijkcell(2) = particle_jcell_arr(i)
            my_particle_list%particles(i)%ijkcell(3) = particle_kcell_arr(i)
            my_particle_list%particles(i)%seed = particle_seed_arr(i)
            my_particle_list%particles(i)%x = particle_x_arr(i)
            my_particle_list%particles(i)%y = particle_y_arr(i)
            my_particle_list%particles(i)%z = particle_z_arr(i)

        END DO

    END SUBROUTINE copy_particle_data_from


    SUBROUTINE finish_particle_arrays()

        !$omp target exit data map(delete: particle_ipart_arr, particle_igrid_arr)
        !$omp target exit data map(delete: particle_icell_arr, particle_jcell_arr, particle_kcell_arr)
        !$omp target exit data map(delete: particle_seed_arr)
        !$omp target exit data map(delete: particle_x_arr, particle_y_arr, particle_z_arr)
        
        DEALLOCATE(particle_ipart_arr)
        DEALLOCATE(particle_igrid_arr)
        DEALLOCATE(particle_icell_arr)
        DEALLOCATE(particle_jcell_arr)
        DEALLOCATE(particle_kcell_arr)
        DEALLOCATE(particle_seed_arr)
        DEALLOCATE(particle_x_arr)
        DEALLOCATE(particle_y_arr)
        DEALLOCATE(particle_z_arr)

    END SUBROUTINE finish_particle_arrays

END MODULE particle_opart_mod
