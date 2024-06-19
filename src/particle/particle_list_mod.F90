! this file is for particle lists

MODULE particle_list_mod

    USE particlecore_mod

    IMPLICIT NONE

    TYPE :: particle_list_t

        INTEGER(intk) :: iproc	            ! REMOVE (obsolete) ?

        INTEGER(intk) :: max_np    			! max number of particles of this process/list
        INTEGER(intk) :: active_np 			! number of active particles of this process/list
        INTEGER(intk) :: ifinal 			! index of last entry of the list which holds an active particle

        TYPE(baseparticle_t), ALLOCATABLE :: particles(:)
        LOGICAL, ALLOCATABLE :: particle_stored(:) 	! each logical value reflects whether a particle is stored in the list
                                                    ! at the respective index. Is this a feasable and good way to keep track
                                                    ! of particle storage (especially as is_init in particle_t carries the same information?

    END TYPE particle_list_t

    !===================================

    ! module variables

    INTEGER(intk), PARAMETER :: default_max_np = 1000
    INTEGER(intk), PARAMETER :: default_initial_np = 100 !ONLY DUMMY VALUE FOR NOW, SHOULD BE SCALED WITH THE SIZE OF THE SPATIAL DOMAIN THAT THE PROCESS HANDLES

    TYPE(particle_list_t) :: my_particle_list ! rather declare in init_particles?

    PUBLIC :: my_particle_list

    !===================================

CONTAINS

    SUBROUTINE init_particles()

        ! local variables
         INTEGER(intk) :: i
         INTEGER(intk), ALLOCATABLE :: ipart_arr(:)
         REAL(realk) :: x, y, z

        my_particle_list%max_np = default_max_np
         my_particle_list%active_np = default_initial_np
        my_particle_list%ifinal = default_initial_np
        my_particle_list%iproc = myid

        ALLOCATE(my_particle_list%particles(my_particle_list%max_np))
        ALLOCATE(my_particle_list%particle_stored(my_particle_list%max_np))

        CALL dist_ipart(ipart_arr)

        my_particle_list%particle_stored = .FALSE.

         DO i = 1, my_particle_list%active_np

             CALL random_ic(x, y, z) ! ONLY DUMMY FOR NOW, DENPENDS ON THE PROCESS SPATIAL DOMAIN

             CALL my_particle_list%particles(i)%init(ipart_arr(i), x, y, z)

             my_particle_list%particle_stored(i) = .TRUE.

         END DO


    END SUBROUTINE init_particles

    !------------------------------

    SUBROUTINE dist_ipart(ipart_arr) ! this routine is supposed to hand out a list of unique particle ids (ipart) to every process ! ONLY NON MPI RUNS FOR NOW

        ! subroutine arguments
        INTEGER(intk), ALLOCATABLE, INTENT(out) :: ipart_arr(:)

        ! local variables
        INTEGER(intk) :: i, count

        ALLOCATE(ipart_arr(default_initial_np))

        count = 1
        DO i = myid * default_initial_np + 1, myid * default_initial_np + default_initial_np
            ipart_arr(count) = i
            count = count + 1
        END DO

        ! MPI barrier

    END SUBROUTINE dist_ipart


END MODULE
