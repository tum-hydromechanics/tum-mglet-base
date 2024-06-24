! this file is for particle lists

MODULE particle_list_mod

    !===================================

    USE particlecore_mod

    IMPLICIT NONE

    !-----------------------------------

    TYPE :: particle_list_t

        INTEGER(intk) :: iproc           ! REMOVE (obsolete) ?

        INTEGER(intk) :: max_np          ! max number of particles of this process/list
        INTEGER(intk) :: active_np       ! number of active particles of this process/list
        INTEGER(intk) :: ifinal          ! index of last entry of the list which holds an active particle

        TYPE(baseparticle_t), ALLOCATABLE :: particles(:)
        !LOGICAL, ALLOCATABLE :: particle_stored(:)  ! each logical value reflects whether a particle is stored in the list
                                                    ! at the respective index. Is this a feasable and good way to keep track
                                                    ! of particle storage (especially as is_active in particle_t carries the same information?

    END TYPE particle_list_t

    !-----------------------------------

    INTEGER(intk), PARAMETER :: default_max_np = 1000
    INTEGER(intk), PARAMETER :: default_initial_np = 100 !ONLY DUMMY VALUE FOR NOW

    TYPE(particle_list_t) :: my_particle_list ! rather declare in init_particle_list? NO

    PUBLIC :: my_particle_list

    !===================================

CONTAINS

    SUBROUTINE init_particle_list()

        ! local variables
        INTEGER(intk) :: i
        INTEGER(intk), ALLOCATABLE :: ipart_arr(:), p_igrid_arr(:)
        REAL(realk), ALLOCATABLE :: x(:), y(:), z(:)

        my_particle_list%max_np = default_max_np
        my_particle_list%active_np = default_initial_np
        my_particle_list%ifinal = default_initial_np
        my_particle_list%iproc = myid

        ALLOCATE(my_particle_list%particles(my_particle_list%max_np))
        !ALLOCATE(my_particle_list%particle_stored(my_particle_list%max_np))

        CALL dist_ipart(ipart_arr)

        CALL dist_part(my_particle_list%active_np, p_igrid_arr, x, y, z)

        !my_particle_list%particle_stored = .FALSE.

         DO i = 1, my_particle_list%active_np

             CALL my_particle_list%particles(i)%init(ipart = ipart_arr(i), x = x(i), y = y(i), z = z(i), igrid = p_igrid_arr(i))

             !my_particle_list%particle_stored(i) = .TRUE.

         END DO


    END SUBROUTINE init_particle_list

    !-----------------------------------

    SUBROUTINE dist_ipart(ipart_arr) ! this routine is supposed to hand out a list of unique particle ids (ipart) to every process ! ONLY NON MPI RUNS FOR NOW

        ! subroutine arguments
        INTEGER(intk), ALLOCATABLE, INTENT(out) :: ipart_arr(:)

        ! local variables
        INTEGER(intk) :: i, counter

        ALLOCATE(ipart_arr(default_initial_np))

        counter = 1
        DO i = myid * default_initial_np + 1, myid * default_initial_np + default_initial_np
            ipart_arr(counter) = i
            counter = counter + 1_intk
        END DO

        ! MPI barrier

    END SUBROUTINE dist_ipart

    !-----------------------------------

    SUBROUTINE dist_part(npart, p_igrid_arr, x, y, z)

        ! subroutine arguments...
        INTEGER(intk), INTENT(in) :: npart
        INTEGER(intk), ALLOCATABLE, INTENT(inout) :: p_igrid_arr(:) !p_igrid_arr(npart)
        REAL(realk), ALLOCATABLE, INTENT(inout) :: x(:), y(:), z(:) !x(npart), y(npart), z(npart)

        ! local variables...
        INTEGER(intk) :: i, j, igrid
        REAL(realk) :: myvolume, volume_fractions(nmygrids), grid_rn
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        ALLOCATE(p_igrid_arr(npart))
        ALLOCATE(x(npart))
        ALLOCATE(y(npart))
        ALLOCATE(z(npart))

        myvolume = 0

        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            myvolume = myvolume + (maxx - minx) * (maxy - miny) * (maxz - minz)

        END DO

        igrid = mygrids(1)
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        volume_fractions(1) = (maxx - minx) * (maxy - miny) * (maxz - minz) / myvolume

        DO i = 2, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            volume_fractions(i) = volume_fractions(i-1) + ((maxx - minx) * (maxy - miny) * (maxz - minz) / myvolume)

        END DO

        DO j = 1, npart

            i = 1_intk

            CALL RANDOM_NUMBER(grid_rn)

            DO WHILE (grid_rn > volume_fractions(i))

                i = i + 1_intk

            END DO

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            p_igrid_arr(j) = igrid

            CALL RANDOM_NUMBER(x(j))
            CALL RANDOM_NUMBER(y(j))
            CALL RANDOM_NUMBER(z(j))

            x(j) = minx + x(j) * (maxx - minx)
            y(j) = miny + y(j) * (maxy - miny)
            z(j) = minz + z(j) * (maxz - minz)

        END DO

    END SUBROUTINE dist_part

    !===================================

END MODULE
