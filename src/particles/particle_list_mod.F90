! this file is for particle lists

MODULE particle_list_mod

    !===================================

    USE MPI_f08

    USE comms_mod

    USE particle_core_mod
    USE particle_dict_mod

    IMPLICIT NONE

    !-----------------------------------

    TYPE :: particle_list_t

        INTEGER(intk) :: iproc           ! REMOVE (obsolete) ?

        INTEGER(intk) :: max_np        ! max number of particles of this process/list
        INTEGER(intk) :: active_np = 0           ! number of active particles of this process/list
        INTEGER(intk) :: ifinal          ! index of last entry of the list which holds an active particle

        TYPE(baseparticle_t), ALLOCATABLE :: particles(:)
        !LOGICAL, ALLOCATABLE :: particle_stored(:)  ! each logical value reflects whether a particle is stored in the list
                                                    ! at the respective index. Is this a feasable and good way to keep track
                                                    ! of particle storage (especially as is_active in particle_t carries the same information?

    END TYPE particle_list_t

    !-----------------------------------

    ! LOGICAL :: dread_particles = .FALSE.

    TYPE(particle_list_t) :: my_particle_list ! rather declare in init_particle_list? NO

    PUBLIC :: my_particle_list

    !===================================

CONTAINS

    SUBROUTINE init_particle_list()

        ! local variables
        INTEGER(intk) :: i, global_np, read_np
        INTEGER(intk), ALLOCATABLE :: ipart_arr(:), p_igrid_arr(:)
        REAL(realk), ALLOCATABLE :: x(:), y(:), z(:)

        my_particle_list%iproc = myid
        my_particle_list%max_np = plist_len

        ALLOCATE(my_particle_list%particles(my_particle_list%max_np))

        IF (dread_particles) THEN

            ALLOCATE(ipart_arr(my_particle_list%max_np))
            ALLOCATE(p_igrid_arr(my_particle_list%max_np))
            ALLOCATE(x(my_particle_list%max_np))
            ALLOCATE(y(my_particle_list%max_np))
            ALLOCATE(z(my_particle_list%max_np))

            ! all arguments directly taken from the particle lsit are input only
            CALL read_particles(dread_particles, my_particle_list%max_np, global_np, ipart_arr, p_igrid_arr, x, y, z, read_np)

            IF (dread_particles) THEN

                IF (myid == 0) THEN
                    SELECT CASE (TRIM(particle_terminal))
                        CASE ("none")
                            CONTINUE
                        CASE ("normal")
                            WRITE(*,*) ' '
                            WRITE(*, '("Initializing ", I0, " Particle(s):")') global_np
                        CASE ("verbose")
                            WRITE(*,*) ' '
                            WRITE(*, '("Initializing ", I0, " Particle(s):")') global_np
                    END SELECT
                END IF

                DO i = 1, read_np

                    my_particle_list%active_np = my_particle_list%active_np + 1_intk

                    CALL set_particle(my_particle_list%particles(i), &
                    ipart = ipart_arr(i), x = x(i), y = y(i), z = z(i), igrid = p_igrid_arr(i))

                    CALL print_particle_status(my_particle_list%particles(i))

                END DO

                my_particle_list%ifinal = my_particle_list%active_np

            END IF

            DEALLOCATE(ipart_arr)
            DEALLOCATE(p_igrid_arr)
            DEALLOCATE(x)
            DEALLOCATE(y)
            DEALLOCATE(z)

        END IF

        IF (.NOT. dread_particles) THEN

            my_particle_list%active_np = init_npart

            ALLOCATE(ipart_arr(my_particle_list%active_np))
            ALLOCATE(p_igrid_arr(my_particle_list%active_np))
            ALLOCATE(x(my_particle_list%active_np))
            ALLOCATE(y(my_particle_list%active_np))
            ALLOCATE(z(my_particle_list%active_np))

            CALL dist_ipart(my_particle_list%active_np, ipart_arr)

            CALL dist_particles(my_particle_list%active_np, p_igrid_arr, x, y, z)

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*,*) ' '
                        WRITE(*, '("Initializing ", I0, " Particle(s):")') my_particle_list%active_np
                    CASE ("verbose")
                        WRITE(*,*) ' '
                        WRITE(*, '("Initializing ", I0, " Particle(s):")') my_particle_list%active_np
                END SELECT
            END IF

            ! BARRIER ONLY FOR DEGUGGING -- TEMPORARY <----------------------------------------------- TODO : remove
            CALL MPI_Barrier(MPI_COMM_WORLD)

            my_particle_list%ifinal = my_particle_list%active_np

            DO i = 1, my_particle_list%active_np

                CALL set_particle( my_particle_list%particles(i), &
                    ipart = ipart_arr(i), x = x(i), y = y(i), z = z(i), &
                    igrid = p_igrid_arr(i))

                CALL print_particle_status( my_particle_list%particles(i) )

            END DO

            DEALLOCATE(ipart_arr)
            DEALLOCATE(p_igrid_arr)
            DEALLOCATE(x)
            DEALLOCATE(y)
            DEALLOCATE(z)

        END IF

        SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) ' '
                    WRITE(*, *) "Particle list of length ", my_particle_list%max_np, " initialized on process ", myid, "."
                CASE ("verbose")
                    WRITE(*, *) ' '
                    WRITE(*, *) "Particle list of length ", my_particle_list%max_np, " initialized on process ", myid, "."
        END SELECT

        ! BARRIER ONLY FOR DEGUGGING -- TEMPORARY <----------------------------------------------- TODO : remove
        CALL MPI_Barrier(MPI_COMM_WORLD)

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*,*) ' '
                    WRITE(*, '("Initialization of Particle(s) finished successfully.")')
                    WRITE(*,*) ' '
                CASE ("verbose")
                    WRITE(*,*) ' '
                    WRITE(*, '("Initialization of Particle(s) finished successfully.")')
                    WRITE(*,*) ' '
            END SELECT
        END IF

    END SUBROUTINE init_particle_list

    !-----------------------------------

    ! maybe this way of adding length to the particle list should be replaced by a smart pointer approach for the particle list type?
    SUBROUTINE reallocate_particle_list(particle_list, add_len)

        ! subroutine arguments
        TYPE(particle_list_t), INTENT(inout) :: particle_list
        INTEGER, INTENT(in) :: add_len ! can also be negative to shorten the particle list

        !local variables
        INTEGER(intk) :: i

        ! local variable / this is a new list of particles (not the particle_list_t) with additional length compared to the particle list length
        TYPE(baseparticle_t), ALLOCATABLE :: particles_tmp(:)

        ALLOCATE(particles_tmp(particle_list%max_np + add_len))

        DO i = 1, particle_list%ifinal
            particles_tmp(i) = particle_list%particles(i) !copy particles into temporary list
        END DO

        ! might be unnessecary, but ensures that all new particle slots are inactive
        DO i = particle_list%ifinal + 1, particle_list%max_np + add_len
            particles_tmp(i)%is_active = 0
        END DO

        particle_list%max_np = particle_list%max_np + add_len
        particle_list%ifinal = MIN(particle_list%ifinal, particle_list%max_np)

        CALL MOVE_ALLOC(particles_tmp, particle_list%particles) ! "rename" / make my

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                WRITE(*, *) ' '
                WRITE(*, *) "Particle list enlarged by ", add_len, " on process ", myid, "."
            CASE ("verbose")
                WRITE(*, *) ' '
                WRITE(*, *) "Particle list enlarged by ", add_len, " on process ", myid, "."
        END SELECT

    END SUBROUTINE reallocate_particle_list

    !-----------------------------------

    SUBROUTINE dist_ipart(npart, ipart_arr) ! this routine is supposed to hand out a list of unique particle ids (ipart) to every process ! ONLY NON MPI RUNS FOR NOW

        ! subroutine arguments
        INTEGER(intk), INTENT(IN) :: npart
        INTEGER(intk), INTENT(out) :: ipart_arr(npart)

        ! local variables
        INTEGER(intk) :: i, counter

        counter = 1
        DO i = myid * npart + 1, myid * npart + npart
            ipart_arr(counter) = i
            counter = counter + 1_intk
        END DO

        ! MPI barrier

    END SUBROUTINE dist_ipart

    !-----------------------------------

    SUBROUTINE dist_particles(npart, p_igrid_arr, x, y, z)

        ! subroutine arguments...
        INTEGER(intk), INTENT(in) :: npart
        INTEGER(intk), INTENT(inout) :: p_igrid_arr(npart)
        REAL(realk), INTENT(inout) :: x(npart), y(npart), z(npart)

        ! local variables...
        INTEGER(intk) :: i, j, igrid
        REAL(realk) :: myvolume, volume_fractions(nmygrids), grid_rn
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        ! compute combined volume of all grids this process owns (for now assumning ther is only one level)...

        myvolume = 0.0_realk

        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            myvolume = myvolume + (maxx - minx) * (maxy - miny) * (maxz - minz)

        END DO

        ! compute fraction of the combined volume of each grids volume (for now assumning ther is only one level)...

        igrid = mygrids(1)
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        volume_fractions(1) = (maxx - minx) * (maxy - miny) * (maxz - minz) / myvolume

        volume_fractions(1) = MAX(0.0_realk, volume_fractions(1))
        volume_fractions(1) = MIN(1.0_realk, volume_fractions(1))

        DO i = 2, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            volume_fractions(i) = volume_fractions(i-1) + ((maxx - minx) * (maxy - miny) * (maxz - minz) / myvolume)

            volume_fractions(i) = MAX(0.0_realk, volume_fractions(i))
            volume_fractions(i) = MIN(1.0_realk, volume_fractions(i))

        END DO

        ! distribute particles along all grids this process owns (uniformely distibuted)

        DO j = 1, npart

            i = 1_intk

            CALL RANDOM_NUMBER(grid_rn)

            DO WHILE (grid_rn > volume_fractions(i))

                i = i + 1_intk

            END DO

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            p_igrid_arr(j) = igrid

            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(x(j))
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(y(j))
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(z(j))

            x(j) = minx + x(j) * (maxx - minx)
            y(j) = miny + y(j) * (maxy - miny)
            z(j) = minz + z(j) * (maxz - minz)

        END DO

    END SUBROUTINE dist_particles

    !===================================

    SUBROUTINE print_list_status(particle_list)

        ! subroutine arguments
        TYPE(particle_list_t), INTENT(in) :: particle_list


        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                WRITE(*, *) "Particle list status on process ", myid, ": max_np = ", &
                 particle_list%max_np, ", active_np = ", particle_list%active_np, ", ifinal = ", particle_list%ifinal
                IF (myid == 0) THEN
                    WRITE(*, *) " "
                END IF
            CASE ("verbose")
                WRITE(*, *) "Particle list status on process ", myid, ": max_np = ", &
                 particle_list%max_np, ", active_np = ", particle_list%active_np, ", ifinal = ", particle_list%ifinal
                IF (myid == 0) THEN
                    WRITE(*, *) " "
                END IF
        END SELECT

    END SUBROUTINE print_list_status

END MODULE
