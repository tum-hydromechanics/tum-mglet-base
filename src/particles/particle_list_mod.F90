MODULE particle_list_mod

    ! This module is responsible for:
    ! Definition of the particle_list_t type.
    ! Storage of Particles in my_particle_list.
    ! Basic list manipulations.
    ! Initialization of partticles.

    USE MPI_f08

    USE comms_mod

    USE particle_core_mod
    USE particle_dict_mod

    IMPLICIT NONE

    TYPE :: particle_list_t

        ! TODO: REMOVE (obsolete) ?
        INTEGER(intk) :: iproc

        ! max number of particles of this process/list
        INTEGER(intk) :: max_np
        ! number of active particles of this process/list
        INTEGER(intk) :: active_np = 0
        ! index of last entry of the list which holds an active particle
        INTEGER(intk) :: ifinal

        ! array that hold the actual particles
        TYPE(baseparticle_t), ALLOCATABLE :: particles(:)

    END TYPE particle_list_t

    TYPE(particle_list_t) :: my_particle_list

    PUBLIC :: my_particle_list

CONTAINS    !===================================

    SUBROUTINE init_particle_list()

        ! local variables
        INTEGER(intk) :: i, global_np, read_np
        INTEGER(intk), ALLOCATABLE :: ipart_arr(:), p_igrid_arr(:)
        REAL(realk), ALLOCATABLE :: x(:), y(:), z(:)

        CALL start_timer(900)
        CALL start_timer(910)

        my_particle_list%iproc = myid
        my_particle_list%max_np = plist_len

        ALLOCATE(my_particle_list%particles(my_particle_list%max_np))

        IF (dread_particles) THEN

            ALLOCATE(ipart_arr(my_particle_list%max_np))
            ALLOCATE(p_igrid_arr(my_particle_list%max_np))
            ALLOCATE(x(my_particle_list%max_np))
            ALLOCATE(y(my_particle_list%max_np))
            ALLOCATE(z(my_particle_list%max_np))

            ! all arguments that are passed as particle list attributes are input only
            CALL read_particles(dread_particles, my_particle_list%max_np, global_np, ipart_arr, p_igrid_arr, x, y, z, read_np)

            IF (dread_particles) THEN

                IF (myid == 0) THEN
                    SELECT CASE (TRIM(particle_terminal))
                        CASE ("none")
                            CONTINUE
                        CASE ("normal")
                            WRITE(*,*) ' '
                            WRITE(*, '("INITIALIZING ", I0, " PARTICLE(S):")') global_np
                            WRITE(*, '()')
                        CASE ("verbose")
                            WRITE(*, '("INITIALIZING ", I0, " PARTICLE(S):")') global_np
                            WRITE(*, '()')
                    END SELECT
                END IF

                DO i = 1, read_np

                    my_particle_list%active_np = my_particle_list%active_np + 1_intk

                    CALL set_particle(my_particle_list%particles(i), &
                    ipart = ipart_arr(i), x = x(i), y = y(i), z = z(i), igrid = p_igrid_arr(i))

                    SELECT CASE (TRIM(particle_terminal))
                        CASE ("none")
                            CONTINUE
                        CASE ("normal")
                            CONTINUE
                        CASE ("verbose")
                            CALL print_particle_status(my_particle_list%particles(i))
                            WRITE(*, '()')
                    END SELECT

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
                        WRITE(*, '("INITIALIZING ", I0, " PARTICLE(S):")') my_particle_list%active_np
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, '("INITIALIZING ", I0, " PARTICLE(S):")') my_particle_list%active_np
                        WRITE(*, '()')
                END SELECT
            END IF

            ! BARRIER ONLY FOR DEGUGGING -- TEMPORARY <----------------------------------------------- TODO : remove
            CALL MPI_Barrier(MPI_COMM_WORLD)

            my_particle_list%ifinal = my_particle_list%active_np

            DO i = 1, my_particle_list%active_np

                CALL set_particle(my_particle_list%particles(i), &
                    ipart = ipart_arr(i), x = x(i), y = y(i), z = z(i), &
                    igrid = p_igrid_arr(i))

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        CALL print_particle_status(my_particle_list%particles(i))
                    WRITE(*, '()')
                END SELECT

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
                    WRITE(*, *) "Particle list of length ", my_particle_list%max_np, " initialized on process ", myid, "."
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) "Particle list of length ", my_particle_list%max_np, " initialized on process ", myid, "."
                    WRITE(*, '()')
        END SELECT

        ! BARRIER ONLY FOR DEGUGGING -- TEMPORARY <----------------------------------------------- TODO : remove
        CALL MPI_Barrier(MPI_COMM_WORLD)

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, '("INITIALIZATION OF PARTICLE(S) SUCCESSFULLY COMPLETED.")')
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, '("INITIALIZATION OF PARTICLE(S) SUCCESSFULLY COMPLETED.")')
                    WRITE(*, '()')
            END SELECT
        END IF

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_list

    !-----------------------------------

    SUBROUTINE reallocate_particle_list(particle_list, add_len)

        ! subroutine arguments
        TYPE(particle_list_t), INTENT(inout) :: particle_list
        ! additional length (can also be negative to shorten the particle list!)
        INTEGER, INTENT(in) :: add_len

        !local variables
        INTEGER(intk) :: i
        !new list of particles (not the particle_list_t) with additional length relative to particle_list%max_np
        TYPE(baseparticle_t), ALLOCATABLE :: particles_tmp(:)

        ALLOCATE(particles_tmp(particle_list%max_np + add_len))

        !copy particles into temporary list
        DO i = 1, particle_list%ifinal
            particles_tmp(i) = particle_list%particles(i)
        END DO

        ! might be unnessecary, but ensures that all new particle slots are inactive
        DO i = particle_list%ifinal + 1, particle_list%max_np + add_len
            particles_tmp(i)%state = -1
        END DO

        particle_list%max_np = particle_list%max_np + add_len
        particle_list%ifinal = MIN(particle_list%ifinal, particle_list%max_np)

        CALL MOVE_ALLOC(particles_tmp, particle_list%particles)

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                WRITE(*, *) "Particle list enlarged by ", add_len, " on process ", myid, "."
                WRITE(*, '()')
            CASE ("verbose")
                WRITE(*, *) "Particle list enlarged by ", add_len, " on process ", myid, "."
                WRITE(*, '()')
        END SELECT

    END SUBROUTINE reallocate_particle_list

    !-----------------------------------

    !this routine is supposed to hand out a list of unique particle ids (ipart) to every process ! ONLY NON MPI RUNS FOR NOW
    SUBROUTINE dist_ipart(npart, ipart_arr)

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

        ! compute combined volume of all grids this process owns (assumning ther is only one level)
        myvolume = 0.0_realk

        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            myvolume = myvolume + (maxx - minx) * (maxy - miny) * (maxz - minz)

        END DO

        ! compute fraction of the combined volume of each grids volume (assumning ther is only one level)
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

        ! distribute particles among all grids this process owns (uniformely distibuted)
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

    !-----------------------------------

    SUBROUTINE print_list_status(particle_list)

        ! subroutine arguments
        TYPE(particle_list_t), INTENT(in) :: particle_list

        WRITE(*, *) "Particle list status on process ", myid, ": max_np = ", &
         particle_list%max_np, ", active_np = ", particle_list%active_np, ", ifinal = ", particle_list%ifinal

        IF (myid == 0) THEN
            WRITE(*, '()')
        END IF

    END SUBROUTINE print_list_status

    SUBROUTINE finish_particle_list()

        DEALLOCATE(my_particle_list%particles)

    END SUBROUTINE finish_particle_list

END MODULE
