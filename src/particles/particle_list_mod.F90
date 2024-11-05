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
    USE particle_utils_mod
    USE particle_boundaries_mod

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

        CONTAINS

            PROCEDURE :: defragment

    END TYPE particle_list_t

    TYPE(particle_list_t) :: my_particle_list

    PUBLIC :: my_particle_list

CONTAINS    !===================================

    SUBROUTINE init_particle_list()

        ! local variables
        INTEGER(intk) :: i, global_np, read_np
        INTEGER(intk), ALLOCATABLE :: npart_arr(:), ipart_arr(:), p_igrid_arr(:)
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

            ALLOCATE(npart_arr(numprocs))

            CALL dist_npart(numprocs, npart_arr)

            my_particle_list%active_np = npart_arr(myid + 1)

            ALLOCATE(ipart_arr(my_particle_list%active_np))
            ALLOCATE(p_igrid_arr(my_particle_list%active_np))
            ALLOCATE(x(my_particle_list%active_np))
            ALLOCATE(y(my_particle_list%active_np))
            ALLOCATE(z(my_particle_list%active_np))

            CALL dist_ipart(numprocs, myid, npart_arr, ipart_arr)

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

            DEALLOCATE(npart_arr)
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

    SUBROUTINE defragment(this)

        ! SIMON: Here just as an idea...

        ! Subroutine arguments
        CLASS(particle_list_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(intk) :: i, j, ifin, n
        LOGICAL :: cont

        ! local copy
        ifin = this%ifinal
        cont = .TRUE.

        IF ( this%active_np == ifin ) THEN
            ! all slots are occupied
            RETURN
        END IF

        DO i = 1, ifin

            ! search for empty slot
            IF ( this%particles(i)%state < 1 ) THEN

                ! search from the end of list and find particle to fill in
                DO j = this%ifinal, 1, -1
                    ! finished if positions before "i" are considered
                    IF ( j < (i+1) ) THEN
                        cont = .FALSE.
                        EXIT
                    END IF
                    ! fill empty slot with last valid particle
                    IF ( this%particles(j)%state >= 1 ) THEN
                        this%particles(i) = this%particles(j)
                        this%particles(j)%state = -1
                        this%ifinal = j - 1
                        EXIT
                    END IF
                END DO

                ! empty slot could not be filled
                IF ( .NOT. cont ) THEN
                    this%ifinal = i - 1
                    EXIT
                END IF

            END IF

        END DO

        ! Debug check
        DO i = 1, this%ifinal
            IF ( this%particles(i)%state < 1 ) THEN
                WRITE(*,*) "defragmented list contains inactive particle at ", i
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        ! Some output (to be silenced later)
        IF ( ifin /= this%ifinal ) THEN
            n = ifin - this%ifinal
            WRITE(*,*) " - defragmentation: ", n, " slots cleared"
        END IF

    END SUBROUTINE defragment

    !-----------------------------------

    ! this routine sets the number of particles each process is to initialize
    ! (each process does the same computation, but due to the relatively small number of computations
    ! this should be sufficient or even better than a MPI based routine)
    SUBROUTINE dist_npart(nprocs, npart_arr)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: nprocs
        INTEGER(intk), INTENT(inout) :: npart_arr(nprocs)

        ! local variables
        INTEGER(intk) :: i, j, igrid, iproc, iobst, counter
        REAL(realk), ALLOCATABLE :: volume_frac(:)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, dist, volume

        npart_arr = 0

        ALLOCATE(volume_frac(nprocs))
        volume_frac = 0.0

        counter = 1
        DO igrid = 1, ngrid

            iproc = idprocofgrd(igrid)

            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            volume_frac(iproc + 1) = volume_frac(iproc + 1) + (maxx - minx) * (maxy - miny) * (maxz - minz)

            ! ESTIMATE the "free" volume by substracting the obstacle volume of all obstacles on this grid
            ! (NOT EXACT, UNDERESTIMATES THE VOLUME TAKEN BY OBSTACLES)
            IF (dread_obstacles) THEN
                DO i = 1, SIZE(my_obstacle_pointers(igrid)%grid_obstacles)

                    iobst = my_obstacle_pointers(igrid)%grid_obstacles(i)

                    volume_frac(iproc + 1) = volume_frac(iproc + 1) - 3.14 * 4.0 / 3.0 * my_obstacles(iobst)%radius ** 3

                    ! correction for the obstacles that are only partly on the grid
                    ! (overestimates the volume of relevant spherers that is outsie the grid)
                    DO j = 1, 3
                        IF (my_obstacles(iobst)%x - my_obstacles(iobst)%radius < minx) THEN
                            dist = (minx - (my_obstacles(iobst)%x - my_obstacles(iobst)%radius))
                            volume_frac(iproc + 1) = volume_frac(iproc + 1) + &
                            3.14 / 3 * dist**2 * (3 * my_obstacles(iobst)%radius - dist)
                        ELSEIF (my_obstacles(iobst)%x + my_obstacles(iobst)%radius > maxx) THEN
                            dist = ((my_obstacles(iobst)%x + my_obstacles(iobst)%radius) - maxx)
                            volume_frac(iproc + 1) = volume_frac(iproc + 1) + &
                            3.14 / 3 * dist**2 * (3 * my_obstacles(iobst)%radius - dist)
                        END IF

                        IF (my_obstacles(iobst)%y - my_obstacles(iobst)%radius < miny) THEN
                            dist = (miny - (my_obstacles(iobst)%y - my_obstacles(iobst)%radius))
                            volume_frac(iproc + 1) = volume_frac(iproc + 1) + &
                            3.14 / 3 * dist**2 * (3 * my_obstacles(iobst)%radius - dist)
                        ELSEIF (my_obstacles(iobst)%y + my_obstacles(iobst)%radius > maxy) THEN
                            dist = ((my_obstacles(iobst)%y + my_obstacles(iobst)%radius) - maxy)
                            volume_frac(iproc + 1) = volume_frac(iproc + 1) + &
                            3.14 / 3 * dist**2 * (3 * my_obstacles(iobst)%radius - dist)
                        END IF

                        IF (my_obstacles(iobst)%z - my_obstacles(iobst)%radius < minz) THEN
                            dist = (minz - (my_obstacles(iobst)%z - my_obstacles(iobst)%radius))
                            volume_frac(iproc + 1) = volume_frac(iproc + 1) + &
                            3.14 / 3 * dist**2 * (3 * my_obstacles(iobst)%radius - dist)
                        ELSEIF (my_obstacles(iobst)%z + my_obstacles(iobst)%radius > maxz) THEN
                            dist = ((my_obstacles(iobst)%z + my_obstacles(iobst)%radius) - maxz)
                            volume_frac(iproc + 1) = volume_frac(iproc + 1) + &
                            3.14 / 3 * dist**2 * (3 * my_obstacles(iobst)%radius - dist)
                        END IF
                    END DO

                END DO
            END IF

        END DO

        volume = SUM(volume_frac)

        DO iproc = 1, nprocs - 1
            volume_frac(iproc + 1) = volume_frac(iproc + 1) / volume
            npart_arr(iproc + 1) = NINT(init_npart * volume_frac(iproc + 1))
        END DO

        npart_arr(1) = init_npart - SUM(npart_arr)

        do iproc = 1, nprocs
            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) npart_arr(iproc), " particles to be assigned to process ", iproc - 1, "."
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) npart_arr(iproc), " particles to be assigned to process ", iproc - 1, "."
                    WRITE(*, '()')
                END SELECT
            END IF
        END DO

    END SUBROUTINE dist_npart

    !-----------------------------------

    !this routine generates a list of unique particle ids (ipart) on every process
    SUBROUTINE dist_ipart(nprocs, iproc, npart_arr, ipart_arr)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: nprocs, iproc
        INTEGER(intk), INTENT(in) :: npart_arr(nprocs)
        INTEGER(intk), INTENT(out) :: ipart_arr(npart_arr(iproc))

        ! local variables
        INTEGER(intk) :: i, counter, offset

        offset = 0
        DO i = 1, iproc
            offset = offset + npart_arr(i)
        END DO

        counter = 1
        DO i = 1, npart_arr(iproc + 1)
            ipart_arr(counter) = offset + i
            counter = counter + 1
        END DO

    END SUBROUTINE dist_ipart

    !-----------------------------------

    SUBROUTINE dist_particles(npart, p_igrid_arr, x, y, z)

        ! subroutine arguments...
        INTEGER(intk), INTENT(in) :: npart
        INTEGER(intk), INTENT(inout) :: p_igrid_arr(npart)
        REAL(realk), INTENT(inout) :: x(npart), y(npart), z(npart)

        ! local variables...
        INTEGER(intk) :: i, j, igrid, iobst, counter
        REAL(realk) :: myvolume, volume_fractions(nmygrids), grid_rn
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, dist
        LOGICAL :: valid_location

        ! compute combined volume of all grids this process owns (assumning there is only one level!)
        myvolume = 0.0_realk

        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            myvolume = myvolume + (maxx - minx) * (maxy - miny) * (maxz - minz)

        END DO

        ! compute fraction of the combined volume of each grids volume (assumning ther is only one level!)
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
        DO i = 1, npart

            valid_location = .FALSE.
            DO WHILE (.NOT. valid_location)

                CALL RANDOM_NUMBER(grid_rn)

                counter = 1
                DO WHILE (grid_rn > volume_fractions(counter))
                    counter = counter + 1
                END DO

                igrid = mygrids(counter)
                CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

                p_igrid_arr(i) = igrid

                CALL RANDOM_SEED()
                CALL RANDOM_NUMBER(x(i))
                CALL RANDOM_SEED()
                CALL RANDOM_NUMBER(y(i))
                CALL RANDOM_SEED()
                CALL RANDOM_NUMBER(z(i))

                x(i) = minx + x(i) * (maxx - minx)
                y(i) = miny + y(i) * (maxy - miny)
                z(i) = minz + z(i) * (maxz - minz)

                valid_location = .TRUE.

                IF (dread_obstacles) THEN
                    DO j = 1, SIZE(my_obstacle_pointers(igrid)%grid_obstacles)

                        iobst = my_obstacle_pointers(igrid)%grid_obstacles(j)

                        dist = SQRT((my_obstacles(iobst)%x - x(i))**2 + &
                        (my_obstacles(iobst)%y - y(i))**2 + &
                        (my_obstacles(iobst)%z - z(i))**2)

                        IF (dist < my_obstacles(iobst)%radius + 100 * EPSILON(dist)) THEN
                            valid_location = .FALSE.
                            EXIT
                        END IF

                    END DO
                END IF

            END DO

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
