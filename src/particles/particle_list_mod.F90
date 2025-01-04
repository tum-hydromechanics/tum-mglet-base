MODULE particle_list_mod

    ! This module is responsible for:
    ! Definition of the particle_list_t type.
    ! Storage of Particles in my_particle_list.
    ! Basic list manipulations.
    ! Initialization of particles.

    USE MPI_f08
    USE comms_mod

    USE particle_dict_mod
    USE particle_boundaries_mod

    IMPLICIT NONE

    TYPE :: particle_list_t

        ! TODO: remove (obsolete)?
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

    INTEGER(intk) :: global_np

    PUBLIC :: my_particle_list

CONTAINS    !===================================

    SUBROUTINE init_particle_list()

        ! local variables
        INTEGER(intk) :: i, read_np, dummy
        INTEGER(intk), ALLOCATABLE :: ipart_arr(:), igrid_arr(:)
        REAL(realk), ALLOCATABLE :: x_arr(:), y_arr(:), z_arr(:)

        CALL start_timer(900)
        CALL start_timer(910)

        my_particle_list%iproc = myid

        IF (dread_particles_dict) THEN

            ! all arguments that are passed as particle list attributes are input only
            CALL read_particles(dread_particles_dict, ipart_arr, igrid_arr, x_arr, y_arr, z_arr, read_np)

            CALL MPI_Allreduce(read_np, global_np, 1, mglet_mpi_int, MPI_SUM, MPI_COMM_WORLD)

            IF (dread_particles_dict) THEN

                IF (myid == 0) THEN
                    IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                        WRITE(*, '("INITIALIZING ", I0, " PARTICLE(S):")') global_np
                        WRITE(*, '()')
                    END IF
                END IF

                my_particle_list%max_np = SIZE(ipart_arr) ! this is either plist_len (if list_limit is true) or dict_len

                ALLOCATE(my_particle_list%particles(my_particle_list%max_np))

                DO i = 1, read_np

                    CALL set_particle(particle = my_particle_list%particles(i), ipart = ipart_arr(i), &
                     x = x_arr(i), y = y_arr(i), z = z_arr(i), igrid = igrid_arr(i))

                    my_particle_list%active_np = my_particle_list%active_np + 1

                    IF (TRIM(particle_terminal) == "verbose") THEN
                        CALL print_particle_status(my_particle_list%particles(i))
                        WRITE(*, '()')
                    END IF

                END DO

                my_particle_list%ifinal = my_particle_list%active_np

            END IF

            DEALLOCATE(ipart_arr)
            DEALLOCATE(igrid_arr)
            DEALLOCATE(x_arr)
            DEALLOCATE(y_arr)
            DEALLOCATE(z_arr)

        END IF

        IF (.NOT. dread_particles_dict) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, '("INITIALIZING PARTICLE(S): ")')
                    WRITE(*, '()')
                END IF
            END IF

            CALL distribute_particles(ipart_arr, igrid_arr, x_arr, y_arr, z_arr)

            IF (list_limit) THEN
                my_particle_list%max_np = plist_len
            ELSE
                my_particle_list%max_np = INT(1.2 * REAL(SIZE(ipart_arr)))
            END IF

            ALLOCATE(my_particle_list%particles(my_particle_list%max_np))

            DO i = 1, SIZE(ipart_arr)

                CALL set_particle(particle = my_particle_list%particles(i), ipart = ipart_arr(i), &
                 x = x_arr(i), y = y_arr(i), z = z_arr(i), igrid = igrid_arr(i))

                    my_particle_list%active_np = my_particle_list%active_np + 1

                IF (TRIM(particle_terminal) == "verbose") THEN
                    CALL print_particle_status(my_particle_list%particles(i))
                    WRITE(*, '()')
                END IF

            END DO

            CALL MPI_Allreduce(my_particle_list%active_np, global_np, 1, mglet_mpi_int, MPI_SUM, MPI_COMM_WORLD)

            my_particle_list%ifinal = my_particle_list%active_np

            DEALLOCATE(ipart_arr)
            DEALLOCATE(igrid_arr)
            DEALLOCATE(x_arr)
            DEALLOCATE(y_arr)
            DEALLOCATE(z_arr)

        END IF

        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
            IF (myid /= 0) THEN
                CALL MPI_Recv(dummy, 1, mglet_mpi_int, myid - 1, 900, &
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE)
            END IF
            CALL print_list_status(my_particle_list)
            IF (myid /= numprocs - 1) THEN
                CALL MPI_Send(dummy, 1, mglet_mpi_int, myid + 1, 900, &
                 MPI_COMM_WORLD)
            END IF
        END IF

        ! TODO : remove ?
        CALL MPI_Barrier(MPI_COMM_WORLD)

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("INITIALIZATION OF PARTICLE(S) SUCCESSFULLY COMPLETED.")')
                WRITE(*, '()')
            END IF
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

        ! probably unnessecary, but ensures that all new particle slots are inactive
        DO i = particle_list%ifinal + 1, particle_list%max_np + add_len
            particles_tmp(i)%state = -1
        END DO

        particle_list%max_np = particle_list%max_np + add_len
        particle_list%ifinal = MIN(particle_list%ifinal, particle_list%max_np)

        CALL MOVE_ALLOC(particles_tmp, particle_list%particles) ! includes deallocation of particles temp

        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
            WRITE(*, *) "Particle list enlarged by ", add_len, " on process ", myid, "."
            WRITE(*, '()')
        END IF

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

    SUBROUTINE distribute_particles(ipart_arr, igrid_arr, x_arr, y_arr, z_arr)

        ! subroutine arguments
        INTEGER(intk), ALLOCATABLE, INTENT(out) :: ipart_arr(:), igrid_arr(:)
        REAL(realk), ALLOCATABLE, INTENT(out) :: x_arr(:), y_arr(:), z_arr(:)

        ! local variables
        INTEGER(intk), ALLOCATABLE :: init_npart_arr(:), npart_arr(:), itm_igrid_arr(:)!, particles_per_grid_counter(:)
        INTEGER(intk) :: i, j, counter, igrid, iproc, iobst, grid_counter, part_counter, itm_npart, offset
        REAL(realk), ALLOCATABLE :: my_grid_volume_fractions(:), proc_volume_fractions(:)
        REAL(realk), ALLOCATABLE :: itm_x_arr(:), itm_y_arr(:), itm_z_arr(:)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, dist, volume, myvolume, grid_rn
        REAL(realk) :: x, y, z
        LOGICAL :: valid_location

        ALLOCATE(init_npart_arr(0:numprocs-1))
        init_npart_arr = 0

        ALLOCATE(proc_volume_fractions(0:numprocs-1))
        proc_volume_fractions = 0.0

        ! every proc computes the preliminary number of particles of all procs respectively from each proc's combined grid volume
        DO igrid = 1, ngrid

            iproc = idprocofgrd(igrid)

            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            proc_volume_fractions(iproc) = proc_volume_fractions(iproc) + (maxx - minx) * (maxy - miny) * (maxz - minz)

        END DO

        volume = SUM(proc_volume_fractions)

        DO iproc = 1, numprocs - 1
            proc_volume_fractions(iproc) = proc_volume_fractions(iproc) / volume
            ! using init_npart from particle_config_mod
            init_npart_arr(iproc) = NINT(init_npart * proc_volume_fractions(iproc))
        END DO

        init_npart_arr(0) = init_npart - SUM(init_npart_arr)

        ! if list_limit is set true, scale all previously computed entries of init_npart_arr such that
        ! they are smaller or equal to plist_len and their relative size to each other is preserved
        IF (list_limit) THEN
            DO iproc = 0, numprocs - 1
                IF (plist_len < init_npart_arr(iproc)) THEN
                    DO i = 0, numprocs - 1
                        init_npart_arr(i) = FLOOR(REAL(init_npart_arr(i)) * REAL(plist_len) / REAL(init_npart_arr(iproc)))
                    END DO
                END IF
            END DO
        END IF

        ! sanity check
        DO iproc = 0, numprocs -1
            IF (init_npart_arr(iproc) < 0) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        ALLOCATE(itm_igrid_arr(init_npart_arr(myid)))
        ALLOCATE(itm_x_arr(init_npart_arr(myid)))
        ALLOCATE(itm_y_arr(init_npart_arr(myid)))
        ALLOCATE(itm_z_arr(init_npart_arr(myid)))

        !ALLOCATE(particles_per_grid_counter(nmygrids))

        ! compute combined volume of all grids this process owns (assumning there is only one level!)
        myvolume = 0.0_realk

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
            myvolume = myvolume + (maxx - minx) * (maxy - miny) * (maxz - minz)
        END DO

        ! compute each grids (normalized) volume as a fraction of the combined proc volume (volume of all grids this process owns)
        ! (assumning ther is only one level!)
        ALLOCATE(my_grid_volume_fractions(nmygrids))

        igrid = mygrids(1)
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        my_grid_volume_fractions(1) = (maxx - minx) * (maxy - miny) * (maxz - minz) / myvolume

        my_grid_volume_fractions(1) = MAX(0.0_realk, my_grid_volume_fractions(1))
        my_grid_volume_fractions(1) = MIN(1.0_realk, my_grid_volume_fractions(1))

        DO i = 2, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            my_grid_volume_fractions(i) = my_grid_volume_fractions(i-1) + ((maxx - minx) * (maxy - miny) * (maxz - minz) / myvolume)

            my_grid_volume_fractions(i) = MAX(0.0_realk, my_grid_volume_fractions(i))
            my_grid_volume_fractions(i) = MIN(1.0_realk, my_grid_volume_fractions(i))

        END DO
        ! until here, no particles have actually been initialized

        ALLOCATE(npart_arr(0:numprocs-1))
        npart_arr = 0

        ! generate and distribute particles among all grids this process owns proportianally to the volume fraction of each grid
        counter = 1
        part_counter = 1
        ! the while loop ensures that at least 1 particle on any process is initialized here (assuming init_npart > 0)
        DO WHILE(part_counter <= MIN(1, init_npart_arr(myid)) .OR. counter <= init_npart_arr(myid))

            valid_location = .TRUE.

            CALL RANDOM_NUMBER(grid_rn)

            grid_counter = 1
            DO WHILE (grid_rn > my_grid_volume_fractions(grid_counter))
                grid_counter = grid_counter + 1
            END DO

            igrid = mygrids(grid_counter)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            CALL RANDOM_NUMBER(x)
            CALL RANDOM_NUMBER(y)
            CALL RANDOM_NUMBER(z)

            x = minx + x * (maxx - minx)
            y = miny + y * (maxy - miny)
            z = minz + z * (maxz - minz)

            IF (dread_obstacles) THEN
                DO j = 1, SIZE(my_obstacle_pointers(igrid)%grid_obstacles)

                    iobst = my_obstacle_pointers(igrid)%grid_obstacles(j)

                    dist = SQRT((my_obstacles(iobst)%x - x)**2 + &
                    (my_obstacles(iobst)%y - y)**2 + &
                    (my_obstacles(iobst)%z - z)**2)

                    IF (dist < (my_obstacles(iobst)%radius + aura)) THEN
                        valid_location = .FALSE.
                        EXIT
                    END IF

                END DO
            END IF

            IF (valid_location) THEN
                !particles_per_grid_counter(grid_counter) = particles_per_grid_counter(grid_counter) + 1
                itm_igrid_arr(part_counter) = igrid
                itm_x_arr(part_counter) = x
                itm_y_arr(part_counter) = y
                itm_z_arr(part_counter) = z
                part_counter = part_counter + 1
            END IF

            counter = counter + 1

        END DO

        ! gather number of particles that have been initialized at a valid location on each process
        CALL MPI_Gather(part_counter - 1, 1, mglet_mpi_int, &
         npart_arr, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)

        ! compute the final number of particles that is to be initialized on each process respectively
        IF (myid == 0) THEN
            itm_npart = SUM(npart_arr)

            IF (itm_npart < 1) CALL errr(__FILE__, __LINE__)

            DO iproc = 1, numprocs - 1
                npart_arr(iproc) = INT(npart_arr(iproc) * init_npart / itm_npart, intk)
            END DO

            npart_arr(0) = 0
            npart_arr(0) = init_npart - SUM(npart_arr)

            ! if list_limit is set true, scale all previously computed entries of init_npart_arr such that
            ! they are smaller or equal to plist_len and their relative size to each other is preserved (once more)
            IF (list_limit) THEN
                DO iproc = 0, numprocs - 1
                    IF (plist_len < npart_arr(iproc)) THEN
                        DO i = 0, numprocs - 1
                            npart_arr(i) = FLOOR(REAL(npart_arr(i)) * REAL(plist_len) / REAL(npart_arr(iproc)))
                        END DO
                    END IF
                END DO
                IF (SUM(npart_arr) /= init_npart) THEN
                    IF (TRIM(particle_terminal) == "verbose") THEN
                        WRITE(*, *) "The specified number of Particles could not be initialized due to Particle List Limitations."
                        WRITE(*, '()')
                    END IF
                END IF
            ELSE
                IF (SUM(npart_arr) /= init_npart) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
            END IF
        END IF

        ! TODO: remove Bcast?
        CALL MPI_Bcast(npart_arr, numprocs, mglet_mpi_int, &
         0, MPI_COMM_WORLD)

        ! move already generated positions to the final output array
        ALLOCATE(ipart_arr(npart_arr(myid)))
        ALLOCATE(igrid_arr(npart_arr(myid)))
        ALLOCATE(x_arr(npart_arr(myid)))
        ALLOCATE(y_arr(npart_arr(myid)))
        ALLOCATE(z_arr(npart_arr(myid)))

        DO i = 1, part_counter - 1
            igrid_arr(i) = itm_igrid_arr(i)
            x_arr(i) = itm_x_arr(i)
            y_arr(i) = itm_y_arr(i)
            z_arr(i) = itm_z_arr(i)
        END DO

        ! distribute particles among all grids this process owns (uniformely distibuted)
        ! to eachieve the predefined number of particles
        DO i = part_counter, npart_arr(myid)

            valid_location = .FALSE.
            DO WHILE (.NOT. valid_location)

                CALL RANDOM_NUMBER(grid_rn)

                grid_counter = 1
                DO WHILE (grid_rn > my_grid_volume_fractions(grid_counter))
                    grid_counter = grid_counter + 1
                END DO

                igrid = mygrids(grid_counter)
                CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

                CALL RANDOM_NUMBER(x)
                CALL RANDOM_NUMBER(y)
                CALL RANDOM_NUMBER(z)

                x = minx + x * (maxx - minx)
                y = miny + y * (maxy - miny)
                z = minz + z * (maxz - minz)

                valid_location = .TRUE.

                IF (dread_obstacles) THEN
                    DO j = 1, SIZE(my_obstacle_pointers(igrid)%grid_obstacles)

                        iobst = my_obstacle_pointers(igrid)%grid_obstacles(j)

                        dist = SQRT((my_obstacles(iobst)%x - x)**2 + &
                        (my_obstacles(iobst)%y - y)**2 + &
                        (my_obstacles(iobst)%z - z)**2)

                        IF (dist < (my_obstacles(iobst)%radius + aura)) THEN
                            valid_location = .FALSE.
                            EXIT
                        END IF

                    END DO
                END IF

                IF (valid_location) THEN
                    igrid_arr(i) = igrid
                    x_arr(i) = x
                    y_arr(i) = y
                    z_arr(i) = z
                END IF

            END DO

        END DO

        offset = 0
        DO i = 0, myid - 1
            offset = offset + npart_arr(i)
        END DO

        DO i = 1, npart_arr(myid)
            ipart_arr(i) = offset + i
        END DO

        DEALLOCATE(itm_igrid_arr)
        DEALLOCATE(itm_x_arr)
        DEALLOCATE(itm_y_arr)
        DEALLOCATE(itm_z_arr)
        DEALLOCATE(my_grid_volume_fractions)
        DEALLOCATE(proc_volume_fractions)
        DEALLOCATE(init_npart_arr)
        DEALLOCATE(npart_arr)

    END SUBROUTINE distribute_particles

    !-----------------------------------

    SUBROUTINE print_list_status(particle_list)

        ! subroutine arguments
        TYPE(particle_list_t), INTENT(in) :: particle_list

        WRITE(*, *) "Particle list status on process ", myid, ": max_np = ", &
         particle_list%max_np, ", active_np = ", particle_list%active_np, ", ifinal = ", particle_list%ifinal
        WRITE(*, '()')

    END SUBROUTINE print_list_status

    !-----------------------------------

    ! for debugging
    SUBROUTINE write_particle_list_txt(itstep, suffix)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        CHARACTER(len = 3), INTENT(in), OPTIONAL :: suffix

        ! local varibales
        CHARACTER(len = mglet_filename_max) :: filename
        INTEGER(intk) :: unit, i
        LOGICAL :: exists

        IF (PRESENT(suffix)) THEN
            WRITE(filename,'("ParticleList-", I0, "-", A, ".txt")') my_particle_list%iproc, suffix
        ELSE
            WRITE(filename,'("ParticleList-", I0, ".txt")') my_particle_list%iproc
        END IF

        INQUIRE(file = TRIM(filename), exist = exists)

        IF (exists) THEN
            OPEN(newunit = unit, file = TRIM(filename), status = 'OLD', action = 'WRITE')
        ELSE
            OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')
        END IF

        WRITE(unit, '("PARTICLE LIST ", I0, " - TIMESTEP ", I0)') my_particle_list%iproc, itstep
        WRITE(unit, '(" ")')
        WRITE(unit, '("active_np: ", I0)') my_particle_list%active_np
        WRITE(unit, '("ifinal: ", I0)') my_particle_list%ifinal
        WRITE(unit, '(" ")')
        WRITE(unit, '("PARTICLES")')

        DO i = 1, SIZE(my_particle_list%particles)

                WRITE(unit, '("i = ", I0)') i
                WRITE(unit, '("ipart = ", I9, ", iproc = ", I3, ", igrid = ", I3, ", state = ", I3)') my_particle_list%particles(i)%ipart, &
                 my_particle_list%particles(i)%iproc, my_particle_list%particles(i)%igrid, my_particle_list%particles(i)%state
                WRITE(unit, '("i/j/k cell :", 3I9)') my_particle_list%particles(i)%ijkcell(1), &
                 my_particle_list%particles(i)%ijkcell(2), my_particle_list%particles(i)%ijkcell(3)
                WRITE(unit, '("x/y/z      :", 3F9.6)') my_particle_list%particles(i)%x, &
                 my_particle_list%particles(i)%y, my_particle_list%particles(i)%z

        END DO

        CLOSE(unit)

    END SUBROUTINE write_particle_list_txt

    SUBROUTINE finish_particle_list()

        DEALLOCATE(my_particle_list%particles)

    END SUBROUTINE finish_particle_list

END MODULE
