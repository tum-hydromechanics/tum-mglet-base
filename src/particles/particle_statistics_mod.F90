MODULE particle_statistics_mod

    USE utils_mod
    USE grids_mod

    USE particle_list_mod

    IMPLICIT NONE

    ! TODO: renaming, merging of collector types into abstract type and definition of type extendions for differentiation

    TYPE :: gridstat_collector_t

        INTEGER(intk) :: igrid

        ! volume of the grid
        REAL(realk) :: grid_volume

        ! array that stores the number of particles on the grid for each timestep
        INTEGER(intk), ALLOCATABLE :: np_counter(:)

        ! array that stores the residence time distribution
        INTEGER(intk), ALLOCATABLE :: rt_counter(:)

    END TYPE gridstat_collector_t

    TYPE :: slicestat_collector_t

        ! indicates if the slice is relevant for the process by which it is owned (true) or not (false)
        LOGICAL :: is_active = .FALSE.

        ! slice bounds (llim = lower limit, ulim = upper limit)
        REAL(realk) :: llim
        REAL(realk) :: ulim

        ! array that stores the number of particles on the grid for each timestep
        INTEGER(intk), ALLOCATABLE :: np_counter(:)

        ! array that stores the residence time distribution
        INTEGER(intk), ALLOCATABLE :: rt_counter(:)

    END TYPE slicestat_collector_t

    TYPE(gridstat_collector_t), ALLOCATABLE :: my_gridcol_list(:)
    TYPE(slicestat_collector_t), ALLOCATABLE :: my_slicecol_list(:)

CONTAINS

    SUBROUTINE init_particle_statistics(mtstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: mtstep

        CALL start_timer(900)
        CALL start_timer(950)

        CALL init_particle_gridstat(mtstep)

        CALL init_particle_slicestat(mtstep)

        CALL stop_timer(950)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_statistics

    SUBROUTINE init_particle_gridstat(mtstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: mtstep

        ! local variables
        INTEGER(intk) :: igrid, i
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        !LOGICAL :: gridstat_exists

        ! TODO: inquire for directories ?
        !INQUIRE(directory = './Particle_Statistics', exist = gridstat_exists)
!
        !IF (gridstat_exists) THEN
!
        !    IF (myid == 0) THEN
        !        SELECT CASE (TRIM(particle_terminal))
        !            CASE ("none")
        !                CONTINUE
        !            CASE ("normal")
        !                WRITE(*, *) ' '
        !                WRITE(*, *) "ERROR: Directory Particle_Statistics already exists. Terminating Process!"
        !            CASE ("verbose")
        !                WRITE(*, *) ' '
        !                WRITE(*, *) "ERROR: Directory Particle_Statistics already exists. Terminating Process!"
        !        END SELECT
        !    END IF
!
        !    CALL errr(__FILE__, __LINE__)
!
        !END IF

        IF (.NOT. dgridstat) THEN
            RETURN
        END IF

        ! allocate collector list
        ALLOCATE(my_gridcol_list(nmygrids))

        ! allocate arrays of the individual collectors
        DO i = 1, nmygrids

            igrid = mygrids(i)

            my_gridcol_list(i)%igrid = igrid

            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            my_gridcol_list(i)%grid_volume = (maxx - minx) * (maxy -miny) * (maxz - minz)

            ALLOCATE(my_gridcol_list(i)%np_counter(mtstep + 1))
            my_gridcol_list(i)%np_counter = 0

            ALLOCATE(my_gridcol_list(i)%rt_counter(rt_tstep_max))
            my_gridcol_list(i)%rt_counter = 0

        END DO

        CALL start_np_counter()

    END SUBROUTINE init_particle_gridstat

    SUBROUTINE start_np_counter()

        ! local variables
        INTEGER(intk) :: igrid, i, j

        IF (.NOT. dgridstat) THEN
            RETURN
        END IF

        DO i = 1, nmygrids
            igrid = mygrids(i)
            ! initialize first entry of np_counter (count the number of particles in each grid at itstep = 0 resprectively)
            DO j = 1, my_particle_list%ifinal
                IF (my_particle_list%particles(i)%gitstep >= 0 .AND. &
                 my_particle_list%particles(i)%gitstep < rt_ittot_start) THEN
                    CALL errr(__FILE__,__LINE__)
                END IF
                IF (my_particle_list%particles(j)%igrid == my_gridcol_list(i)%igrid) THEN
                    my_gridcol_list(i)%np_counter(1) = my_gridcol_list(i)%np_counter(1) + 1
                END IF
            END DO
        END DO

    END SUBROUTINE start_np_counter

    SUBROUTINE advance_np_counter(itstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep

        ! local variables
        INTEGER(intk) :: i

        IF (.NOT. dgridstat) THEN
            RETURN
        END IF

        CALL start_timer(950)

        DO i = 1, nmygrids
            my_gridcol_list(i)%np_counter(itstep + 1) = my_gridcol_list(i)%np_counter(itstep)
        END DO

        CALL stop_timer(950)

    END SUBROUTINE advance_np_counter

    SUBROUTINE deregister_particle(particle, ittot, itstep)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: itstep

        ! local variables
        INTEGER(intk) :: igrid, i, rt_tstep

        IF (.NOT. dgridstat) THEN
            RETURN
        END IF

        CALL start_timer(950)

        DO i = 1, nmygrids

            igrid = mygrids(i)

            IF (particle%igrid == igrid) THEN

                my_gridcol_list(i)%np_counter(itstep + 1) = my_gridcol_list(i)%np_counter(itstep + 1) - 1

                IF (ittot >= rt_ittot_start) THEN

                    IF (particle%gitstep >= 0) THEN
                        rt_tstep = ittot - particle%gitstep

                        IF (rt_tstep > SIZE(my_gridcol_list(i)%rt_counter) .OR. rt_tstep < 1) THEN
                            WRITE(*, '("WARNING in Particle Statistics (Grids)")')
                            WRITE(*, '("On Proc ", I0, ": Residence time of particle ", I0," is too large or too small!")') &
                             myid, particle%ipart
                            WRITE(*, '("rt_tstep ", I0, "; size(rt_counter)", I0)') rt_tstep, SIZE(my_gridcol_list(i)%rt_counter)
                        ELSE
                            my_gridcol_list(i)%rt_counter(rt_tstep) = my_gridcol_list(i)%rt_counter(rt_tstep) + 1
                        END IF
                    END IF

                    particle%gitstep = ittot
                    EXIT

                END IF

            END IF

        END DO

        CALL stop_timer(950)

    END SUBROUTINE deregister_particle

    SUBROUTINE register_particle(particle, itstep)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: itstep

        ! local variables
        INTEGER(intk) :: igrid, i

        IF (.NOT. dgridstat) THEN
            RETURN
        END IF

        CALL start_timer(950)

        DO i = 1, nmygrids

            igrid = mygrids(i)

            IF (particle%igrid == igrid) THEN
                my_gridcol_list(i)%np_counter(itstep + 1) = my_gridcol_list(i)%np_counter(itstep + 1) + 1
                EXIT
            END IF

        END DO

        CALL stop_timer(950)

    END SUBROUTINE register_particle

    SUBROUTINE init_particle_slicestat(mtstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: mtstep

        ! local variables
        INTEGER(intk) :: igrid, i, j, k, counter, dummy
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, grid_min, grid_max
        REAL(realk) :: global_x0, global_x1, global_y0, global_y1, global_z0, global_z1
        REAL(realk) :: global_min, global_max, lim
        !TYPE(MPI_Request) :: request

        IF (slice_dir == "N") THEN
            ALLOCATE(my_slicecol_list(0))
            RETURN
        END IF

        CALL get_bbox(global_x0, global_x1, global_y0, global_y1, global_z0, global_z1, igrid = igrdoflevel(1, particle_level))

        DO i = 2, noflevel(particle_level)

            igrid = igrdoflevel(i, particle_level)

            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            IF (minx < global_x0) global_x0 = minx
            IF (maxx > global_x1) global_x1 = maxx
            IF (miny < global_y0) global_y0 = miny
            IF (maxy > global_y1) global_y1 = maxy
            IF (minz < global_z0) global_z0 = minz
            IF (maxz > global_z1) global_z1 = maxz

        END DO

        SELECT CASE(slice_dir)
            CASE("X")
                    global_min = global_x0
                    global_max = global_x1
            CASE("Y")
                    global_min = global_y0
                    global_max = global_y1
            CASE("Z")
                    global_min = global_z0
                    global_max = global_z1
        END SELECT

        ALLOCATE(my_slicecol_list(SUM(nslices)))

        ! determine slice limits
        counter = 1
        lim = global_min
        DO i = 1, SIZE(nslices)
            DO j = 1, nslices(i)

                my_slicecol_list(counter)%llim = lim
                lim = lim + (global_max - global_min) * slice_levels(i) / nslices(i)

                IF (counter == SUM(nslices)) THEN
                    my_slicecol_list(counter)%ulim = global_max
                ELSE
                    my_slicecol_list(counter)%ulim = lim
                END IF

                ! check if this slice does overlap with any grid on this process
                DO k = 1, nmygridslvl(particle_level)

                    igrid = mygridslvl(k, particle_level)
                    CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

                    SELECT CASE(slice_dir)
                        CASE("X")
                                grid_min = minx
                                grid_max = maxx
                        CASE("Y")
                                grid_min = miny
                                grid_max = maxy
                        CASE("Z")
                                grid_min = minz
                                grid_max = maxz
                    END SELECT

                    IF (grid_min > my_slicecol_list(counter)%ulim .OR. grid_max < my_slicecol_list(counter)%llim) THEN
                        CYCLE
                    ELSE
                        my_slicecol_list(counter)%is_active = .TRUE.
                        ! neighbouring slices that do not overlap with any grid also have to be active for residence time tracking!
                        IF (counter == 1) THEN
                            my_slicecol_list(SUM(nslices))%is_active = .TRUE.
                            my_slicecol_list(counter + 1)%is_active = .TRUE.
                        ELSEIF (counter == SUM(nslices)) THEN
                            my_slicecol_list(counter - 1)%is_active = .TRUE.
                            my_slicecol_list(1)%is_active = .TRUE.
                        ELSE
                            my_slicecol_list(counter - 1)%is_active = .TRUE.
                            my_slicecol_list(counter + 1)%is_active = .TRUE.
                        END IF

                        EXIT
                    END IF
                END DO

                counter = counter + 1

            END DO
        END DO

        DO i = 1, SIZE(my_slicecol_list)
            IF (my_slicecol_list(i)%is_active) THEN
                ALLOCATE(my_slicecol_list(i)%np_counter(mtstep + 1))
                my_slicecol_list(i)%np_counter = 0
                ALLOCATE(my_slicecol_list(i)%rt_counter(rt_tstep_max))
                my_slicecol_list(i)%rt_counter = 0
            ELSE
                IF (myid /= 0) THEN
                    ALLOCATE(my_slicecol_list(i)%np_counter(0))
                    ALLOCATE(my_slicecol_list(i)%rt_counter(0))
                END IF
            END IF
        END DO

        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN

            IF (myid /= 0) THEN
                CALL MPI_Recv(dummy, 1, mglet_mpi_int, myid - 1, 900, &
                MPI_COMM_WORLD, MPI_STATUS_IGNORE)
            END IF

            WRITE(*, '("Slice Statistics Collector of Length ", I0, " allocated on Proccess ", I0)') SIZE(my_slicecol_list), myid
            WRITE(*, '()')

            IF (myid /= numprocs - 1) THEN
                CALL MPI_Send(dummy, 1, mglet_mpi_int, myid + 1, 900, &
                MPI_COMM_WORLD)
            END IF

        END IF

        CALL assign_initial_slice(my_particle_list)

        CALL MPI_Barrier(MPI_COMM_WORLD)

    END SUBROUTINE init_particle_slicestat

    SUBROUTINE assign_initial_slice(particle_list)

        ! subroutine arguments
        TYPE(particle_list_t), INTENT(inout) :: particle_list

        ! local variables
        INTEGER(intk) :: i, islice
        REAL(realk) :: p_coord

        IF (slice_dir == "N") THEN
            RETURN
        END IF

        DO i = 1, particle_list%ifinal
            IF (particle_list%particles(i)%sitstep >= 0 .AND. &
             particle_list%particles(i)%sitstep < rt_ittot_start) THEN
                CALL errr(__FILE__,__LINE__)
            END IF
            SELECT CASE(slice_dir)
                CASE("X")
                        p_coord = particle_list%particles(i)%x
                CASE("Y")
                        p_coord = particle_list%particles(i)%y
                CASE("Z")
                        p_coord = particle_list%particles(i)%z
            END SELECT
            sliceloop: DO islice = 1, SIZE(my_slicecol_list)
                IF (my_slicecol_list(islice)%llim <= p_coord .AND. p_coord <= my_slicecol_list(islice)%ulim) THEN
                    IF (particle_list%particles(i)%islice <= 0) THEN
                        particle_list%particles(i)%islice = islice
                    ELSEIF (particle_list%particles(i)%islice /= islice) THEN
                        WRITE(*, '("ERROR on proc ", I0, ": Unexpected initial value of islice for particle ", I0)') &
                             myid, particle_list%particles(i)%ipart
                        CALL errr(__FILE__,__LINE__)
                    END IF
                    my_slicecol_list(islice)%np_counter(1) = my_slicecol_list(islice)%np_counter(1) + 1
                    EXIT sliceloop
                END IF
            END DO sliceloop
        END DO

    END SUBROUTINE assign_initial_slice

    SUBROUTINE associate_new_slice(particle, ittot, itstep)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: itstep

        ! local variables
        INTEGER(intk) :: rt_tstep
        REAL(realk) :: p_coord

        IF (slice_dir == "N") THEN
            RETURN
        END IF

        CALL start_timer(950)

        SELECT CASE(slice_dir)
            CASE("X")
                    p_coord = particle%x
            CASE("Y")
                    p_coord = particle%y
            CASE("Z")
                    p_coord = particle%z
        END SELECT

        IF (p_coord < my_slicecol_list(particle%islice)%llim) THEN

            IF (particle%islice == 1) THEN
                particle%islice = SIZE(my_slicecol_list)
            ELSE
                particle%islice = particle%islice - 1
            END IF

            IF (ittot >= rt_ittot_start) THEN
                IF (particle%sitstep >= 0) THEN
                    rt_tstep = ittot - particle%sitstep
                    IF (rt_tstep > SIZE(my_slicecol_list(particle%islice)%rt_counter) .OR. rt_tstep < 1) THEN
                        WRITE(*, '("WARNING in Particle Statistics (Slices)")')
                        WRITE(*, '("On Proc ", I0, ": Residence time of particle ", I0," is too large or too small!")') &
                         myid, particle%ipart
                        WRITE(*, '("rt_tstep ", I0, "; size(rt_counter)", I0)') rt_tstep, SIZE(my_slicecol_list(particle%islice)%rt_counter)
                    ELSE
                        my_slicecol_list(particle%islice)%rt_counter(rt_tstep) = &
                        my_slicecol_list(particle%islice)%rt_counter(rt_tstep) + 1
                    END IF
                END IF

                particle%sitstep = ittot
            END IF

        ELSEIF (p_coord > my_slicecol_list(particle%islice)%ulim) THEN

            IF (particle%islice == SIZE(my_slicecol_list)) THEN
                particle%islice = 1
            ELSE
                particle%islice = particle%islice + 1
            END IF

            IF (ittot >= rt_ittot_start) THEN
                IF (particle%sitstep >= 0) THEN
                    rt_tstep = ittot - particle%sitstep
                    IF (rt_tstep > SIZE(my_slicecol_list(particle%islice)%rt_counter) .OR. rt_tstep < 1) THEN
                        WRITE(*, '("WARNING in Particle Statistics (Slices)")')
                        WRITE(*, '("On Proc ", I0, ": Residence time of particle ", I0," is too large or too small!")') &
                         myid, particle%ipart
                        WRITE(*, '("rt_tstep ", I0, "; size(rt_counter)", I0)') rt_tstep, SIZE(my_slicecol_list(particle%islice)%rt_counter)
                    ELSE
                        my_slicecol_list(particle%islice)%rt_counter(rt_tstep) = &
                        my_slicecol_list(particle%islice)%rt_counter(rt_tstep) + 1
                    END IF
                END IF

                particle%sitstep = ittot
            END IF

        END IF

        my_slicecol_list(particle%islice)%np_counter(itstep + 1) = &
         my_slicecol_list(particle%islice)%np_counter(itstep + 1) + 1

        CALL stop_timer(950)

    END SUBROUTINE associate_new_slice

    SUBROUTINE merge_slicestat(islice, mtstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: islice
        INTEGER(intk), INTENT(in) :: mtstep

        ! local variables
        INTEGER(intk) :: i, j, counter, buffer_len, dsend
        INTEGER(intk), ALLOCATABLE :: drecv(:)
        INTEGER(intk), ALLOCATABLE :: buffer(:)
        TYPE(MPI_Request) :: send_req
        TYPE(MPI_Request), ALLOCATABLE :: recv_req(:)

        !TODO: optimize this dirty routine

        IF (slice_dir == "N") THEN
            RETURN
        END IF

        IF (myid == 0 .AND. .NOT. my_slicecol_list(islice)%is_active) THEN
            ALLOCATE(my_slicecol_list(islice)%np_counter(mtstep + 1))
            my_slicecol_list(islice)%np_counter = 0
            ALLOCATE(my_slicecol_list(islice)%rt_counter(rt_tstep_max))
            my_slicecol_list(islice)%rt_counter = 0
        END IF

        ALLOCATE(recv_req(numprocs -1))

        ALLOCATE(drecv(numprocs - 1))
        drecv = 0
        dsend = 0

        IF (my_slicecol_list(islice)%is_active) THEN
            dsend = 1
        END IF

        IF (myid /= 0) THEN
            CALL MPI_Isend(dsend, 1, mglet_mpi_int, 0, 9401, &
             MPI_COMM_WORLD, send_req)
            CALL MPI_Wait(send_req, MPI_STATUS_IGNORE)
        ELSEIF (myid == 0) THEN
            DO i = 1, numprocs - 1
                CALL MPI_Irecv(drecv(i), 1, mglet_mpi_int, i, 9401, &
                     MPI_COMM_WORLD, recv_req(i))
            END DO
            CALL MPI_Waitall(numprocs - 1, recv_req, MPI_STATUSES_IGNORE)
        END IF


        buffer_len = mtstep + 1 + rt_tstep_max
        ALLOCATE(buffer(buffer_len))
        buffer = 0

        IF (myid /= 0) THEN

            IF (dsend == 1) THEN
                ! fill buffer
                counter = 1
                DO j = 1, mtstep + 1
                    buffer(j) = my_slicecol_list(islice)%np_counter(counter)
                    counter = counter + 1
                END DO

                counter = 1
                DO j = mtstep + 2, buffer_len
                    buffer(j) = my_slicecol_list(islice)%rt_counter(counter)
                    counter = counter + 1
                END DO

                CALL MPI_Send(buffer, buffer_len, mglet_mpi_int, 0, 9402, &
                     MPI_COMM_WORLD)

                !CALL MPI_Wait(send_req, MPI_STATUS_IGNORE)

            END IF

        ELSEIF (myid == 0) THEN

            DO i = 1, numprocs - 1

                IF (drecv(i) == 1) THEN

                    buffer = 0

                    CALL MPI_Recv(buffer, buffer_len, mglet_mpi_int, i, 9402, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE)

                    !CALL MPI_Wait(recv_req(i), MPI_STATUS_IGNORE)

                    counter = 1
                    DO j = 1, mtstep + 1
                        my_slicecol_list(islice)%np_counter(counter) = my_slicecol_list(islice)%np_counter(counter) + buffer(j)
                        counter = counter + 1
                    END DO

                    counter = 1
                    DO j = mtstep + 2, buffer_len
                        my_slicecol_list(islice)%rt_counter(counter) = my_slicecol_list(islice)%rt_counter(counter) + buffer(j)
                        counter = counter + 1
                    END DO

                ELSE

                    CYCLE

                END IF

            END DO

        END IF

        DEALLOCATE(buffer)
        DEALLOCATE(recv_req)
        DEALLOCATE(drecv)

    END SUBROUTINE merge_slicestat

    SUBROUTINE deallocate_slicestat(islice)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: islice

        IF (.NOT. ALLOCATED(my_slicecol_list)) RETURN

        IF (islice < 1 .OR. islice > SIZE(my_slicecol_list)) RETURN

        IF (ALLOCATED(my_slicecol_list(islice)%np_counter)) DEALLOCATE(my_slicecol_list(islice)%np_counter)
        IF (ALLOCATED(my_slicecol_list(islice)%rt_counter)) DEALLOCATE(my_slicecol_list(islice)%rt_counter)

    END SUBROUTINE deallocate_slicestat

    SUBROUTINE deallocate_gridstat(i)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: i

        IF (.NOT. ALLOCATED(my_gridcol_list)) RETURN

        IF (i < 1 .OR. i > SIZE(my_gridcol_list)) RETURN

        IF (ALLOCATED(my_gridcol_list(i)%np_counter)) DEALLOCATE(my_gridcol_list(i)%np_counter)
        IF (ALLOCATED(my_gridcol_list(i)%rt_counter)) DEALLOCATE(my_gridcol_list(i)%rt_counter)

    END SUBROUTINE deallocate_gridstat

    SUBROUTINE write_particle_statistics(mtstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: mtstep

        ! local variables
        INTEGER(intk) :: islice

        CALL start_timer(900)
        CALL start_timer(950)

        CALL write_gridstat_folder()

        CALL write_gridstat_files()

        DO islice = 1, SIZE(my_slicecol_list)
            CALL merge_slicestat(islice, mtstep)
            CALL write_slicestat_file(islice)
        END DO

        CALL stop_timer(950)
        CALL stop_timer(900)

    END SUBROUTINE write_particle_statistics

    SUBROUTINE write_gridstat_folder()

        IF (myid == 0) THEN
            CALL create_directory("Particle_Statistics") ! ! ! realtive to working directory ! ! !
        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

    END SUBROUTINE write_gridstat_folder

    SUBROUTINE write_gridstat_files()

        INTEGER(intk) :: i, j, unit
        CHARACTER(len = mglet_filename_max) :: filename

        DO i = 1, nmygrids

            WRITE(filename,'("Particle_Statistics/Grid", I0, ".txt")') my_gridcol_list(i)%igrid

            OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

            WRITE(unit, '("IGRID: ", I0)') my_gridcol_list(i)%igrid

            WRITE(unit, '("GRID VOLUME: ", F12.6)') my_gridcol_list(i)%grid_volume

            WRITE(unit, '("NPART TIMELINE:")')

            DO j = 1, SIZE(my_gridcol_list(i)%np_counter)

                WRITE(unit, '(I0)', advance = "no") my_gridcol_list(i)%np_counter(j)

                IF (j /= SIZE(my_gridcol_list(i)%np_counter)) THEN
                    WRITE(unit, '(", ")', advance = "no")
                    IF (MOD(j, 10) == 0) THEN
                        WRITE(unit, '("")', advance = "yes")
                    END IF
                END IF

            END DO

            WRITE(unit, '("")', advance = "yes")

            WRITE(unit, '("RT DISTRIBUTION:")')

            DO j = 1, SIZE(my_gridcol_list(i)%rt_counter)

                WRITE(unit, '(I0)', advance = "no") my_gridcol_list(i)%rt_counter(j)

                IF (j /= SIZE(my_gridcol_list(i)%rt_counter)) THEN
                    WRITE(unit, '(", ")', advance = "no")
                    IF (MOD(j, 10) == 0) THEN
                        WRITE(unit, '("")', advance = "yes")
                    END IF
                END IF

            END DO

            CLOSE(unit)

        END DO

    END SUBROUTINE write_gridstat_files

    SUBROUTINE write_slicestat_file(islice)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: islice

        ! local varibales
        INTEGER(intk) :: j, unit
        CHARACTER(len = mglet_filename_max) :: filename

        IF (myid == 0) THEN

            WRITE(filename,'("Particle_Statistics/Slice", I0, ".txt")') islice

            OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

            WRITE(unit, '("ISLICE: ", I0)') islice

            WRITE(unit, '("SLICE LIMITS: ", 2F20.17)') my_slicecol_list(islice)%llim, my_slicecol_list(islice)%ulim

            WRITE(unit, '("NPART TIMELINE:")')

            DO j = 1, SIZE(my_slicecol_list(islice)%np_counter)
                WRITE(unit, '(I0)', advance = "no") my_slicecol_list(islice)%np_counter(j)
                IF (j /= SIZE(my_slicecol_list(islice)%np_counter)) THEN
                    WRITE(unit, '(", ")', advance = "no")
                    IF (MOD(j, 10) == 0) THEN
                        WRITE(unit, '("")', advance = "yes")
                    END IF
                END IF
            END DO

            WRITE(unit, '("")', advance = "yes")

            WRITE(unit, '("RT DISTRIBUTION:")')

            DO j = 1, SIZE(my_slicecol_list(islice)%rt_counter)
                WRITE(unit, '(I0)', advance = "no") my_slicecol_list(islice)%rt_counter(j)

                IF (j /= SIZE(my_slicecol_list(islice)%rt_counter)) THEN
                    WRITE(unit, '(", ")', advance = "no")
                    IF (MOD(j, 10) == 0) THEN
                        WRITE(unit, '("")', advance = "yes")
                    END IF
                END IF
            END DO

            CLOSE(unit)

        END IF

    END SUBROUTINE write_slicestat_file


    SUBROUTINE finish_particle_statistics()

        ! local varibales
        INTEGER(intk) :: i

        CALL start_timer(900)
        CALL start_timer(950)

        IF (ALLOCATED(my_gridcol_list)) THEN
            DO i = 1, SIZE(my_gridcol_list)
                CALL deallocate_gridstat(i)
            END DO
            DEALLOCATE(my_gridcol_list)
        END IF

        IF (ALLOCATED(my_slicecol_list)) THEN
            DO i = 1, SIZE(my_slicecol_list)
                CALL deallocate_slicestat(i)
            END DO
            DEALLOCATE(my_slicecol_list)
        END IF

        CALL stop_timer(950)
        CALL stop_timer(900)

    END SUBROUTINE finish_particle_statistics

END MODULE particle_statistics_mod