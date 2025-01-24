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

        ! arrays that stores the travel distance distributions
        INTEGER(intk), ALLOCATABLE :: tdx_counter(:)
        INTEGER(intk), ALLOCATABLE :: tdy_counter(:)
        INTEGER(intk), ALLOCATABLE :: tdz_counter(:)

        ! arrays that stores the average travel velocity distributions
        INTEGER(intk), ALLOCATABLE :: avu_counter(:)
        INTEGER(intk), ALLOCATABLE :: avv_counter(:)
        INTEGER(intk), ALLOCATABLE :: avw_counter(:)

    END TYPE slicestat_collector_t

    REAL(realk) :: global_min_ddx, global_min_ddy, global_min_ddz

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
                IF (my_particle_list%particles(j)%gitstep >= 0 .AND. &
                 my_particle_list%particles(j)%gitstep < rt_ittot_start) THEN
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
        REAL(realk) :: global_min, global_max, min_ddx, min_ddy, min_ddz, lim
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx, ddy, ddz
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f

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

        ! probably unneccessary
        CALL MPI_Barrier(MPI_COMM_WORLD)

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        min_ddx = global_x1 - global_x0
        min_ddy = global_x1 - global_x0
        min_ddz = global_x1 - global_x0

        DO i = 1, nmygridslvl(particle_level)
            igrid = mygridslvl(i, particle_level)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)
            min_ddx = MIN(MINVAL(ddx), min_ddx)
            min_ddy = MIN(MINVAL(ddy), min_ddy)
            min_ddz = MIN(MINVAL(ddz), min_ddz)
        END DO

        CALL MPI_Allreduce(min_ddx, global_min_ddx, 1, mglet_mpi_real, MPI_MIN, MPI_COMM_WORLD)
        CALL MPI_Allreduce(min_ddy, global_min_ddy, 1, mglet_mpi_real, MPI_MIN, MPI_COMM_WORLD)
        CALL MPI_Allreduce(min_ddz, global_min_ddz, 1, mglet_mpi_real, MPI_MIN, MPI_COMM_WORLD)

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
            IF (my_slicecol_list(i)%is_active .OR. myid == 0) THEN

                ALLOCATE(my_slicecol_list(i)%np_counter(mtstep + 1))
                my_slicecol_list(i)%np_counter = 0

                ALLOCATE(my_slicecol_list(i)%rt_counter(rt_tstep_max))
                my_slicecol_list(i)%rt_counter = 0

                IF (tdx_step_max - tdx_step_min /= 0) THEN
                    ALLOCATE(my_slicecol_list(i)%tdx_counter(tdx_step_min:tdx_step_max))
                    my_slicecol_list(i)%tdx_counter = 0
                ELSE
                    ALLOCATE(my_slicecol_list(i)%tdx_counter(0))
                END IF

                IF (tdy_step_max - tdy_step_min /= 0) THEN
                    ALLOCATE(my_slicecol_list(i)%tdy_counter(tdy_step_min:tdy_step_max))
                    my_slicecol_list(i)%tdy_counter = 0
                ELSE
                    ALLOCATE(my_slicecol_list(i)%tdy_counter(0))
                END IF

                IF (tdz_step_max - tdz_step_min /= 0) THEN
                    ALLOCATE(my_slicecol_list(i)%tdz_counter(tdz_step_min:tdz_step_max))
                    my_slicecol_list(i)%tdz_counter = 0
                ELSE
                    ALLOCATE(my_slicecol_list(i)%tdz_counter(0))
                END IF

                IF (avu_step_max - avu_step_min /= 0) THEN
                    ALLOCATE(my_slicecol_list(i)%avu_counter(avu_step_min:avu_step_max))
                    my_slicecol_list(i)%avu_counter = 0
                ELSE
                    ALLOCATE(my_slicecol_list(i)%avu_counter(0))
                END IF

                IF (avv_step_max - avv_step_min /= 0) THEN
                    ALLOCATE(my_slicecol_list(i)%avv_counter(avv_step_min:avv_step_max))
                    my_slicecol_list(i)%avv_counter = 0
                ELSE
                    ALLOCATE(my_slicecol_list(i)%avv_counter(0))
                END IF

                IF (avw_step_max - avw_step_min /= 0) THEN
                    ALLOCATE(my_slicecol_list(i)%avw_counter(avw_step_min:avw_step_max))
                    my_slicecol_list(i)%avw_counter = 0
                ELSE
                    ALLOCATE(my_slicecol_list(i)%avw_counter(0))
                END IF

            ELSE
                ALLOCATE(my_slicecol_list(i)%np_counter(0))
                ALLOCATE(my_slicecol_list(i)%rt_counter(0))
                ALLOCATE(my_slicecol_list(i)%tdx_counter(0))
                ALLOCATE(my_slicecol_list(i)%tdy_counter(0))
                ALLOCATE(my_slicecol_list(i)%tdz_counter(0))
                ALLOCATE(my_slicecol_list(i)%avu_counter(0))
                ALLOCATE(my_slicecol_list(i)%avv_counter(0))
                ALLOCATE(my_slicecol_list(i)%avw_counter(0))
            END IF
        END DO

        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN

            IF (myid /= 0) THEN
                CALL MPI_Recv(dummy, 1, mglet_mpi_int, myid - 1, 900, &
                MPI_COMM_WORLD, MPI_STATUS_IGNORE)
            END IF

            WRITE(*, '("Slice Statistics Collector of Length ", I0, " allocated on Proccess ", I0)') SIZE(my_slicecol_list), myid

            DO i = 1, SIZE(my_slicecol_list)
                WRITE(*, '("Slice ", I0, " activity: ", L2)') i, my_slicecol_list(i)%is_active
            END DO

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
        REAL(realk) :: p_coord

        IF (slice_dir == "N") THEN
            RETURN
        END IF

        CALL start_timer(950)

        IF (particle%islice < 1 .OR. particle%islice > SIZE(my_slicecol_list)) CALL errr(__FILE__, __LINE__)

        SELECT CASE(slice_dir)
            CASE("X")
                    p_coord = particle%x
            CASE("Y")
                    p_coord = particle%y
            CASE("Z")
                    p_coord = particle%z
        END SELECT

        IF (p_coord < my_slicecol_list(particle%islice)%llim) THEN

            IF (ittot >= rt_ittot_start) THEN
                IF (particle%sitstep >= 0) THEN
                    CALL register_slicestat(particle, ittot)
                END IF
                particle%sitstep = ittot
                particle%xyz_sentry(1) = particle%xyz_abs(1)
                particle%xyz_sentry(2) = particle%xyz_abs(2)
                particle%xyz_sentry(3) = particle%xyz_abs(3)
            END IF

            IF (particle%islice == 1) THEN
                particle%islice = SIZE(my_slicecol_list)
            ELSE
                particle%islice = particle%islice - 1
            END IF

        ELSEIF (p_coord > my_slicecol_list(particle%islice)%ulim) THEN

            IF (ittot >= rt_ittot_start) THEN
                IF (particle%sitstep >= 0) THEN
                    CALL register_slicestat(particle, ittot)
                END IF
                particle%sitstep = ittot
                particle%xyz_sentry(1) = particle%xyz_abs(1)
                particle%xyz_sentry(2) = particle%xyz_abs(2)
                particle%xyz_sentry(3) = particle%xyz_abs(3)
            END IF

            IF (particle%islice == SIZE(my_slicecol_list)) THEN
                particle%islice = 1
            ELSE
                particle%islice = particle%islice + 1
            END IF

        END IF

        my_slicecol_list(particle%islice)%np_counter(itstep + 1) = &
         my_slicecol_list(particle%islice)%np_counter(itstep + 1) + 1

        CALL stop_timer(950)

    END SUBROUTINE associate_new_slice

    SUBROUTINE register_slicestat(particle, ittot)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: ittot

        ! local variables
        INTEGER(intk) :: rt_tstep, itdx, itdy, itdz, iavu, iavv, iavw

        rt_tstep = ittot - particle%sitstep
        IF (rt_tstep > SIZE(my_slicecol_list(particle%islice)%rt_counter) .OR. rt_tstep < 1) THEN
            WRITE(*, '("WARNING in Particle Statistics (Slices)")')
            WRITE(*, '("On Proc ", I0, ": Residence time of particle ", I0," is too large or too small to be registered!")') &
             myid, particle%ipart
            WRITE(*, '("rt_tstep ", I0, "; size(rt_counter)", I0)') &
             rt_tstep, SIZE(my_slicecol_list(particle%islice)%rt_counter)
        ELSE
            my_slicecol_list(particle%islice)%rt_counter(rt_tstep) = &
             my_slicecol_list(particle%islice)%rt_counter(rt_tstep) + 1
        END IF

        IF (SIZE(my_slicecol_list(particle%islice)%tdx_counter) > 0) THEN
            itdx = NINT((particle%xyz_abs(1) - particle%xyz_sentry(1)) / global_min_ddx)
            IF (itdx > tdx_step_max .OR. itdx < tdx_step_min) THEN
                !WRITE(*, '("WARNING in Particle Statistics (Slices)")')
                !WRITE(*, '("On Proc ", I0, ": Travel Distance (X) of particle ", I0," is too large or too small to be registered!")') &
                ! myid, particle%ipart
            ELSE
                my_slicecol_list(particle%islice)%tdx_counter(itdx) = &
                 my_slicecol_list(particle%islice)%tdx_counter(itdx) + 1
            END IF
        END IF

        IF (SIZE(my_slicecol_list(particle%islice)%tdy_counter) > 0) THEN
            itdy = NINT((particle%xyz_abs(2) - particle%xyz_sentry(2)) / global_min_ddy)
            IF (itdy > tdy_step_max .OR. itdy < tdy_step_min) THEN
                !WRITE(*, '("WARNING in Particle Statistics (Slices)")')
                !WRITE(*, '("On Proc ", I0, ": Travel Distance (Y) of particle ", I0," is too large or too small to be registered!")') &
                ! myid, particle%ipart
            ELSE
                my_slicecol_list(particle%islice)%tdy_counter(itdy) = &
                 my_slicecol_list(particle%islice)%tdy_counter(itdy) + 1
            END IF
        END IF

        IF (SIZE(my_slicecol_list(particle%islice)%tdz_counter) > 0) THEN
            itdz = NINT((particle%xyz_abs(3) - particle%xyz_sentry(3)) / global_min_ddz)
            IF (itdz > tdz_step_max .OR. itdz < tdz_step_min) THEN
                !WRITE(*, '("WARNING in Particle Statistics (Slices)")')
                !WRITE(*, '("On Proc ", I0, ": Travel Distance (Z) of particle ", I0," is too large or too small to be registered!")') &
                ! myid, particle%ipart
            ELSE
                my_slicecol_list(particle%islice)%tdz_counter(itdz) = &
                 my_slicecol_list(particle%islice)%tdz_counter(itdz) + 1
            END IF
        END IF

        IF (SIZE(my_slicecol_list(particle%islice)%avu_counter) > 0) THEN
            iavu = NINT((particle%xyz_abs(1) - particle%xyz_sentry(1)) / REAL(rt_tstep) &
             / REAL(global_min_ddx) * REAL(MAX(ABS(avu_step_min), ABS(avu_step_max))))
            IF (iavu > avu_step_max .OR. iavu < avu_step_min) THEN
                !WRITE(*, '("WARNING in Particle Statistics (Slices)")')
                !WRITE(*, '("On Proc ", I0, ": Travel Velocity (U) of particle ", I0," is too large or too small to be registered!")') &
                ! myid, particle%ipart
            ELSE
                my_slicecol_list(particle%islice)%avu_counter(iavu) = &
                 my_slicecol_list(particle%islice)%avu_counter(iavu) + 1
            END IF
        END IF

        IF (SIZE(my_slicecol_list(particle%islice)%avv_counter) > 0) THEN
            iavv = NINT((particle%xyz_abs(2) - particle%xyz_sentry(2)) / REAL(rt_tstep) &
             / REAL(global_min_ddy) * REAL(MAX(ABS(avv_step_min), ABS(avv_step_max))))
            IF (iavv > avv_step_max .OR. iavv < avv_step_min) THEN
                !WRITE(*, '("WARNING in Particle Statistics (Slices)")')
                !WRITE(*, '("On Proc ", I0, ": Travel Velocity (V) of particle ", I0," is too large or too small to be registered!")') &
                ! myid, particle%ipart
            ELSE
                my_slicecol_list(particle%islice)%avv_counter(iavv) = &
                 my_slicecol_list(particle%islice)%avv_counter(iavv) + 1
            END IF
        END IF

        IF (SIZE(my_slicecol_list(particle%islice)%avw_counter) > 0) THEN
            iavw = NINT((particle%xyz_abs(3) - particle%xyz_sentry(3)) / REAL(rt_tstep) &
             / REAL(global_min_ddz) * REAL(MAX(ABS(avw_step_min), ABS(avw_step_max))))
            IF (iavw > avw_step_max .OR. iavw < avw_step_min) THEN
                !WRITE(*, '("WARNING in Particle Statistics (Slices)")')
                !WRITE(*, '("On Proc ", I0, ": Travel Velocity (U) of particle ", I0," is too large or too small to be registered!")') &
                ! myid, particle%ipart
            ELSE
                my_slicecol_list(particle%islice)%avw_counter(iavw) = &
                 my_slicecol_list(particle%islice)%avw_counter(iavw) + 1
            END IF
        END IF

    END SUBROUTINE register_slicestat


    SUBROUTINE merge_slicestat(islice, mtstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: islice
        INTEGER(intk), INTENT(in) :: mtstep

        ! local variables
        INTEGER(intk) :: i, j, counter, buffer_len, dsend, displ
        INTEGER(intk), ALLOCATABLE :: drecv(:)
        INTEGER(intk), ALLOCATABLE :: buffer(:)
        TYPE(MPI_Request) :: send_req
        TYPE(MPI_Request), ALLOCATABLE :: recv_req(:)

        !TODO: optimize this dirty routine

        IF (slice_dir == "N") THEN
            RETURN
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

        IF (myid /= 0) THEN

            IF (dsend == 1) THEN

                buffer_len = mtstep + 1 + rt_tstep_max &
                + SIZE(my_slicecol_list(islice)%tdx_counter) + SIZE(my_slicecol_list(islice)%tdy_counter) &
                + SIZE(my_slicecol_list(islice)%tdz_counter) + SIZE(my_slicecol_list(islice)%avu_counter) &
                + SIZE(my_slicecol_list(islice)%avv_counter) + SIZE(my_slicecol_list(islice)%avw_counter)

                ALLOCATE(buffer(buffer_len))
                buffer = 0

                ! FILL BUFFER
                ! npart
                displ = 1
                counter = 1
                DO j = displ, displ + mtstep
                    buffer(j) = my_slicecol_list(islice)%np_counter(counter)
                    counter = counter + 1
                END DO

                ! residence time
                displ = j
                counter = 1
                DO j = displ, displ + rt_tstep_max - 1
                    buffer(j) = my_slicecol_list(islice)%rt_counter(counter)
                    counter = counter + 1
                END DO

                ! travel distance x
                IF (SIZE(my_slicecol_list(islice)%tdx_counter) > 0) THEN
                    displ = j
                    counter = tdx_step_min
                    DO j = displ, displ + SIZE(my_slicecol_list(islice)%tdx_counter) - 1
                        buffer(j) = my_slicecol_list(islice)%tdx_counter(counter)
                        counter = counter + 1
                    END DO
                END IF

                ! travel distance y
                IF (SIZE(my_slicecol_list(islice)%tdy_counter) > 0) THEN
                    displ = j
                    counter = tdy_step_min
                    DO j = displ, displ + SIZE(my_slicecol_list(islice)%tdy_counter) - 1
                        buffer(j) = my_slicecol_list(islice)%tdy_counter(counter)
                        counter = counter + 1
                    END DO
                END IF

                ! travel distance z
                IF (SIZE(my_slicecol_list(islice)%tdz_counter) > 0) THEN
                    displ = j
                    counter = tdz_step_min
                    DO j = displ, displ + SIZE(my_slicecol_list(islice)%tdz_counter) - 1
                        buffer(j) = my_slicecol_list(islice)%tdz_counter(counter)
                        counter = counter + 1
                    END DO
                END IF

                ! average velocity u
                IF (SIZE(my_slicecol_list(islice)%avu_counter) > 0) THEN
                    displ = j
                    counter = avu_step_min
                    DO j = displ, displ + SIZE(my_slicecol_list(islice)%avu_counter) - 1
                        buffer(j) = my_slicecol_list(islice)%avu_counter(counter)
                        counter = counter + 1
                    END DO
                END IF

                ! average velocity v
                IF (SIZE(my_slicecol_list(islice)%avv_counter) > 0) THEN
                    displ = j
                    counter = avv_step_min
                    DO j = displ, displ + SIZE(my_slicecol_list(islice)%avv_counter) - 1
                        buffer(j) = my_slicecol_list(islice)%avv_counter(counter)
                        counter = counter + 1
                    END DO
                END IF

                ! average velocity w
                IF (SIZE(my_slicecol_list(islice)%avw_counter) > 0) THEN
                    displ = j
                    counter = avw_step_min
                    DO j = displ, displ + SIZE(my_slicecol_list(islice)%avw_counter) - 1
                        buffer(j) = my_slicecol_list(islice)%avw_counter(counter)
                        counter = counter + 1
                    END DO
                END IF

                CALL MPI_Send(buffer, buffer_len, mglet_mpi_int, 0, 9402, &
                 MPI_COMM_WORLD)

                !CALL MPI_Wait(send_req, MPI_STATUS_IGNORE)

            END IF

        ELSEIF (myid == 0) THEN

            buffer_len = mtstep + 1 + rt_tstep_max &
            + SIZE(my_slicecol_list(islice)%tdx_counter) + SIZE(my_slicecol_list(islice)%tdy_counter) &
            + SIZE(my_slicecol_list(islice)%tdz_counter) + SIZE(my_slicecol_list(islice)%avu_counter) &
            + SIZE(my_slicecol_list(islice)%avv_counter) + SIZE(my_slicecol_list(islice)%avw_counter)

            ALLOCATE(buffer(buffer_len))
            buffer = 0

            DO i = 1, numprocs - 1

                IF (drecv(i) == 1) THEN

                    buffer = 0

                    CALL MPI_Recv(buffer, buffer_len, mglet_mpi_int, i, 9402, &
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE)

                    displ = 1
                    counter = 1
                    DO j = displ, displ + mtstep
                        my_slicecol_list(islice)%np_counter(counter) = my_slicecol_list(islice)%np_counter(counter) + buffer(j)
                        counter = counter + 1
                    END DO

                    displ = j
                    counter = 1
                    DO j = displ, displ + rt_tstep_max - 1
                        my_slicecol_list(islice)%rt_counter(counter) = my_slicecol_list(islice)%rt_counter(counter) + buffer(j)
                        counter = counter + 1
                    END DO

                    IF (SIZE(my_slicecol_list(islice)%tdx_counter) > 0) THEN
                        displ = j
                        counter = tdx_step_min
                        DO j = displ, displ + SIZE(my_slicecol_list(islice)%tdx_counter) - 1
                            my_slicecol_list(islice)%tdx_counter(counter) = my_slicecol_list(islice)%tdx_counter(counter) + buffer(j)
                            counter = counter + 1
                        END DO
                    END IF

                    IF (SIZE(my_slicecol_list(islice)%tdy_counter) > 0) THEN
                        displ = j
                        counter = tdy_step_min
                        DO j = displ, displ + SIZE(my_slicecol_list(islice)%tdy_counter) - 1
                            my_slicecol_list(islice)%tdy_counter(counter) = my_slicecol_list(islice)%tdy_counter(counter) + buffer(j)
                            counter = counter + 1
                        END DO
                    END IF

                    IF (SIZE(my_slicecol_list(islice)%tdz_counter) > 0) THEN
                        displ = j
                        counter = tdz_step_min
                        DO j = displ, displ + SIZE(my_slicecol_list(islice)%tdz_counter) - 1
                            my_slicecol_list(islice)%tdz_counter(counter) = my_slicecol_list(islice)%tdz_counter(counter) + buffer(j)
                            counter = counter + 1
                        END DO
                    END IF

                    IF (SIZE(my_slicecol_list(islice)%avu_counter) > 0) THEN
                        displ = j
                        counter = avu_step_min
                        DO j = displ, displ + SIZE(my_slicecol_list(islice)%avu_counter) - 1
                            my_slicecol_list(islice)%avu_counter(counter) = my_slicecol_list(islice)%avu_counter(counter) + buffer(j)
                            counter = counter + 1
                        END DO
                    END IF

                    IF (SIZE(my_slicecol_list(islice)%avv_counter) > 0) THEN
                        displ = j
                        counter = avv_step_min
                        DO j = displ, displ + SIZE(my_slicecol_list(islice)%avv_counter) - 1
                            my_slicecol_list(islice)%avv_counter(counter) = my_slicecol_list(islice)%avv_counter(counter) + buffer(j)
                            counter = counter + 1
                        END DO
                    END IF

                    IF (SIZE(my_slicecol_list(islice)%avw_counter) > 0) THEN
                        displ = j
                        counter = avw_step_min
                        DO j = displ, displ + SIZE(my_slicecol_list(islice)%avw_counter) - 1
                            my_slicecol_list(islice)%avw_counter(counter) = my_slicecol_list(islice)%avw_counter(counter) + buffer(j)
                            counter = counter + 1
                        END DO
                    END IF

                ELSE

                    CYCLE

                END IF

            END DO

        END IF

        IF (ALLOCATED(buffer)) DEALLOCATE(buffer)
        IF (ALLOCATED(recv_req)) DEALLOCATE(recv_req)
        IF (ALLOCATED(drecv)) DEALLOCATE(drecv)

    END SUBROUTINE merge_slicestat

    SUBROUTINE deallocate_slicestat(islice)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: islice

        IF (.NOT. ALLOCATED(my_slicecol_list)) RETURN

        IF (islice < 1 .OR. islice > SIZE(my_slicecol_list)) RETURN

        IF (ALLOCATED(my_slicecol_list(islice)%np_counter)) DEALLOCATE(my_slicecol_list(islice)%np_counter)
        IF (ALLOCATED(my_slicecol_list(islice)%rt_counter)) DEALLOCATE(my_slicecol_list(islice)%rt_counter)
        IF (ALLOCATED(my_slicecol_list(islice)%rt_counter)) DEALLOCATE(my_slicecol_list(islice)%tdx_counter)
        IF (ALLOCATED(my_slicecol_list(islice)%rt_counter)) DEALLOCATE(my_slicecol_list(islice)%tdy_counter)
        IF (ALLOCATED(my_slicecol_list(islice)%rt_counter)) DEALLOCATE(my_slicecol_list(islice)%tdz_counter)
        IF (ALLOCATED(my_slicecol_list(islice)%rt_counter)) DEALLOCATE(my_slicecol_list(islice)%avu_counter)
        IF (ALLOCATED(my_slicecol_list(islice)%rt_counter)) DEALLOCATE(my_slicecol_list(islice)%avv_counter)
        IF (ALLOCATED(my_slicecol_list(islice)%rt_counter)) DEALLOCATE(my_slicecol_list(islice)%avw_counter)

    END SUBROUTINE deallocate_slicestat

    SUBROUTINE deallocate_gridstat(i)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: i

        IF (.NOT. ALLOCATED(my_gridcol_list)) RETURN

        IF (i < 1 .OR. i > SIZE(my_gridcol_list)) RETURN

        IF (ALLOCATED(my_gridcol_list(i)%np_counter)) DEALLOCATE(my_gridcol_list(i)%np_counter)
        IF (ALLOCATED(my_gridcol_list(i)%rt_counter)) DEALLOCATE(my_gridcol_list(i)%rt_counter)

    END SUBROUTINE deallocate_gridstat

    SUBROUTINE write_particle_statistics(mtstep, dt)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: mtstep
        REAL(realk), INTENT(in) :: dt

        ! local variables
        INTEGER(intk) :: islice

        CALL start_timer(900)
        CALL start_timer(950)

        CALL write_partstat_folder()

        CALL write_gridstat_files()

        DO islice = 1, SIZE(my_slicecol_list)
            CALL merge_slicestat(islice, mtstep)
            CALL write_slicestat_file(islice, dt)
        END DO

        CALL stop_timer(950)
        CALL stop_timer(900)

    END SUBROUTINE write_particle_statistics

    SUBROUTINE write_partstat_folder()

        IF (myid == 0) THEN
            CALL create_directory("Particle_Statistics") ! ! ! realtive to working directory ! ! !
        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

    END SUBROUTINE write_partstat_folder

    SUBROUTINE write_gridstat_files()

        INTEGER(intk) :: i, j, unit
        CHARACTER(len = mglet_filename_max) :: filename

        IF (.NOT. dgridstat) THEN
            RETURN
        END IF

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

    SUBROUTINE write_slicestat_file(islice, dt)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: islice
        REAL(realk), INTENT(in) :: dt

        ! local varibales
        INTEGER(intk) :: j, counter, unit
        REAL(realk) :: avu_increment, avv_increment, avw_increment
        CHARACTER(len = mglet_filename_max) :: filename

        IF (myid == 0) THEN

            WRITE(filename,'("Particle_Statistics/Slice", I0, ".txt")') islice

            OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

            WRITE(unit, '("ISLICE:         ", I0)') islice

            WRITE(unit, '("SLICE LIMITS:   ", 2F10.7)') my_slicecol_list(islice)%llim, my_slicecol_list(islice)%ulim

            WRITE(unit, '("TIME INCREMENT: ", 1F10.7)') dt

            IF (SIZE(my_slicecol_list(islice)%tdx_counter) > 0) THEN
                WRITE(unit, '("TDX STENCIL:    ", 1F10.7, A, 1F10.7, A, 1F10.7)') &
                 global_min_ddx*REAL(tdx_step_min), " : ", global_min_ddx, " : ", global_min_ddx*REAL(tdx_step_max)
            END IF
            IF (SIZE(my_slicecol_list(islice)%tdy_counter) > 0) THEN
                WRITE(unit, '("TDY STENCIL:    ", 1F10.7, A, 1F10.7, A, 1F10.7)') &
                 global_min_ddy*REAL(tdy_step_min), " : ", global_min_ddy, " : ", global_min_ddy*REAL(tdy_step_max)
            END IF
            IF (SIZE(my_slicecol_list(islice)%tdz_counter) > 0) THEN
                WRITE(unit, '("TDZ STENCIL:    ", 1F10.7, A, 1F10.7, A, 1F10.7)') &
                 global_min_ddz*REAL(tdz_step_min), " : ", global_min_ddz, " : ", global_min_ddz*REAL(tdz_step_max)
            END IF

            IF (SIZE(my_slicecol_list(islice)%avu_counter) > 0) THEN
                avu_increment = global_min_ddx/dt/REAL(MAX(ABS(avu_step_min), ABS(avu_step_max)))
                WRITE(unit, '("AVU STENCIL:    ", 1F10.7, A, 1F10.7, A, 1F10.7)') &
                 avu_increment*REAL(avu_step_min), " : ", avu_increment, " : ", avu_increment*REAL(avu_step_max)
            END IF
            IF (SIZE(my_slicecol_list(islice)%avv_counter) > 0) THEN
                avv_increment = global_min_ddy/dt/MAX(ABS(avv_step_min), ABS(avv_step_max))
                WRITE(unit, '("AVV STENCIL:    ", 1F10.7, A, 1F10.7, A, 1F10.7)') &
                 avv_increment*REAL(avv_step_min), " : ", avv_increment, " : ", avv_increment*REAL(avv_step_max)
            END IF
            IF (SIZE(my_slicecol_list(islice)%avw_counter) > 0) THEN
                avw_increment = global_min_ddz/dt/MAX(ABS(avw_step_min), ABS(avw_step_max))
                WRITE(unit, '("AVW STENCIL:    ", 1F10.7, A, 1F10.7, A, 1F10.7)') &
                avw_increment*REAL(avw_step_min), " : ", avw_increment, " : ", avw_increment*REAL(avw_step_max)
            END IF

            WRITE(unit, '("")', advance = "yes")

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

            WRITE(unit, '("")', advance = "yes")

            IF (SIZE(my_slicecol_list(islice)%tdx_counter) > 0) THEN
                WRITE(unit, '("TDX DISTRIBUTION:")')
                counter = 1
                DO j = tdx_step_min, tdx_step_max
                    WRITE(unit, '(I0)', advance = "no") my_slicecol_list(islice)%tdx_counter(j)
                    IF (j /= tdx_step_max) THEN
                        WRITE(unit, '(", ")', advance = "no")
                        IF (MOD(counter, 10) == 0) THEN
                            WRITE(unit, '("")', advance = "yes")
                        END IF
                    END IF
                    counter = counter + 1
                END DO
                WRITE(unit, '("")', advance = "yes")
            END IF

            IF (SIZE(my_slicecol_list(islice)%tdy_counter) > 0) THEN
                WRITE(unit, '("TDY DISTRIBUTION:")')
                counter = 1
                DO j = tdy_step_min, tdy_step_max
                    WRITE(unit, '(I0)', advance = "no") my_slicecol_list(islice)%tdy_counter(j)
                    IF (j /= tdy_step_max) THEN
                        WRITE(unit, '(", ")', advance = "no")
                        IF (MOD(counter, 10) == 0) THEN
                            WRITE(unit, '("")', advance = "yes")
                        END IF
                    END IF
                    counter = counter + 1
                END DO
                WRITE(unit, '("")', advance = "yes")
            END IF

            IF (SIZE(my_slicecol_list(islice)%tdz_counter) > 0) THEN
                WRITE(unit, '("TDZ DISTRIBUTION:")')
                counter = 1
                DO j = tdz_step_min, tdz_step_max
                    WRITE(unit, '(I0)', advance = "no") my_slicecol_list(islice)%tdz_counter(j)
                    IF (j /= tdz_step_max) THEN
                        WRITE(unit, '(", ")', advance = "no")
                        IF (MOD(counter, 10) == 0) THEN
                            WRITE(unit, '("")', advance = "yes")
                        END IF
                    END IF
                    counter = counter + 1
                END DO
                WRITE(unit, '("")', advance = "yes")
            END IF

            IF (SIZE(my_slicecol_list(islice)%avu_counter) > 0) THEN
                WRITE(unit, '("AVU DISTRIBUTION:")')
                counter = 1
                DO j = avu_step_min, avu_step_max
                    WRITE(unit, '(I0)', advance = "no") my_slicecol_list(islice)%avu_counter(j)
                    IF (j /= avu_step_max) THEN
                        WRITE(unit, '(", ")', advance = "no")
                        IF (MOD(counter, 10) == 0) THEN
                            WRITE(unit, '("")', advance = "yes")
                        END IF
                    END IF
                    counter = counter + 1
                END DO
                WRITE(unit, '("")', advance = "yes")
            END IF

            IF (SIZE(my_slicecol_list(islice)%avv_counter) > 0) THEN
                WRITE(unit, '("AVV DISTRIBUTION:")')
                counter = 1
                DO j = avv_step_min, avv_step_max
                    WRITE(unit, '(I0)', advance = "no") my_slicecol_list(islice)%avv_counter(j)
                    IF (j /= avv_step_max) THEN
                        WRITE(unit, '(", ")', advance = "no")
                        IF (MOD(counter, 10) == 0) THEN
                            WRITE(unit, '("")', advance = "yes")
                        END IF
                    END IF
                    counter = counter + 1
                END DO
                WRITE(unit, '("")', advance = "yes")
            END IF

            IF (SIZE(my_slicecol_list(islice)%avw_counter) > 0) THEN
                WRITE(unit, '("AVW DISTRIBUTION:")')
                counter = 1
                DO j = avw_step_min, avw_step_max
                    WRITE(unit, '(I0)', advance = "no") my_slicecol_list(islice)%avw_counter(j)
                    IF (j /= avw_step_max) THEN
                        WRITE(unit, '(", ")', advance = "no")
                        IF (MOD(counter, 10) == 0) THEN
                            WRITE(unit, '("")', advance = "yes")
                        END IF
                    END IF
                    counter = counter + 1
                END DO
            END IF

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