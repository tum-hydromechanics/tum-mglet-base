MODULE particle_gridstat_mod

    USE utils_mod

    USE particle_list_mod

    IMPLICIT NONE

    TYPE :: gridstat_collector_t

        INTEGER(intk) :: igrid

        ! volume of the grid
        REAL(realk) :: grid_volume

        ! array that stores the number of particles on the grid for each timestep
        INTEGER(intk), ALLOCATABLE :: np_counter(:)

        ! array that stores the residence time distribution
        INTEGER(intk), ALLOCATABLE :: rt_counter(:)

    END TYPE gridstat_collector_t

    TYPE(gridstat_collector_t), ALLOCATABLE :: my_collector_list(:)

CONTAINS

    SUBROUTINE init_particle_gridstat(mtstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: mtstep

        ! local variables
        INTEGER(intk) :: igrid, i, j
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        LOGICAL :: gridstat_exists

        CALL start_timer(900)
        CALL start_timer(940)

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

        ! allocate collector list
        ALLOCATE(my_collector_list(nmygrids))

        ! allocate arrays of the individual collectors
        DO i = 1, nmygrids

            igrid = mygrids(i)

            my_collector_list(i)%igrid = igrid

            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            my_collector_list(i)%grid_volume = (maxx - minx) * (maxy -miny) * (maxz - minz)

            ALLOCATE(my_collector_list(i)%np_counter(mtstep + 1))
            my_collector_list(i)%np_counter = 0

            ALLOCATE(my_collector_list(i)%rt_counter(mtstep))
            my_collector_list(i)%rt_counter = 0

            ! initialize first entry of np_counter (count the number of particles in each grid at t = 0 resprectively)
            DO j = 1, my_particle_list%ifinal

                IF (my_particle_list%particles(j)%igrid == igrid) THEN

                    my_collector_list(i)%np_counter(1) = my_collector_list(i)%np_counter(1) + 1

                END IF

            END DO

        END DO

        CALL stop_timer(940)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_gridstat

    SUBROUTINE advance_np_counter(itstep)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep

        ! local variables
        INTEGER(intk) :: i

        CALL start_timer(900)
        CALL start_timer(940)

        DO i = 1, nmygrids

            my_collector_list(i)%np_counter(itstep + 1) = my_collector_list(i)%np_counter(itstep)

        END DO

        CALL stop_timer(940)
        CALL stop_timer(900)

    END SUBROUTINE advance_np_counter

    SUBROUTINE deregister_particle(particle, itstep)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: itstep

        ! local variables
        INTEGER(intk) :: igrid, i, irt

        CALL start_timer(940)

        DO i = 1, nmygrids

            igrid = mygrids(i)

            IF (particle%igrid == igrid) THEN

                my_collector_list(i)%np_counter(itstep + 1) = my_collector_list(i)%np_counter(itstep + 1) - 1

                irt = itstep - particle%itstep

                my_collector_list(i)%rt_counter(irt) = my_collector_list(i)%rt_counter(irt) + 1

            END IF

        END DO

        particle%itstep = itstep

        CALL stop_timer(940)

    END SUBROUTINE deregister_particle

    SUBROUTINE register_particle(particle, itstep)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: itstep

        ! local variables
        INTEGER(intk) :: igrid, i

        CALL start_timer(940)

        DO i = 1, nmygrids

            igrid = mygrids(i)

            IF (particle%igrid == igrid) THEN

                my_collector_list(i)%np_counter(itstep + 1) = my_collector_list(i)%np_counter(itstep + 1) + 1

            END IF

        END DO

        CALL stop_timer(940)

    END SUBROUTINE register_particle

    SUBROUTINE write_gridstat()

        CALL start_timer(900)
        CALL start_timer(940)

        CALL write_gridstat_folder()

        CALL write_gridstat_files()

        CALL stop_timer(940)
        CALL stop_timer(900)

    END SUBROUTINE write_gridstat

    SUBROUTINE write_gridstat_folder()

        IF (myid == 0) THEN
            CALL create_directory("Particle_Statistics") ! ! ! realtive to working directory ! ! !
        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

    END SUBROUTINE write_gridstat_folder

    SUBROUTINE write_gridstat_files()

        INTEGER(intk) :: i, j, unit = 164
        CHARACTER(len = mglet_filename_max) :: filename

        DO i = 1, nmygrids

            WRITE(filename,'("Particle_Statistics/Grid", I0, ".txt")') my_collector_list(i)%igrid

            OPEN(unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

            WRITE(unit, '("IGRID: ", I0)') my_collector_list(i)%igrid

            WRITE(unit, '("GRID VOLUME: ", F12.6)') my_collector_list(i)%grid_volume

            WRITE(unit, '("NPART TIMELINE:")')

            DO j = 1, SIZE(my_collector_list(i)%np_counter)

                WRITE(unit, '(I0)', advance = "no") my_collector_list(i)%np_counter(j)

                IF (j /= SIZE(my_collector_list(i)%np_counter)) THEN
                    WRITE(unit, '(", ")', advance = "no")
                END IF

            END DO

            WRITE(unit, '("")', advance = "yes")

            DO j = 1, SIZE(my_collector_list(i)%rt_counter)

                WRITE(unit, '(I0)', advance = "no") my_collector_list(i)%rt_counter(j)

                IF (j /= SIZE(my_collector_list(i)%rt_counter)) THEN
                    WRITE(unit, '(", ")', advance = "no")
                END IF

            END DO

        END DO

    END SUBROUTINE write_gridstat_files

    SUBROUTINE finish_particle_gridstat()

        ! local variables

        INTEGER(intk) :: i

        DO i = 1, nmygrids
            DEALLOCATE(my_collector_list(i)%np_counter)
            DEALLOCATE(my_collector_list(i)%rt_counter)
        END DO

        DEALLOCATE(my_collector_list)

    END SUBROUTINE finish_particle_gridstat

END MODULE particle_gridstat_mod