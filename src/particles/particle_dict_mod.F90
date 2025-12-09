MODULE particle_dict_mod

    ! This module is responsible for:
    ! Reading of initial particle coordinates from a ParticleDict.txt file.

    USE grids_mod

    USE particle_config_mod
    USE particle_basetype_mod

    IMPLICIT NONE

    PUBLIC :: read_particles

CONTAINS

    SUBROUTINE read_particles(dread_particles_dict, ipart_arr, igrid_arr, x_arr, y_arr, z_arr, read_np)

        !subroutine arguments
        LOGICAL, INTENT(inout) :: dread_particles_dict
        INTEGER(intk), INTENT(out) :: read_np
        INTEGER(intk), ALLOCATABLE, INTENT(inout) :: ipart_arr(:), igrid_arr(:)
        REAL(realk), ALLOCATABLE, INTENT(inout) :: x_arr(:), y_arr(:), z_arr(:)

        !local variables
        LOGICAL :: grid_found
        INTEGER(intk) :: i, j, ipart, igrid, unit, dict_np, ntemp, global_np
        REAL(realk) :: xtemp, ytemp, ztemp
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        INQUIRE(file = 'ParticleDict.txt', exist = dread_particles_dict)

        IF (.NOT. dread_particles_dict) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: No ParticleDict.txt detected! Using automated initial particle distribution instead."
                    WRITE(*, '()')
                END IF
            END IF

            RETURN

        END IF

        ! CAUTION: the following is not optimized for multiple processes !

        OPEN(newunit = unit, file = 'ParticleDict.txt', status = 'OLD', action = 'READ')

        READ(unit, fmt = *) dict_np

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("READING ", I0, " PARTICLE(S) ON ", I0, " PROCESSES.")') dict_np, numprocs
                WRITE(*, '()')
            END IF
        END IF

        IF (list_limit) THEN
            ALLOCATE(ipart_arr(plist_len))
            ALLOCATE(igrid_arr(plist_len))
            ALLOCATE(x_arr(plist_len))
            ALLOCATE(y_arr(plist_len))
            ALLOCATE(z_arr(plist_len))
        ELSE
            ALLOCATE(ipart_arr(dict_np))
            ALLOCATE(igrid_arr(dict_np))
            ALLOCATE(x_arr(dict_np))
            ALLOCATE(y_arr(dict_np))
            ALLOCATE(z_arr(dict_np))
        END IF

        ! ParticleDict.txt is screened from top to bottom.
        ! If a particle is found to lie on a grid of this process, the particle is stored on this process.
        ! Once the particle list length or dict length has been reached, no more particles are read on this process.
        ! Hence, depending on the parameterization of the particle list and the dict length, some particles might not be registered!

        read_np = 0
        ipart = 0

        DO WHILE (ipart < dict_np .AND. read_np < SIZE(ipart_arr))

            READ(unit, fmt = *) ntemp, xtemp, ytemp, ztemp

            grid_found = .FALSE.

            DO i = 1, nmygridslvl(particle_level)

                igrid = mygridslvl(i, particle_level)
                CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

                IF (xtemp < minx) THEN
                    CYCLE
                END IF

                IF (xtemp >= maxx) THEN
                    CYCLE
                END IF

                IF (ytemp < miny) THEN
                    CYCLE
                END IF

                IF (ytemp >= maxy) THEN
                    CYCLE
                END IF

                IF (ztemp < minz) THEN
                    CYCLE
                END IF

                IF (ztemp >= maxz) THEN
                    CYCLE
                END IF

                grid_found = .TRUE.

                DO j = 1, ntemp

                    read_np = read_np + 1
                    ipart = ipart + 1
                    ipart_arr(read_np) = ipart
                    igrid_arr(read_np) = igrid
                    x_arr(read_np) = xtemp
                    y_arr(read_np) = ytemp
                    z_arr(read_np) = ztemp

                    IF (TRIM(particle_terminal) == "verbose") THEN
                        WRITE(*,'("Particle read on proc ", I0, ": ID = ", I0, " | x/y/z = ", 3F12.6)') myid, ipart, xtemp, ytemp, ztemp
                        WRITE(*, '()')
                    END IF

                    IF (SIZE(ipart_arr) == read_np .OR. ipart == dict_np) THEN
                        EXIT
                    END IF

                END DO

                EXIT

            END DO

            IF (.NOT. grid_found) THEN
                ipart = ipart + ntemp
            END IF

        END DO

        CALL MPI_Allreduce(read_np, global_np, 1, mglet_mpi_int, MPI_SUM, MPI_COMM_WORLD)

        IF (myid == 0) THEN
            IF (global_np < dict_np) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("Warning: The number of registered particles is smaller than the specified number of particles in ParticleDict.txt.")')
                    WRITE(*,'("Warning: This is likely caused by a limited list length or incosistencies/invalid positions in the ParticleDict.txt file.")')
                    WRITE(*, '()')
                END IF
            END IF
        END IF

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("READING OF PARTICLES SUCCESSFULLY COMPLETED.")')
                WRITE(*, '()')
            END IF
        END IF

        CLOSE(unit)

        ! Deallocation is performed in calling routine

    END SUBROUTINE read_particles

END MODULE particle_dict_mod
