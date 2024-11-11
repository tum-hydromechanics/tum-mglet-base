MODULE particle_dict_mod

    ! This module is responsible for:
    ! Reading of initial particle coordinates from the ParticleDict.txt file.

    USE grids_mod

    USE particle_config_mod

    IMPLICIT NONE

    PUBLIC :: read_particles

CONTAINS    !===================================

    SUBROUTINE read_particles(dread_particles, dict_len, ipart_arr, igrid_arr, x_arr, y_arr, z_arr, read_np)

        !subroutine arguments
        LOGICAL, INTENT(inout) :: dread_particles
        INTEGER(intk), INTENT(out) :: dict_len, read_np
        INTEGER(intk), ALLOCATABLE, INTENT(inout) :: ipart_arr(:), igrid_arr(:)
        REAL(realk), ALLOCATABLE, INTENT(inout) :: x_arr(:), y_arr(:), z_arr(:)

        !local variables
        INTEGER(intk) :: i, ipart, igrid, unit
        REAL(realk) :: xtemp, ytemp, ztemp
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        INQUIRE(file = 'ParticleDict.txt', exist = dread_particles)

        IF (.NOT. dread_particles) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: No file for reading particles detected! Using automated initial particle distribution instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) ' '
                        WRITE(*, *) "WARNING: No file for reading particles detected! Using automated initial particle distribution instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            RETURN

        END IF

        ! CAUTION: the following is not optimized for multiple processes !

        OPEN(newunit = unit, file = 'ParticleDict.txt', status = 'OLD', action = 'READ')

        READ(unit, fmt = *) dict_len

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, '("READING ", I0, " PARTICLE(S) ON ", I0, " PROCESSES.")') dict_len, numprocs
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, '("READING ", I0, " PARTICLE(S) ON ", I0, " PROCESSES.")') dict_len, numprocs
                    WRITE(*, '()')
            END SELECT
        END IF

        IF (list_limit) THEN
            ALLOCATE(ipart_arr(plist_len))
            ALLOCATE(igrid_arr(plist_len))
            ALLOCATE(x_arr(plist_len))
            ALLOCATE(y_arr(plist_len))
            ALLOCATE(z_arr(plist_len))
        ELSE
            ALLOCATE(ipart_arr(dict_len))
            ALLOCATE(igrid_arr(dict_len))
            ALLOCATE(x_arr(dict_len))
            ALLOCATE(y_arr(dict_len))
            ALLOCATE(z_arr(dict_len))
        END IF

        ! ParticleDict.txt is screened from top to bottom.
        ! If a particle is found to lie on a grid of this process, the particle is stored on this process.
        ! Once the particle list length or dict length has been reached, no more particles are read on this process.
        ! Hence, depending on the parameterization of the particle list and the dict length, some particles might not be registered!

        read_np = 0

        DO ipart = 1, dict_len

            READ(unit, fmt = *) xtemp, ytemp, ztemp

            DO i = 1, nmygrids

                igrid = mygrids(i)
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

                read_np = read_np + 1

                ipart_arr(read_np) = ipart
                igrid_arr(read_np) = igrid
                x_arr(read_np) = xtemp
                y_arr(read_np) = ytemp
                z_arr(read_np) = ztemp

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*,'("Particle read on proc ", I0, ": ID = ", I0, " | x/y/z = ", 3F12.6)') myid, ipart, xtemp, ytemp, ztemp
                END SELECT

                EXIT

            END DO

            IF (SIZE(ipart_arr) == read_np) THEN
                IF (ipart < dict_len) THEN
                    SELECT CASE (TRIM(particle_terminal))
                        CASE ("none")
                            CONTINUE
                        CASE ("normal")
                            CONTINUE
                        CASE ("verbose")
                            WRITE(*,'("Warning on proc ", I0, ": Maximum Number of Particles has been registered on this Proccess.", &
                             "Stopped reading ParticleDict.txt, so specified Particles might be unregistered." )') myid
                    END SELECT
                END IF
                EXIT
            END IF

        END DO

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, '()')
                    WRITE(*, '("READING OF PARTICLES SUCCESSFULLY COMPLETED.")')
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, '()')
                    WRITE(*, '("READING OF PARTICLES SUCCESSFULLY COMPLETED.")')
                    WRITE(*, '()')
            END SELECT
        END IF

        CLOSE(unit)

    END SUBROUTINE read_particles

END MODULE particle_dict_mod
