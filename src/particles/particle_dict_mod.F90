MODULE particle_dict_mod

    ! This module is responsible for:
    ! Reading of initial particle coordinates from a ParticleDict.txt file.

    USE grids_mod

    USE particle_config_mod
    USE particle_core_mod

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
        INTEGER(intk) :: i, ipart, igrid, unit, dict_len
        REAL(realk) :: xtemp, ytemp, ztemp
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        INQUIRE(file = 'ParticleDict.txt', exist = dread_particles_dict)

        IF (.NOT. dread_particles_dict) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: No file for reading particles detected! Using automated initial particle distribution instead."
                    WRITE(*, '()')
                END IF
            END IF

            RETURN

        END IF

        ! CAUTION: the following is not optimized for multiple processes !

        OPEN(newunit = unit, file = 'ParticleDict.txt', status = 'OLD', action = 'READ')

        READ(unit, fmt = *) dict_len

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("READING ", I0, " PARTICLE(S) ON ", I0, " PROCESSES.")') dict_len, numprocs
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

                read_np = read_np + 1

                ipart_arr(read_np) = ipart
                igrid_arr(read_np) = igrid
                x_arr(read_np) = xtemp
                y_arr(read_np) = ytemp
                z_arr(read_np) = ztemp

                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("Particle read on proc ", I0, ": ID = ", I0, " | x/y/z = ", 3F12.6)') myid, ipart, xtemp, ytemp, ztemp
                    WRITE(*, '()')
                END IF

                EXIT

            END DO

            IF (SIZE(ipart_arr) == read_np) THEN
                IF (ipart < dict_len) THEN
                    IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                        WRITE(*,'("Warning on proc ", I0, ": Maximum Number of Particles has been registered on this Proccess.")') myid
                        WRITE(*, '("Stopped reading ParticleDict.txt, so specified Particles might be unregistered.")')
                        WRITE(*, '()')
                    END IF
                END IF
                EXIT
            END IF

        END DO

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
