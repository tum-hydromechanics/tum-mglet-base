MODULE particle_dict_mod

    USE core_mod
    USE particle_config_mod
    IMPLICIT NONE

CONTAINS

    SUBROUTINE read_particles(dread_particles, max_npart, npart, ipart_arr, p_igrid_arr, x, y, z)

    !subroutine arguments
    LOGICAL, INTENT(inout) :: dread_particles
    INTEGER(intk), INTENT(in) :: max_npart
    INTEGER(intk), INTENT(out) :: npart
    INTEGER(intk), ALLOCATABLE, INTENT(out) :: ipart_arr(:)
    INTEGER(intk), ALLOCATABLE, INTENT(out) :: p_igrid_arr(:)
    REAL(realk), ALLOCATABLE, INTENT(out) :: x(:), y(:), z(:)

    !local variables
    INTEGER(intk) :: i, ipart, igrid
    INTEGER :: unit
    REAL(realk) :: xtemp, ytemp, ztemp
    REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

    INQUIRE(file = 'ParticleDict.txt', exist = dread_particles)

    IF (.NOT. dread_particles) THEN

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                WRITE(*, *) ' '
                WRITE(*, *) "WARNING: No file for reading particles detected! Using automated initial particle distribution instead."
            CASE ("verbose")
                WRITE(*, *) ' '
                WRITE(*, *) "WARNING: No file for reading particles detected! Using automated initial particle distribution instead."
        END SELECT

        RETURN

    END IF

    ! the following is not optimized for multiple processes !

	OPEN(unit, file = 'ParticleDict.txt', status = 'OLD', action = 'READ') ! can file be opened by more than 1 process at the same time?

	READ(unit, fmt = *) npart

    npart = MIN(npart, max_npart)

    ALLOCATE(ipart_arr(npart))
    ALLOCATE(p_igrid_arr(npart))
    ALLOCATE(x(npart))
    ALLOCATE(y(npart))
    ALLOCATE(z(npart))

    SELECT CASE (TRIM(particle_terminal))
        CASE ("none")
            CONTINUE
        CASE ("normal")
            WRITE(*,*) ' '
            WRITE(*, '("Reading ", I0, "Particles:")') npart
        CASE ("verbose")
            WRITE(*,*) ' '
            WRITE(*, '("Reading ", I0, "Particles:")') npart
    END SELECT

    DO ipart = 1, npart

        READ(unit, fmt = *) xtemp, ytemp, ztemp

        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            IF (xtemp < minx) THEN
                CYCLE
            END IF

            IF (xtemp > maxx) THEN
                CYCLE
            END IF

            IF (ytemp < miny) THEN
                CYCLE
            END IF

            IF (ytemp > maxy) THEN
                CYCLE
            END IF

            IF (ztemp < minz) THEN
                CYCLE
            END IF

            IF (ztemp > maxz) THEN
                CYCLE
            END IF

            ipart_arr(ipart) = ipart
            p_igrid_arr(ipart) = igrid
            x(ipart) = xtemp
            y(ipart) = ytemp
            z(ipart) = ztemp

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("Particle read: ID = ", I0, " x/y/z = ", 3F6.2)') ipart, xtemp, ytemp, ztemp
            END SELECT

        END DO

    END DO

    SELECT CASE (TRIM(particle_terminal))
        CASE ("none")
            CONTINUE
        CASE ("normal")
            WRITE(*, *) ' '
            WRITE(*, *) "Reading Particles finished successfully."
        CASE ("verbose")
            WRITE(*, *) ' '
            WRITE(*, *) "Reading Particles finished successfully."
    END SELECT

    END SUBROUTINE read_particles

END MODULE particle_dict_mod
