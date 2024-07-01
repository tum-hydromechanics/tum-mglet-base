MODULE particle_dict_mod

    USE core_mod

CONTAINS

    SUBROUTINE read_particles(dread_particles, npart, ipart_arr, p_igrid_arr, x, y, z)

    !subroutine arguments
    LOGICAL, INTENT(inout) :: dread_particles
    INTEGER(intk), INTENT(in) :: npart
    INTEGER(intk), INTENT(out) :: ipart_arr(npart)
    INTEGER(intk), INTENT(out) :: p_igrid_arr(npart)
    REAL(realk), INTENT(out) :: x(npart), y(npart), z(npart)

    !local variables
    INTEGER(intk) :: unit, n, i, ipart, igrid
    REAL(realk) :: xtemp, ytemp, ztemp
    REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

    INQUIRE(file = 'ParticleDict.txt', exist = dread_particles)

    IF (.NOT. dread_particles) THEN

        !WRITE(*,*) 'No file for reading detected! Using automated initial particle distribution instead.'

        RETURN

    END IF

    ! the following is not optimized for multiple processes !

	OPEN(unit, file = 'ParticleDict.txt', status = 'OLD', action = 'READ') ! can file be opened by more than 1 process at the same time?

	READ(unit, fmt = *) n

    DO ipart = 1, MIN(n, npart)

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

            WRITE(*,*) 'Particle ', ipart, 'read: ', 'x = ', xtemp, 'y = ', ytemp, 'z = ', ztemp

        END DO

    END DO

    END SUBROUTINE read_particles

END MODULE particle_dict_mod
