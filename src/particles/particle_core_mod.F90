MODULE particle_core_mod

    ! This module is responsible for:
    ! Definition of the baseparticle_t type.
    ! Basic operations on individual particles.

    USE grids_mod
    USE field_mod
    USE fields_mod

    USE particle_config_mod

    IMPLICIT NONE

    INTEGER(c_intk), PARAMETER :: particle_mpi_elems = 11

    ! C binding for MPI compatability!
    TYPE, BIND(C) :: baseparticle_t

        ! Particle state definitions:
        ! ...-2/-1/0 : particle is inactive
        !  1 : particle is active, no further information specified
        !  2 : particle is active, cell info has to be updated after displacement
        !  3 : particle is active, grid and cell info have to be updated after displacement
        !  4 : particle is active, particle has to be passed to other process after displacement
        INTEGER(c_intk) :: state = -1

        INTEGER(c_intk) :: ipart = -1
        INTEGER(c_intk) :: iproc = -1
        INTEGER(c_intk) :: igrid = -1
        INTEGER(c_intk) :: islice = -1

        !gitstep is the timestep at which a particle entered its current grid (for residence time tracking)
        INTEGER(c_intk) :: gitstep = 0
        !sitstep is the timestep at which a particle entered its current slice (for residence time tracking)
        INTEGER(c_intk) :: sitstep = 0

        INTEGER(c_intk) :: ijkcell(3)

        REAL(c_realk) :: x = 0.0
        REAL(c_realk) :: y = 0.0
        REAL(c_realk) :: z = 0.0

    END TYPE baseparticle_t

    PUBLIC :: set_particle, set_particle_igrid, set_particle_cell, &
              update_particle_cell, print_particle_status

CONTAINS    !===================================

    SUBROUTINE set_particle(particle, ipart, x, y, z, iproc, igrid, ijkcell, gitstep, sitstep)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: ipart
        REAL(realk), INTENT(in) :: x, y, z
        INTEGER(intk), INTENT(in), OPTIONAL :: iproc
        INTEGER(intk), INTENT(in), OPTIONAL :: igrid
        INTEGER(intk), INTENT(in), OPTIONAL :: ijkcell(3)
        INTEGER(intk), INTENT(in), OPTIONAL :: gitstep, sitstep

        particle%state = 1
        particle%ipart = ipart
        particle%iproc = myid
        particle%x = x
        particle%y = y
        particle%z = z

        IF (PRESENT(iproc)) THEN
            particle%iproc = iproc
        ELSE
            particle%iproc = myid
        END IF

        IF (PRESENT(igrid)) THEN
            particle%igrid = igrid
        ELSE
            CALL set_particle_igrid(particle)
        END IF

        IF (PRESENT(ijkcell)) THEN
            particle%ijkcell = ijkcell
        ELSEIF (particle%state >= 1) THEN
            CALL set_particle_cell(particle)
        END IF

        IF (PRESENT(gitstep)) THEN
            particle%gitstep = gitstep
        ELSE
            particle%gitstep = 0
        END IF

        IF (PRESENT(sitstep)) THEN
            particle%sitstep = sitstep
        ELSE
            particle%sitstep = 0
        END IF

    END SUBROUTINE set_particle

    !-----------------------------------

    SUBROUTINE set_particle_igrid(particle)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle

        ! local variables:
        INTEGER(intk) :: found = -1
        INTEGER(intk) :: i, igrid
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        IF (particle%state < 1) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, *) "WARNING: In get_p_ijkcell: Tried to locate particle that is not active!"
                    WRITE(*, '()')
            END SELECT
            RETURN
        END IF

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
            IF (particle%x < minx) THEN
                CYCLE
            END IF
            IF (particle%x > maxx) THEN
                CYCLE
            END IF
            IF (particle%y < miny) THEN
                CYCLE
            END IF
            IF (particle%y > maxy) THEN
                CYCLE
            END IF
            IF (particle%z < minz) THEN
                CYCLE
            END IF
            IF (particle%z > maxz) THEN
                CYCLE
            END IF
            found = 1
            particle%igrid = igrid
            EXIT
        END DO

        ! if particle is not found, deactivate particle
        particle%state = found

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                WRITE(*, '("WARNING: In get_p_igrid: Particle ", I0, " could not be found on any grid, process ", I0, " owns!")') &
                 particle%ipart, particle%iproc
                WRITE(*, '()')
        END SELECT

    END SUBROUTINE set_particle_igrid

    !-----------------------------------

    SUBROUTINE set_particle_cell(particle)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle

        ! local variables
        INTEGER(intk) :: k, j, i, kk, jj, ii, counter, max_iterations = 10

        TYPE(field_t), POINTER :: x_f, y_f, z_f
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk) :: diff_old, diff_new
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        IF (particle%state < 1) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,*) ' '
                    WRITE(*, *) "WARNING: In get_p_ijkcell: Tried to locate particle that is not active!"
            END SELECT

            RETURN

        END IF

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, particle%igrid)

        ! check if particle is located on igrid, KEEP particle check up?
        IF (particle%x < minx .OR. particle%x > maxx .OR. particle%y < miny .OR. particle%y > maxy .OR. particle%z < minz .OR. particle%z > maxz) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,*) ' '
                    WRITE(*, '("WARNING: Coordinates of particle ", I0 ," do not lie within boundaries of particle grid ", I0)') &
                     particle%ipart, particle%igrid
            END SELECT

            CALL set_particle_igrid(particle)

        END IF

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")
        CALL x_f%get_ptr(x, particle%igrid)
        CALL y_f%get_ptr(y, particle%igrid)
        CALL z_f%get_ptr(z, particle%igrid)
        CALL get_mgdims(kk, jj, ii, particle%igrid)

        ! the following assumes that the x/y/z values are sorted such that for any i < j and any direction x_k, x_k(i) < x_k(j) !
        ! the following procedure is capable of handling stretched grids!

        ! find nearest x(i):
        counter = 1

        i = 3 + NINT((ii - 5) * (particle%x - minx) / (maxx - minx), intk)
        particle%ijkcell(1) = i
        diff_old = ABS(x(i) - particle%x)
        diff_new = 0

        DO WHILE (diff_new < diff_old .AND. counter <= max_iterations)
            particle%ijkcell(1) = i
            diff_old = ABS(x(i) - particle%x)
            IF (x(i) <= particle%x) THEN
                ! the denominator of the fraction will NOT be zero because x(i) /= maxx for all i
                i = i + CEILING((ii - 2 - i) * (particle%x - x(i)) / (maxx - x(i)), intk)
            ELSEIF (x(i) > particle%x) THEN
                i = 3 + FLOOR((i - 2) * (particle%x - minx) / (x(i) - minx), intk)
            ELSE
                EXIT
            END IF
            diff_new = ABS(x(i) - particle%x)

            counter = counter + 1
        END DO

        IF (i < 3) THEN
            i = 3
            particle%ijkcell(1) = i
        ELSEIF (i > (ii - 2)) THEN
            i = ii - 2
            particle%ijkcell(1) = i
        END IF

        ! find nearest y(j):
        counter = 1

        j = 3 + NINT((jj - 5) * (particle%y - miny) / (maxy - miny), intk)
        particle%ijkcell(2) = j
        diff_old = ABS(y(j) - particle%y)
        diff_new = 0

        DO WHILE (diff_new < diff_old .AND. counter <= max_iterations)
            particle%ijkcell(2) = j
            diff_old = ABS(y(j) - particle%y)
            IF (y(j) <= particle%y) THEN
            ! the denominator of the fraction will NOT be zero because y(j) /= maxx for all j
                j = j + CEILING((jj - 2 - j) * (particle%y - y(j)) / (maxy - y(j)), intk)
            ELSEIF (y(j) > particle%y) THEN
                j = 3 + FLOOR((j - 2) * (particle%y - miny) / (y(j) - miny), intk)
            ELSE
                EXIT
            END IF
            diff_new = ABS(y(j) - particle%y)

            counter = counter + 1
        END DO

        IF (j < 3) THEN
            j = 3
            particle%ijkcell(2) = j
        ELSEIF (j > (jj - 2)) THEN
            j = jj - 2
            particle%ijkcell(2) = j
        END IF

        ! find nearest z(k):
        counter = 1

        k = 3 + NINT((kk - 5) * (particle%z - minz) / (maxz - minz), intk)
        particle%ijkcell(3) = k
        diff_old = ABS(z(k) - particle%z)
        diff_new = 0

        DO WHILE (diff_new < diff_old .AND. counter <= max_iterations)
            particle%ijkcell(3) = k
            diff_old = ABS(z(k) - particle%z)
            IF (z(k) <= particle%z) THEN
                ! the denominator of the fraction will NOT be zero because z(k) /= maxx for all k
                k = k + CEILING((kk - 2 - k) * (particle%z - z(k)) / (maxz - z(k)), intk)
            ELSEIF (z(k) > particle%z) THEN
                k = 3 + FLOOR((k - 2) * (particle%z - minz) / (z(k) - minz), intk)
            ELSE
                EXIT
            END IF
            diff_new = ABS(z(k) - particle%z)

            counter = counter + 1
        END DO

        IF (k < 3) THEN
            k = 3
            particle%ijkcell(3) = k
        ELSEIF (i > (ii - 2)) THEN
            k = kk - 2
            particle%ijkcell(3) = k
        END IF

        ! TODO: rethink this safety operation
        CALL update_particle_cell(particle)

    END SUBROUTINE set_particle_cell

    !-----------------------------------

    SUBROUTINE update_particle_cell(particle)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle

        ! local variables
        TYPE(field_t), POINTER :: x_f, y_f, z_f, dx_f, dy_f, dz_f
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:), dx(:), dy(:), dz(:)

        REAL(realk) :: diff_old, diff_new
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        INTEGER(intk) :: k, j, i, kk, jj, ii, istep, jstep, kstep

        IF (particle%state < 1) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, *) "WARNING: In update_p_ijkcell: Tried to locate particle that is not active!"
                    WRITE(*, '()')
            END SELECT
            RETURN
        END IF

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL x_f%get_ptr(x, particle%igrid)
        CALL y_f%get_ptr(y, particle%igrid)
        CALL z_f%get_ptr(z, particle%igrid)

        CALL dx_f%get_ptr(dx, particle%igrid)
        CALL dy_f%get_ptr(dy, particle%igrid)
        CALL dz_f%get_ptr(dz, particle%igrid)

        CALL get_mgdims(kk, jj, ii, particle%igrid)
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, particle%igrid)

        ! the following assumes that the x/y/z grid coordinates are sorted such that for any i < j and any direction x, x(i) < x(j) !
        ! the following procedure is capable of handling stretched grids!

        ! find nearest x:
        istep = INT(SIGN(1.0_realk, particle%x - x(particle%ijkcell(1))), intk)

        i = particle%ijkcell(1) + istep

        diff_old = ABS(x(particle%ijkcell(1)) - particle%x)
        diff_new = ABS(x(i) - particle%x)

        DO WHILE (diff_new < diff_old)
            i = i + istep
            diff_old = diff_new
            diff_new = ABS(x(i) - particle%x)
        END DO

        particle%ijkcell(1) = MIN(MAX(i - istep, 1), ii)

        ! find nearest y:
        jstep = INT(SIGN(1.0_realk, particle%y - y(particle%ijkcell(2))), intk)

        j = particle%ijkcell(2) + jstep

        diff_old = ABS(y(particle%ijkcell(2)) - particle%y)
        diff_new = ABS(y(j) - particle%y)

        DO WHILE (diff_new < diff_old)
            j = j + jstep
            diff_old = diff_new
            diff_new = ABS(y(j) - particle%y)
        END DO

        particle%ijkcell(2) = MIN(MAX(j - jstep, 1), jj)

        ! find nearest z:
        kstep = INT(SIGN(1.0_realk, particle%z - z(particle%ijkcell(3))), intk)

        k = particle%ijkcell(3) + kstep

        diff_old = ABS(z(particle%ijkcell(3)) - particle%z)
        diff_new = ABS(z(k) - particle%z)

        DO WHILE (diff_new < diff_old)
            k = k + kstep
            diff_old = diff_new
            diff_new = ABS(z(k) - particle%z)
        END DO

        particle%ijkcell(3) = MIN(MAX(k - kstep, 1), kk)

    END SUBROUTINE update_particle_cell

    !-----------------------------------

    SUBROUTINE print_particle_status(particle)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(in) :: particle

        IF (particle%state >= 1) THEN
                WRITE(*, '("Particle ", I0, " - Status:")') particle%ipart
                WRITE(*, '("iproc       = ", I20)') particle%iproc
                WRITE(*, '("igrid       = ", I20)') particle%igrid
                WRITE(*, '("islice      = ", I20)') particle%islice
                WRITE(*, '("gitstep     = ", I20)') particle%gitstep
                WRITE(*, '("sitstep     = ", I20)') particle%sitstep
                WRITE(*, '("x/y/z       = ", 3F120.17)') particle%x, particle%y, particle%z
                WRITE(*, '("i/j/k cell  = ", 3I20)') particle%ijkcell(1), particle%ijkcell(2), particle%ijkcell(3)
                WRITE(*, '()')
        ELSE
            RETURN
        END IF

    END SUBROUTINE print_particle_status

END MODULE particle_core_mod
