! particle file is for TYPE definitions and general handling of particles

MODULE particle_core_mod

    USE grids_mod
    USE field_mod
    USE fields_mod


    USE particle_config_mod

    IMPLICIT NONE

    INTEGER(c_intk), PARAMETER :: particle_mpi_elems = 8

    TYPE, BIND(C) :: baseparticle_t

        ! state definitions:
        !-1/0 : particle is inactive
        !  1 : particle is active
        !  2 : particle is active, cell info has to be updated after displacement
        !  3 : particle is active, grid and cell info have to be updated after displacement
        !  4 : particle is active, particle has to be passed to other process after displacement
        INTEGER(c_intk) :: state = -1

        INTEGER(c_intk) :: ipart = -1
        INTEGER(c_intk) :: iproc = -1
        INTEGER(c_intk) :: igrid = -1

        INTEGER(c_intk) :: ijkcell(3)

        REAL(c_realk) :: x = 0.0
        REAL(c_realk) :: y = 0.0
        REAL(c_realk) :: z = 0.0

    END TYPE baseparticle_t

    PUBLIC :: set_particle, set_particle_igrid, set_particle_cell, &
              update_particle_cell, print_particle_status

CONTAINS

    SUBROUTINE set_particle(particle, ipart, x, y, z, iproc, igrid, ijkcell)

        ! Subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle

        INTEGER(intk), INTENT(in) :: ipart
        REAL(realk), INTENT(in) :: x, y, z
        INTEGER(intk), INTENT(in), OPTIONAL :: iproc
        INTEGER(intk), INTENT(in), OPTIONAL :: igrid
        INTEGER(intk), INTENT(in), OPTIONAL :: ijkcell(3)

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

    END SUBROUTINE set_particle



    SUBROUTINE set_particle_igrid(particle)

        ! Subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle

        !local variables:
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
                    WRITE(*,*) ' '
                    WRITE(*, *) "WARNING: In get_p_ijkcell: Tried to locate particle that is not active!"
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

        particle%state = found ! if particle is not found, deactivate particle

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                WRITE(*,*) ' '
                WRITE(*, '("WARNING: In get_p_igrid: Particle ", I0, " could not be found on any grid, process ", I0, " owns!")') particle%ipart, particle%iproc
        END SELECT

    END SUBROUTINE set_particle_igrid



    SUBROUTINE set_particle_cell(particle)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle

        ! local variables
        INTEGER(intk) :: k, j, i, kk, jj, ii

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
                    WRITE(*, '("WARNING: Coordinates of particle ", I0 ," do not lie within boundaries of particle grid ", I0)') particle%ipart, particle%igrid
            END SELECT

            CALL set_particle_igrid( particle )
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
        i = 1 + NINT((ii - 1) * (particle%x - minx) / (maxx - minx), intk) ! particle expression avoids errors for particlex = minx and particlex = maxx
        diff_old = ABS(x(i) - particle%x)
        diff_new = 0_realk

        DO WHILE (diff_new < diff_old)
            particle%ijkcell(1) = i
            diff_old = ABS(x(i) - particle%x)
            IF (x(i) <= particle%x) THEN
                i = i + CEILING((ii - i) * (particle%x - x(i)) / (maxx - x(i)), intk)
            ELSEIF (x(i) > particle%x) THEN
                i = 1_intk + FLOOR(i * (particle%x - minx) / (x(i) - minx), intk)
            ELSE
                EXIT
            END IF
            diff_new = ABS(x(i) - particle%x)
        END DO

        ! find nearest y(j):
        j = 1 + NINT((jj - 1) * (particle%y - miny) / (maxy - miny), intk) ! particle expression avoids errors for particley = miny and particley = maxy
        diff_old = ABS(y(j) - particle%y)
        diff_new = 0_realk

        DO WHILE (diff_new < diff_old)
            particle%ijkcell(2) = j
            diff_old = ABS(y(j) - particle%y)
            IF (y(j) <= particle%y) THEN
                j = j + CEILING((jj - j) * (particle%y - y(j)) / (maxy - y(j)), intk)
            ELSEIF (y(j) > particle%y) THEN
                j = 1_intk + FLOOR(j * (particle%y - miny) / (y(j) - miny), intk)
            ELSE
                EXIT
            END IF
            diff_new = ABS(y(j) - particle%y)
        END DO

        ! find nearest z(k):
        k = 1 + NINT((kk - 1) * (particle%z - minz) / (maxz - minz), intk) ! particle expression avoids errors for particlez = minz and particlez = maxz
        diff_old = ABS(z(k) - particle%z)
        diff_new = 0_realk

        DO WHILE (diff_new < diff_old)
            particle%ijkcell(3) = k
            diff_old = ABS(z(k) - particle%z)
            IF (z(k) <= particle%z) THEN
                k = k + CEILING((kk - k) * (particle%z - z(k)) / (maxz - z(k)), intk)
            ELSEIF (z(k) > particle%z) THEN
                k = 1_intk + FLOOR(k * (particle%z - minz) / (z(k) - minz), intk)
            ELSE
                EXIT
            END IF
            diff_new = ABS(z(k) - particle%z)
        END DO

    END SUBROUTINE set_particle_cell



    SUBROUTINE update_particle_cell(particle)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle

        ! local variables
        TYPE(field_t), POINTER :: x_f, y_f, z_f, dx_f, dy_f, dz_f
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:), dx(:), dy(:), dz(:)

        REAL(realk) :: diff_old, diff_new
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        INTEGER(intk) :: k, j, i, kk, jj, ii
        INTEGER(intk) :: istart, iend, istep, jstart, jend, jstep, kstart, kend, kstep

        IF (particle%state < 1) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,*) ' '
                    WRITE(*, *) "WARNING: In update_p_ijkcell: Tried to locate particle that is not active!"
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

        ! the following assumes that the x/y/z values are sorted such that for any i < j and any direction x, x(i) < x(j) !
        ! the following procedure is capable of handling stretched grids!

        ! find nearest x:
        istart = particle%ijkcell(1) + NINT((particle%x - x(particle%ijkcell(1))) / dx(particle%ijkcell(1)))
        istart = MAX(istart, 1_intk)
        istart = MIN(istart, ii)
        istep = SIGN(INT(1, intk), NINT((particle%x - x(istart)), intk))
        iend = 1_intk + (istep + 1_intk) / 2_intk * (ii - 1_intk)

        diff_old = x(istart) - particle%x

        DO i = istart + istep, iend, istep
            diff_new = x(i) - particle%x
            IF (diff_new > diff_old) THEN
                particle%ijkcell(1) = i - istep
                EXIT
            ELSEIF (i == iend) THEN
                particle%ijkcell(1) = i
            ELSE
                diff_old = diff_new
            END IF
        END DO

        ! find nearest y:
        jstart = particle%ijkcell(2) + NINT((particle%y - y(particle%ijkcell(2))) / dy(particle%ijkcell(2)))
        jstart = MAX(jstart, 1_intk)
        jstart = MIN(jstart, jj)
        jstep = SIGN(INT(1, intk), NINT((particle%y - y(jstart)), intk))
        jend = 1_intk + (jstep + 1_intk) / 2_intk * (jj - 1_intk)

        diff_old = y(jstart) - particle%y

        DO j = jstart + jstep, jend, jstep
            diff_new = y(j) - particle%y
            IF (diff_new > diff_old) THEN
                particle%ijkcell(2) = j - jstep
                EXIT
            ELSEIF (j == jend) THEN
                particle%ijkcell(2) = j
            ELSE
                diff_old = diff_new
                CYCLE
            END IF
        END DO

        ! find nearest z:
        kstart = particle%ijkcell(3) + NINT((particle%z - z(particle%ijkcell(3))) / dz(particle%ijkcell(3)))
        kstart = MAX(kstart, 1_intk)
        kstart = MIN(kstart, kk)
        kstep = SIGN(INT(1, intk), NINT((particle%z - z(kstart)), intk))
        kend = 1_intk + (kstep + 1_intk) / 2_intk * (kk - 1_intk)

        diff_old = z(kstart) - particle%z

        DO k = kstart + kstep, kend, kstep
            diff_new = z(k) - particle%z
            IF (diff_new > diff_old) THEN
                particle%ijkcell(3) = k - kstep
                EXIT
            ELSEIF (k == kend) THEN
                particle%ijkcell(3) = k
            ELSE
                diff_old = diff_new
                CYCLE
            END IF
        END DO

    END SUBROUTINE update_particle_cell



    SUBROUTINE print_particle_status(particle)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(in) :: particle

        IF (particle%state >= 1) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, '("Particle ", I0, " - Status:")') particle%ipart
                    WRITE(*, '("iproc       = ", I12)') particle%iproc
                    WRITE(*, '("igrid       = ", I12)') particle%igrid
                    WRITE(*, '("x/y/z       = ", 3F12.6)') particle%x, particle%y, particle%z
                    WRITE(*, '("i/j/k cell  = ", 3I12)') particle%ijkcell(1), particle%ijkcell(2), particle%ijkcell(3)
                    !WRITE(*, '("facepath    = ", 3I12)') particle%facepath(1), particle%facepath(2), particle%facepath(3)
                    WRITE(*, *) ' '
            END SELECT
        ELSE
            RETURN
        END IF

    END SUBROUTINE print_particle_status

END MODULE particle_core_mod
