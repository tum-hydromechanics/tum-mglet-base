! this file is for TYPE definitions and general handling of particles

MODULE particlecore_mod

    !===================================

    USE core_mod
    USE particle_config_mod

    IMPLICIT NONE

    !-----------------------------------

    TYPE :: baseparticle_t ! could be extended by a particle type that includes the mass

        LOGICAL :: is_active = .FALSE. ! if a list of particles is handled, this can be used to identify if an entry is an actual particle

        INTEGER(intk) :: ipart   ! particle ID in global scope (of all processes)
        INTEGER(intk) :: iproc   ! id of process that currently handles the particle
        INTEGER(intk) :: igrid   ! stores the index of grid, which the particle is currently on

        INTEGER(intk) :: ijkcell(3) ! stores indices of the PRESSURE grid cell the particle is in

        REAL(realk) :: x, y, z

        ! time via timekeeper

        CONTAINS

        PROCEDURE :: init
        PROCEDURE :: get_p_igrid
        PROCEDURE :: get_p_ijkcell ! RENAME ?
        PROCEDURE :: update_p_ijkcell ! how to get one pointer for get_p_ijkcell and update_p_ijkcell ? RENAME ?
        PROCEDURE :: print_status

    END TYPE baseparticle_t

    !===================================

CONTAINS

    SUBROUTINE init(this, ipart, x, y, z, iproc, igrid, ijkcell)

        ! Subroutine arguments
        CLASS(baseparticle_t), INTENT(out) :: this

        INTEGER(intk), INTENT(in) :: ipart
        REAL(realk), INTENT(in) :: x, y, z

        INTEGER(intk), INTENT(in), OPTIONAL :: iproc
        INTEGER(intk), INTENT(in), OPTIONAL :: igrid
        INTEGER(intk), INTENT(in), OPTIONAL :: ijkcell(3)

        this%is_active = .TRUE.
        this%ipart = ipart
        this%iproc = myid
        this%x = x
        this%y = y
        this%z = z

        IF (PRESENT(iproc)) THEN

            this%iproc = iproc

        ELSE

            this%iproc = myid

        END IF

        IF (PRESENT(igrid)) THEN

            this%igrid = igrid

        ELSE

            CALL this%get_p_igrid()

        END IF

        IF (PRESENT(ijkcell)) THEN

            this%ijkcell = ijkcell

        ELSEIF (this%is_active) THEN

            CALL this%get_p_ijkcell()

        END IF

    END SUBROUTINE init

    !-----------------------------------

    SUBROUTINE get_p_igrid(this)

        ! Subroutine arguments
        CLASS(baseparticle_t), INTENT(inout) :: this

        !local variables:
        LOGICAL :: found = .FALSE.
        INTEGER(intk) :: i, igrid
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            IF (this%x < minx) THEN
                CYCLE
            END IF

            IF (this%x > maxx) THEN
                CYCLE
            END IF

            IF (this%y < miny) THEN
                CYCLE
            END IF

            IF (this%y > maxy) THEN
                CYCLE
            END IF

            IF (this%z < minz) THEN
                CYCLE
            END IF

            IF (this%z > maxz) THEN
                CYCLE
            END IF

            found = .TRUE.
            this%igrid = igrid

            EXIT

        END DO

        this%is_active = found ! if particle is not found, deactivate particle

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                WRITE(*,*) ' '
                WRITE(*, '("WARNING: In get_p_igrid: Particle ", I0, " could not be found on any grid, process ", I0, " owns!")') this%ipart, this%iproc
        END SELECT

    END SUBROUTINE get_p_igrid

    !-----------------------------------

    SUBROUTINE get_p_ijkcell(this)

        ! subroutine arguments
        CLASS(baseparticle_t), INTENT(inout) :: this

        ! local variables
        INTEGER(intk) :: k, j, i, kk, jj, ii

        TYPE(field_t), POINTER :: x_f, y_f, z_f
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)

        REAL(realk) :: diff_old, diff_new
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz


        IF (.NOT. this%is_active) THEN

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

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, this%igrid)

        ! check if particle is located on igrid, KEEP this check up?

        IF (this%x < minx .OR. this%x > maxx .OR. this%y < miny .OR. this%y > maxy .OR. this%z < minz .OR. this%z > maxz) THEN
            CALL this%get_p_igrid()
        END IF


        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL x_f%get_ptr(x, this%igrid)
        CALL y_f%get_ptr(y, this%igrid)
        CALL z_f%get_ptr(z, this%igrid)

        CALL get_mgdims(kk, jj, ii, this%igrid)

        ! the following assumes that the x/y/z values are sorted such that for any i < j and any direction x_k, x_k(i) < x_k(j) !
        ! the following procedure is capable of handling stretched grids!

        ! find nearest x(i):

        i = 1 + NINT((ii - 1) * (this%x - minx) / (maxx - minx), intk) ! this expression avoids errors for thisx = minx and thisx = maxx

        diff_old = ABS(x(i) - this%x)
        diff_new = 0_realk

        DO WHILE (diff_new < diff_old)

            this%ijkcell(1) = i

            diff_old = ABS(x(i) - this%x)

            IF (x(i) <= this%x) THEN

                i = i + CEILING((ii - i) * (this%x - x(i)) / (maxx - x(i)), intk)

            ELSEIF (x(i) > this%x) THEN

                i = 1_intk + FLOOR(i * (this%x - minx) / (x(i) - minx), intk)

            ELSE

                EXIT

            END IF

            diff_new = ABS(x(i) - this%x)

        END DO

        ! find nearest y(j):

        j = 1 + NINT((jj - 1) * (this%y - miny) / (maxy - miny), intk) ! this expression avoids errors for thisy = miny and thisy = maxy

        diff_old = ABS(y(j) - this%y)
        diff_new = 0_realk

        DO WHILE (diff_new < diff_old)

            this%ijkcell(2) = j

            diff_old = ABS(y(j) - this%y)

            IF (y(j) <= this%y) THEN

                j = j + CEILING((jj - j) * (this%y - y(j)) / (maxy - y(j)), intk)

            ELSEIF (y(j) > this%y) THEN

                j = 1_intk + FLOOR(j * (this%y - miny) / (y(j) - miny), intk)

            ELSE

                EXIT

            END IF

            diff_new = ABS(y(j) - this%y)

        END DO

        ! find nearest z(k):

        k = 1 + NINT((kk - 1) * (this%z - minz) / (maxz - minz), intk) ! this expression avoids errors for thisz = minz and thisz = maxz

        diff_old = ABS(z(k) - this%z)
        diff_new = 0_realk

        DO WHILE (diff_new < diff_old)

            this%ijkcell(3) = k

            diff_old = ABS(z(k) - this%z)

            IF (z(k) <= this%z) THEN

                k = k + CEILING((kk - k) * (this%z - z(k)) / (maxz - z(k)), intk)

            ELSEIF (z(k) > this%z) THEN

                k = 1_intk + FLOOR(k * (this%z - minz) / (z(k) - minz), intk)

            ELSE

                EXIT

            END IF

            diff_new = ABS(z(k) - this%z)

        END DO

    END SUBROUTINE get_p_ijkcell

    !-----------------------------------

    SUBROUTINE update_p_ijkcell(this, pdx, pdy, pdz)

        ! subroutine arguments
        CLASS(baseparticle_t), INTENT(inout) :: this
        REAL(realk), INTENT(in) :: pdx, pdy, pdz

        ! local variables
        TYPE(field_t), POINTER :: x_f, y_f, z_f, dx_f, dy_f, dz_f
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:), dx(:), dy(:), dz(:)

        REAL(realk) :: diff_old, diff_new
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        INTEGER(intk) :: k, j, i, kk, jj, ii
        INTEGER(intk) :: istart, iend, istep, jstart, jend, jstep, kstart, kend, kstep

        IF (.NOT. this%is_active) THEN

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

        CALL x_f%get_ptr(x, this%igrid)
        CALL y_f%get_ptr(y, this%igrid)
        CALL z_f%get_ptr(z, this%igrid)

        CALL dx_f%get_ptr(dx, this%igrid)
        CALL dy_f%get_ptr(dy, this%igrid)
        CALL dz_f%get_ptr(dz, this%igrid)

        CALL get_mgdims(kk, jj, ii, this%igrid)

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, this%igrid)

        ! the following assumes that the x/y/z values are sorted such that for any i < j and any direction x, x(i) < x(j) !
        ! the following procedure is capable of handling stretched grids!

        ! find nearest x:
        istart = this%ijkcell(1) + NINT((this%x - x(this%ijkcell(1))) / dx(this%ijkcell(1)))
        istart = MAX(istart, 1_intk)
        istart = MIN(istart, ii)

        istep = SIGN(INT(1, intk), NINT((this%x - x(istart)), intk))
        iend = 1_intk + (istep + 1_intk) / 2_intk * (ii - 1_intk)

        diff_old = x(istart) - this%x

        DO i = istart + istep, iend, istep

            diff_new = x(i) - this%x

            IF (diff_new > diff_old) THEN

                this%ijkcell(1) = i - istep
                EXIT

            ELSEIF (i == iend) THEN

                this%ijkcell(1) = i

            ELSE

                diff_old = diff_new

            END IF
        END DO

        ! find nearest y:

        jstart = this%ijkcell(2) + NINT((this%y - y(this%ijkcell(2))) / dy(this%ijkcell(2)))
        jstart = MAX(jstart, 1_intk)
        jstart = MIN(jstart, jj)

        jstep = SIGN(INT(1, intk), NINT((this%y - y(jstart)), intk))
        jend = 1_intk + (jstep + 1_intk) / 2_intk * (jj - 1_intk)

        diff_old = y(jstart) - this%y

        DO j = jstart + jstep, jend, jstep

            diff_new = y(j) - this%y

            IF (diff_new > diff_old) THEN

                this%ijkcell(2) = j - jstep
                EXIT

            ELSEIF (j == jend) THEN

                this%ijkcell(2) = j

            ELSE
                diff_old = diff_new
                CYCLE

            END IF
        END DO

        ! find nearest z:

        kstart = this%ijkcell(3) + NINT((this%z - z(this%ijkcell(3))) / dz(this%ijkcell(3)))
        kstart = MAX(kstart, 1_intk)
        kstart = MIN(kstart, kk)

        kstep = SIGN(INT(1, intk), NINT((this%z - z(kstart)), intk))
        kend = 1_intk + (kstep + 1_intk) / 2_intk * (kk - 1_intk)

        diff_old = z(kstart) - this%z

        DO k = kstart + kstep, kend, kstep

            diff_new = z(k) - this%z

            IF (diff_new > diff_old) THEN

                this%ijkcell(3) = k - kstep
                EXIT

            ELSEIF (k == kend) THEN

                this%ijkcell(3) = k

            ELSE

                diff_old = diff_new
                CYCLE

            END IF
        END DO

    END SUBROUTINE update_p_ijkcell

    SUBROUTINE print_status(this, time)

    ! subroutine arguments
    CLASS(baseparticle_t), INTENT(inout) :: this
    REAL(realk) :: time

        IF (this%is_active) THEN

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("Particle ", I0, " - Status at t = ", 3F12.6, ":")') this%ipart, time
                    WRITE(*,'("iproc       = ", I12,   " igrid = ", I12)') this%iproc, this%igrid
                    WRITE(*,'("x/y/z       = ", F12.6, " / ", F12.6, " / ", F12.6)') this%x, this%y, this%z
                    WRITE(*,'("i/j/k cell  = ", I12  , " / ", I12  , " / ", I12)') this%ijkcell(1), this%ijkcell(2), this%ijkcell(3)
            END SELECT

        ELSE

            RETURN

        END IF

    END SUBROUTINE print_status

    !===================================

END MODULE particlecore_mod
