MODULE particle_diffusion_mod

    USE precision_mod
    USE charfunc_mod
    USE fort7_mod
    USE grids_mod
    USE field_mod
    USE fields_mod
    USE connect2_mod

    USE particle_config_mod
    USE particle_core_mod
    USE particle_interpolation_mod

    IMPLICIT NONE

    ! REFERENCE VALUES FOR THE STANDART NORMAL DISTRIBUTIONS:
    ! realization values
    REAL(intk) :: sn_x(22) = [-3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, &
       1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,  2.8,  3.0]
    ! corresponding cumulative probabilities
    REAL(intk) :: sn_pc(22) = [0.00135, 0.00256, 0.00466, 0.00820, 0.01390, 0.02275, 0.03593, 0.05480, 0.08076, 0.11507, 0.15866, &
      0.84134, 0.88493, 0.91924, 0.94520, 0.96407, 0.97725, 0.98610, 0.99180, 0.99534, 0.99744, 0.99865]
    ! corresponding probabilites
    REAL(intk) :: sn_p(22)

    ! truncation limit stored in config mod
    REAL(realk) :: truncation_factor

CONTAINS

    SUBROUTINE init_particle_diffusion()

        ! local_variables
        INTEGER(intk) :: i
        INTEGER(intk), PARAMETER :: units_diff(7) = [0, 2, -1, 0, 0, 0, 0]

        CALL start_timer(900)
        CALL start_timer(910)

        DO i = 1, SIZE(sn_p)
            sn_p(i) = 1 / SQRT(2 * pi) * EXP(- (sn_x(i)** 2) / 2)
        END DO

        ! given a truncation narrower than APPROXIMATELY [-truncation_limit = -1.75, truncation_limit = 1.75]
        ! of the (stretched) parent pdf, the target standart deviation of 1 cannot be achieved for the truncated pdf
        IF (truncation_limit <= 1.8) THEN
            WRITE(*, *) "Truncation set too narrow. Invalid truncation limit of ", truncation_limit, "."
            WRITE(*, *) "Try a truncation limit above 1.8, or use another random walk mode"
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL get_truncation_factor(truncation_limit, truncation_factor)

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("TRUNCATION CORRECTION FACTOR: ", F12.3)') truncation_factor
                WRITE(*, '()')
            END IF
        END IF

        IF (.NOT. dturb_diff) THEN
            CALL stop_timer(910)
            CALL stop_timer(900)
            RETURN
        END IF

        ! fields that hold the total (trubulent + diffusive) particle diffusion
        ! if trub_diff == .TRUE. fields i/o will try to read the diffusion fields,
        ! if they cannot be read, then generate_turbulent_diffusion will try to generate the diffsuion fields
        ! if diffsuion fields cannot be read nor generated, the simulation is terminated
        CALL set_field("P_DIFF_X", istag = 1, units = units_diff, &
         dread = dcont, required = dcont, dwrite = .TRUE., buffers = .TRUE.)
        CALL set_field("P_DIFF_Y", jstag = 1, units = units_diff, &
         dread = dcont, required = dcont, dwrite = .TRUE., buffers = .TRUE.)
        CALL set_field("P_DIFF_Z", kstag = 1, units = units_diff, &
         dread = dcont, required = dcont, dwrite = .TRUE., buffers = .TRUE.)

        IF (.NOT. dcont) THEN
            CALL generate_diffusion_field
        END IF

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_diffusion

    SUBROUTINE generate_diffusion_field()

        ! local_variables
        INTEGER(intk) :: igrid, ii, jj, kk, g, i, j, k, ilevel

        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz
        TYPE(field_t), POINTER :: t1_avg_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: t1_avg
        TYPE(field_t), POINTER :: u_avg_f, v_avg_f, w_avg_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u_avg, v_avg, w_avg
        TYPE(field_t), POINTER :: ut1_avg_f, vt1_avg_f, wt1_avg_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ut1_avg, vt1_avg, wt1_avg
        TYPE(field_t), POINTER :: diffx_f, diffy_f, diffz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: diffx, diffy, diffz

        ! if one of the following fields does not exist, get_field will call an error and terminate the process
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(t1_avg_f, "T1_AVG")

        CALL get_field(u_avg_f, "U_AVG")
        CALL get_field(v_avg_f, "V_AVG")
        CALL get_field(w_avg_f, "W_AVG")

        CALL get_field(ut1_avg_f, "UT1_AVG")
        CALL get_field(vt1_avg_f, "VT1_AVG")
        CALL get_field(wt1_avg_f, "WT1_AVG")

        CALL get_field(diffx_f, "P_DIFF_X")
        CALL get_field(diffy_f, "P_DIFF_Y")
        CALL get_field(diffz_f, "P_DIFF_Z")

        DO g = 1, nmygrids

            igrid = mygrids(g)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL t1_avg_f%get_ptr(t1_avg, igrid)

            CALL u_avg_f%get_ptr(u_avg, igrid)
            CALL v_avg_f%get_ptr(v_avg, igrid)
            CALL w_avg_f%get_ptr(w_avg, igrid)

            CALL ut1_avg_f%get_ptr(ut1_avg, igrid)
            CALL vt1_avg_f%get_ptr(vt1_avg, igrid)
            CALL wt1_avg_f%get_ptr(wt1_avg, igrid)

            CALL diffx_f%get_ptr(diffx, igrid)
            CALL diffy_f%get_ptr(diffy, igrid)
            CALL diffz_f%get_ptr(diffz, igrid)

            DO i = 3, ii - 2
                DO j = 3, jj - 2
                    DO k = 3, kk - 2
                        diffx(k, j, i) = (ut1_avg(k, j, i) - u_avg(k, j, i) * 0.5 * (t1_avg(k, j, i + 1) + t1_avg(k, j, i))) &
                         * dx(i) / (t1_avg(k, j, i + 1) - t1_avg(k, j, i)) + D(1)
                    END DO
                END DO
            END DO

            DO i = 3, ii - 2
                DO j = 3, jj - 2
                    DO k = 3, kk - 2
                        diffy(k, j, i) = (vt1_avg(k, j, i) - v_avg(k, j, i) * 0.5 * (t1_avg(k, j + 1, i) + t1_avg(k, j, i))) &
                         * dy(j) / (t1_avg(k, j + 1, i) - t1_avg(k, j, i)) + D(2)
                    END DO
                END DO
            END DO

            DO i = 3, ii - 2
                DO j = 3, jj - 2
                    DO k = 3, kk - 2
                        diffz(k, j, i) = (wt1_avg(k, j, i) - w_avg(k, j, i) * 0.5 * (t1_avg(k + 1, j, i) + t1_avg(k, j, i))) &
                         * dz(k) / (t1_avg(k + 1, j, i) - t1_avg(k, j, i)) + D(3)
                    END DO
                END DO
            END DO

            ! TODO: check this approach for buffers; esp. levels ...
            DO ilevel = minlevel, maxlevel
                CALL connect(ilevel, 1, v1 = diffx_f, v2 = diffy_f, v3 = diffz_f, corners = .TRUE.)
            END DO
        END DO

    END SUBROUTINE generate_diffusion_field

    SUBROUTINE get_pdiff(particle, dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff, dinterp_pdiffsuion)

            ! subroutine arguments
            TYPE(baseparticle_t), INTENT(in) :: particle
            REAL(realk), INTENT(in) :: dt
            LOGICAL, INTENT(in) :: dinterp_pdiffsuion
            REAL(realk), INTENT(out) :: pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff

            ! local variables
            INTEGER(intk) :: kk, jj, ii
            TYPE(field_t), POINTER :: x_f, y_f, z_f
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
            TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz
            TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx, ddy, ddz
            TYPE(field_t), POINTER :: diffx_f, diffy_f, diffz_f
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: diffx, diffy, diffz

            ! local variables
            REAL(realk) :: p_diffx, p_diffy, p_diffz, ranx, rany, ranz
            LOGICAL :: dusefield, dinterp

            IF (.NOT. ddiffusion) THEN
                RETURN
            END IF

            CALL get_mgdims(kk, jj, ii, particle%igrid)

            IF (dturb_diff) THEN

                CALL get_field(x_f, "X")
                CALL get_field(y_f, "Y")
                CALL get_field(z_f, "Z")

                CALL x_f%get_ptr(x, particle%igrid)
                CALL y_f%get_ptr(y, particle%igrid)
                CALL z_f%get_ptr(z, particle%igrid)

                CALL get_field(diffx_f, "P_DIFF_X")
                CALL get_field(diffy_f, "P_DIFF_Y")
                CALL get_field(diffz_f, "P_DIFF_Z")

                CALL diffx_f%get_ptr(diffx, particle%igrid)
                CALL diffy_f%get_ptr(diffy, particle%igrid)
                CALL diffz_f%get_ptr(diffz, particle%igrid)

                IF (dinterp_pdiffsuion) THEN

                    CALL get_field(dx_f, "DX")
                    CALL get_field(dy_f, "DY")
                    CALL get_field(dz_f, "DZ")

                    CALL dx_f%get_ptr(dx, particle%igrid)
                    CALL dy_f%get_ptr(dy, particle%igrid)
                    CALL dz_f%get_ptr(dz, particle%igrid)

                    CALL get_field(ddx_f, "DDX")
                    CALL get_field(ddy_f, "DDY")
                    CALL get_field(ddz_f, "DDZ")

                    CALL ddx_f%get_ptr(ddx, particle%igrid)
                    CALL ddy_f%get_ptr(ddy, particle%igrid)
                    CALL ddz_f%get_ptr(ddz, particle%igrid)

                    CALL interpolate_lincon(particle, kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                     diffx, diffy, diffz, p_diffx, p_diffy, p_diffz)

                ELSE

                    CALL get_nearest_value(particle, kk, jj, ii, x, y, z, &
                    diffx, diffy, diffz, p_diffx, p_diffy, p_diffz)

                END IF

            ELSE

                P_diffx = D(1)
                p_diffy = D(2)
                p_diffz = D(3)

            END IF

    END SUBROUTINE get_pdiff

    SUBROUTINE generate_diffusive_displacement(dt, D_x, D_y, D_z, pdx, pdy, pdz)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: D_x, D_y, D_z
        REAL(realk), INTENT(out) :: pdx, pdy, pdz

        ! local variables
        REAL(realk) :: sigx, sigy, sigz, ranx, rany, ranz

        IF (D_x > 0) THEN
            sigx = SQRT(2 * D_x * dt)

            SELECT CASE (lower(TRIM(random_walk_mode)))
            CASE ("rademacher")
                CALL rademacher_dist(sigx, ranx)
            CASE ("uniform")
                CALL uniform_dist(sigx, ranx)
            CASE ("gaussian2")
                CALL gaussian_dist2(0.0_realk, sigx, ranx)
            END SELECT

            pdx = ranx ! diffusion length

        END IF

        IF (D_y > 0) THEN

            sigy = SQRT(2 * D_y * dt)

            SELECT CASE (lower(TRIM(random_walk_mode)))
            CASE ("rademacher")
                CALL rademacher_dist(sigy, rany)
            CASE ("uniform")
                CALL uniform_dist(sigy, rany)
            CASE ("gaussian2")
                CALL gaussian_dist2(0.0_realk, sigy, rany)
            END SELECT

            pdy = rany ! diffusion length

        END IF

        IF (D_z > 0) THEN

            sigz = SQRT(2 * D_z * dt)

            SELECT CASE (lower(TRIM(random_walk_mode)))
            CASE ("rademacher")
                CALL rademacher_dist(sigz, ranz)
            CASE ("uniform")
                CALL uniform_dist(sigz, ranz)
            CASE ("gaussian2")
                CALL gaussian_dist2(0.0_realk, sigz, ranz)
            END SELECT

            pdz = ranz ! diffusion length

        END IF

    END SUBROUTINE generate_diffusive_displacement

    SUBROUTINE rademacher_dist(sigma, R)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: sigma
        REAL(realk), INTENT(out) :: R

        R = 0.0
        !CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(R)
        R = R - 0.5_realk
        R = SIGN(1.0_realk, R) * sigma

    END SUBROUTINE rademacher_dist

    SUBROUTINE uniform_dist(sigma, R)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: sigma
        REAL(realk), INTENT(out) :: R

        !CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(R)
        R = 2 * SQRT(3.0) * sigma * (R - 0.5)

    END SUBROUTINE uniform_dist

    ! TODO: implement polar method for gaussian distribution (https://de.wikipedia.org/wiki/Polar-Methode)

    ! from: Simulation of truncated normal variables, Christian Robert, Statistics and Computing (1995) 5, 121-125
    ! TODO: fully understand paper and potentially optimize this
    SUBROUTINE gaussian_dist2(mu, sigma, R)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: mu, sigma
        REAL(realk), INTENT(out) :: R

        ! local variables
        REAL(realk) :: rand1, rand2, P
        LOGICAL :: found

        found = .FALSE.
        !CALL RANDOM_SEED()

        DO WHILE (.NOT. found)

            CALL RANDOM_NUMBER(rand1)
            rand1 = truncation_limit / truncation_factor * (rand1 - 0.5) * 2

            P = EXP(-(rand1 ** 2) / 2)

            CALL RANDOM_NUMBER(rand2)

            IF (rand2 <= P) THEN
                ! linear transformation to match given mean and standard deviation
                R = mu + sigma * truncation_factor * rand1
                found = .TRUE.
            END IF

        END DO

    END SUBROUTINE gaussian_dist2

    SUBROUTINE get_truncation_factor(sym_limit, tcf)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: sym_limit
        REAL(realk), INTENT(out) :: tcf

        ! local variables
        INTEGER(intk) :: i, counter
        REAL(realk) :: alpha, beta, prob_alpha, cprob_alpha, prob_beta, cprob_beta, rhs, diff, eps
        LOGICAL :: found_tcf, abort

        ! given a truncation narrower than APPROXIMATELY [-sym_limit = -1.75, sym_limit = 1.75]
        ! of the (stretched) parent pdf, the target standart deviation of 1 cannot be achieved for the truncated pdf
        IF (sym_limit <= 1.8) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! IMPLICIT
        found_tcf = .FALSE.
        abort = .FALSE.
        eps = 10**(-3)

        ! initialization value
        tcf = 1.0

        ! iteration to numerically determine tcf
        counter = 1
        DO WHILE (.NOT. found_tcf .AND. .NOT. abort)

            alpha = - sym_limit / tcf
            beta = sym_limit / tcf

            DO i = 1, SIZE(sn_x)
                IF (alpha <= sn_x(i)) THEN
                    prob_alpha = sn_p(i -1) + (sn_p(i) - sn_p(i-1)) * (alpha - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
                    cprob_alpha = sn_pc(i -1) + (sn_pc(i) - sn_pc(i-1)) * (alpha - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
                    EXIT
                END IF
            END DO

            DO i = 1, SIZE(sn_x)
                IF (beta <= sn_x(i)) THEN
                    prob_beta = sn_p(i -1) + (sn_p(i) - sn_p(i-1)) * (beta - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
                    cprob_beta = sn_pc(i -1) + (sn_pc(i) - sn_pc(i-1)) * (beta - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
                    EXIT
                END IF
            END DO

            rhs = 1 / SQRT(1 - sym_limit / tcf * ((prob_beta + prob_alpha) / (cprob_beta - cprob_alpha)))

            IF (ABS(tcf - rhs) / tcf < eps) THEN
                found_tcf = .TRUE.
            END IF

            IF (counter > 10**(6)) THEN
                abort = .TRUE.
            END IF

            IF (rhs < 1.0) THEN
                abort = .TRUE.
            END IF

            IF (counter > 1 .AND. ABS(tcf - rhs) > diff) THEN
                abort = .TRUE.
            END IF

            diff = ABS(tcf - rhs)
            tcf = rhs
            counter = counter + 1

        END DO

        ! EXPLICIT
        !alpha = - sym_limit
        !beta = sym_limit

        !DO i = 1, SIZE(sn_x)
        !    IF (alpha <= sn_x(i)) THEN
        !        prob_alpha = sn_p(i-1) + (sn_p(i) - sn_p(i-1)) * (alpha - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
        !        cprob_alpha = sn_pc(i-1) + (sn_pc(i) - sn_pc(i-1)) * (alpha - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
        !        EXIT
        !    END IF
        !END DO

        !DO i = 1, SIZE(sn_x)
        !    IF (beta <= sn_x(i)) THEN
        !        prob_beta = sn_p(i-1) + (sn_p(i) - sn_p(i-1)) * (beta - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
        !        cprob_beta = sn_pc(i-1) + (sn_pc(i) - sn_pc(i-1)) * (beta - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
        !        EXIT
        !    END IF
        !END DO

        !tcf = 1 / SQRT(1 - sym_limit * ((prob_beta + prob_alpha) / (cprob_beta - cprob_alpha)))

    END SUBROUTINE get_truncation_factor

    SUBROUTINE finish_particle_diffusion()

    END SUBROUTINE finish_particle_diffusion

END MODULE