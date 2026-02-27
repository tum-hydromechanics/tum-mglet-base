MODULE particle_diffusion_mod

    USE, INTRINSIC :: ISO_FORTRAN_ENV
    USE, INTRINSIC :: ISO_C_BINDING

    USE omp_lib

    USE MPI_f08
    USE precision_mod
    USE charfunc_mod
    USE comms_mod
    USE fort7_mod
    USE grids_mod
    USE field_mod
    USE fields_mod
    USE connect2_mod

    USE particle_config_mod
    USE particle_rng_mod
    USE particle_basetype_mod
    USE particle_interpolation_mod

    IMPLICIT NONE

    ! REFERENCE VALUES FOR THE STANDART NORMAL DISTRIBUTIONS:
    ! realization values
    REAL(realk) :: sn_x(22) = [-3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, &
       1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,  2.8,  3.0]
    ! corresponding cumulative probabilities
    REAL(realk) :: sn_pc(22) = [0.00135, 0.00256, 0.00466, 0.00820, 0.01390, 0.02275, 0.03593, 0.05480, 0.08076, 0.11507, 0.15866, &
      0.84134, 0.88493, 0.91924, 0.94520, 0.96407, 0.97725, 0.98610, 0.99180, 0.99534, 0.99744, 0.99865]
    ! corresponding probabilites
    REAL(realk) :: sn_p(22)

    ! truncation limit stored in config mod
    REAL(realk) :: truncation_factor

    !$omp declare target link(truncation_factor)

CONTAINS

    SUBROUTINE init_particle_diffusion()

        ! local_variables
        INTEGER(intk) :: i

        CALL start_timer(900)
        CALL start_timer(910)

        DO i = 1, SIZE(sn_p)
            sn_p(i) = 1.0_realk / SQRT(2.0_realk * pi) * EXP(- (sn_x(i)**2) / 2.0_realk)
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

#ifdef _MGLET_OPENMP_
        CALL init_parallel_lcg()
#endif

        !$omp target enter data map(to: truncation_factor)

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_diffusion


    SUBROUTINE generate_diffusive_displacement(dt, D_x, D_y, D_z, pdx, pdy, pdz)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: D_x, D_y, D_z
        REAL(realk), INTENT(out) :: pdx, pdy, pdz

        ! local variables
        REAL(realk) :: sigx, sigy, sigz, ranx, rany, ranz

        IF (D_x > 0.0_realk) THEN
            sigx = SQRT(2 * D_x * dt)

            SELECT CASE (lower(TRIM(random_walk_mode)))
            CASE ("rademacher")
                CALL rademacher_dist(sigx, ranx)
            CASE ("uniform")
                CALL uniform_dist(sigx, ranx)
            CASE ("gaussian2")
                CALL gaussian_dist2(0.0_realk, sigx, truncation_limit, truncation_factor, ranx)
            END SELECT

            pdx = ranx ! diffusion length

        END IF

        IF (D_y > 0.0_realk) THEN

            sigy = SQRT(2 * D_y * dt)

            SELECT CASE (lower(TRIM(random_walk_mode)))
            CASE ("rademacher")
                CALL rademacher_dist(sigy, rany)
            CASE ("uniform")
                CALL uniform_dist(sigy, rany)
            CASE ("gaussian2")
                CALL gaussian_dist2(0.0_realk, sigy, truncation_limit, truncation_factor, rany)
            END SELECT

            pdy = rany ! diffusion length

        END IF

        IF (D_z > 0.0_realk) THEN

            sigz = SQRT(2 * D_z * dt)

            SELECT CASE (lower(TRIM(random_walk_mode)))
            CASE ("rademacher")
                CALL rademacher_dist(sigz, ranz)
            CASE ("uniform")
                CALL uniform_dist(sigz, ranz)
            CASE ("gaussian2")
                CALL gaussian_dist2(0.0_realk, sigz, truncation_limit, truncation_factor, ranz)
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
    ! TODO: potentially optimize this
    SUBROUTINE gaussian_dist2(mu, sigma, trunc_limit, trunc_factor, R)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: mu, sigma, trunc_limit, trunc_factor
        REAL(realk), INTENT(out) :: R

        ! local variables
        REAL(realk) :: rand1, rand2, P
        LOGICAL :: found

        found = .FALSE.

        DO WHILE (.NOT. found)

            CALL RANDOM_NUMBER(rand1)
            rand1 = trunc_limit / trunc_factor * (rand1 - 0.5) * 2.0

            P = EXP(-(rand1 ** 2) / 2)

            CALL RANDOM_NUMBER(rand2)

            IF (rand2 <= P) THEN
                ! linear transformation to match given mean and standard deviation
                R = mu + sigma * trunc_factor * rand1
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
        eps = 0.001

        ! initialization value
        tcf = 1.0

        ! iteration to numerically determine tcf
        counter = 1
        DO WHILE (.NOT. found_tcf .AND. .NOT. abort)

            alpha = - sym_limit / tcf
            beta = sym_limit / tcf

            IF (alpha < sn_x(1) .OR. alpha > sn_x(SIZE(sn_x)) .OR. beta < sn_x(1) .OR. beta > sn_x(SIZE(sn_x))) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            DO i = 2, SIZE(sn_x)
                IF (alpha <= sn_x(i)) THEN
                    prob_alpha = sn_p(i-1) + (sn_p(i) - sn_p(i-1)) * (alpha - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
                    cprob_alpha = sn_pc(i-1) + (sn_pc(i) - sn_pc(i-1)) * (alpha - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
                    EXIT
                END IF
            END DO

            DO i = 2, SIZE(sn_x)
                IF (beta <= sn_x(i)) THEN
                    prob_beta = sn_p(i-1) + (sn_p(i) - sn_p(i-1)) * (beta - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
                    cprob_beta = sn_pc(i-1) + (sn_pc(i) - sn_pc(i-1)) * (beta - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
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

        !IF (alpha < sn_x(1) .OR. alpha > sn_x(SIZE(sn_x)) .OR. beta < sn_x(1) .OR. beta > sn_x(SIZE(sn_x))) THEN
        !    CALL errr(__FILE__, __LINE__)
        !END IF

        !DO i = 2, SIZE(sn_x)
        !    IF (alpha <= sn_x(i)) THEN
        !        prob_alpha = sn_p(i-1) + (sn_p(i) - sn_p(i-1)) * (alpha - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
        !        cprob_alpha = sn_pc(i-1) + (sn_pc(i) - sn_pc(i-1)) * (alpha - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
        !        EXIT
        !    END IF
        !END DO

        !DO i = 2, SIZE(sn_x)
        !    IF (beta <= sn_x(i)) THEN
        !        prob_beta = sn_p(i-1) + (sn_p(i) - sn_p(i-1)) * (beta - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
        !        cprob_beta = sn_pc(i-1) + (sn_pc(i) - sn_pc(i-1)) * (beta - sn_x(i-1)) / (sn_x(i) - sn_x(i-1))
        !        EXIT
        !    END IF
        !END DO

        !tcf = 1 / SQRT(1 - sym_limit * ((prob_beta + prob_alpha) / (cprob_beta - cprob_alpha)))

    END SUBROUTINE get_truncation_factor

    SUBROUTINE generate_diffusive_displacement_target(dt, D_x, D_y, D_z, pdx, pdy, pdz, seed)

        !$omp declare target

        ! subroutine arguments
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: D_x, D_y, D_z
        REAL(realk), INTENT(out) :: pdx, pdy, pdz
        INTEGER(c_int), INTENT(inout) :: seed

        ! local variables
        REAL(realk) :: sigx, sigy, sigz, ranx, rany, ranz

        sigx = SQRT(2 * D_x * dt)
        !CALL gaussian_dist_target(0.0_realk, sigx, truncation_limit, truncation_factor, seed, ranx)
        CALL uniform_dist_target(sigx, ranx, seed)
        pdx = ranx ! diffusion length

        sigy = SQRT(2 * D_y * dt)
        !CALL gaussian_dist_target(0.0_realk, sigy, truncation_limit, truncation_factor, seed, rany)
        CALL uniform_dist_target(sigy, rany, seed)
        pdy = rany ! diffusion length

        sigz = SQRT(2 * D_z * dt)
        !CALL gaussian_dist_target(0.0_realk, sigz, truncation_limit, truncation_factor, seed, ranz)
        CALL uniform_dist_target(sigz, ranz, seed)
        pdz = ranz ! diffusion length

    END SUBROUTINE generate_diffusive_displacement_target

    SUBROUTINE gaussian_dist_target(mu, sigma, trunc_limit, trunc_factor, seed, R)

        !$omp declare target

        ! subroutine arguments
        REAL(realk), INTENT(in) :: mu, sigma, trunc_limit, trunc_factor
        INTEGER(c_int), INTENT(inout) :: seed
        REAL(realk), INTENT(out) :: R
        
        ! local variables
        REAL(realk) :: rand1, rand2, P
        LOGICAL :: found

        found = .FALSE.

        DO WHILE (.NOT. found)

#if defined __GFORTRAN__
            CALL RANDOM_NUMBER(rand1)
#else
            CALL lcg(seed, rand1)
#endif

            rand1 = trunc_limit / trunc_factor * (rand1 - 0.5) * 2.0

            P = EXP(-(rand1 ** 2) / 2)

#if defined __GFORTRAN__
            CALL RANDOM_NUMBER(rand2)
#else
            CALL lcg(seed, rand2)
#endif

            IF (rand2 <= P) THEN
                ! linear transformation to match given mean and standard deviation
                R = mu + sigma * trunc_factor * rand1
                found = .TRUE.
            END IF

        END DO

    END SUBROUTINE gaussian_dist_target

    SUBROUTINE uniform_dist_target(sigma, R, seed)
        
        !$omp declare target
        
        ! subroutine arguments
        
        REAL(realk), INTENT(in) :: sigma
        REAL(realk), INTENT(out) :: R
        INTEGER(c_int), INTENT(inout) :: seed

#if defined __GFORTRAN__
        CALL RANDOM_NUMBER(R)
#else
        CALL lcg(seed, R)
#endif
        
        R = 2 * SQRT(3.0) * sigma * (R - 0.5)

    END SUBROUTINE uniform_dist_target


    SUBROUTINE init_custom_prng()

        

    END SUBROUTINE init_custom_prng


    SUBROUTINE finish_particle_diffusion()

        !$omp target exit data map(delete: truncation_factor)

    END SUBROUTINE finish_particle_diffusion

END MODULE