MODULE particle_timeintegration_mod

    USE fields_mod
    USE particle_interpolation_mod
    USE particle_exchange_mod
    USE particle_list_mod
    USE particle_boundaries_mod

    IMPLICIT NONE

    TYPE(rk_2n_t) :: prkscheme

CONTAINS

    SUBROUTINE init_particle_timeintegration()

        CALL start_timer(900)
        CALL start_timer(910)

        ! init rk scheme
        CALL prkscheme%init(prkmethod)

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_timeintegration

    SUBROUTINE timeintegrate_particles(itstep, dt)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt

        ! local variables
        !TYPE(baseparticle_t) :: particle_clone
        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: u_f, v_f, w_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w

        REAL(realk), ALLOCATABLE :: pdx(:), pdy(:), pdz(:)

        INTEGER(intk) :: igrid, i, j, ii, jj, kk, gfound, ig, destgrid
        REAL(realk) :: pu_adv, pv_adv, pw_adv, pu_diff, pv_diff, pw_diff
        REAL(realk) :: pdx_diff, pdy_diff, pdz_diff, pdx_eff, pdy_eff, pdz_eff

        INTEGER(intk) :: irk
        REAL(realk) :: A, B

        CALL start_timer(900)
        CALL start_timer(920)

        ALLOCATE(pdx(my_particle_list%ifinal))
        ALLOCATE(pdy(my_particle_list%ifinal))
        ALLOCATE(pdz(my_particle_list%ifinal))

        pdx = 0.0
        pdy = 0.0
        pdz = 0.0

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        IF (dinterp_particles) THEN

            CALL get_field(dx_f, "DX")
            CALL get_field(dy_f, "DY")
            CALL get_field(dz_f, "DZ")

            CALL get_field(ddx_f, "DDX")
            CALL get_field(ddy_f, "DDY")
            CALL get_field(ddz_f, "DDZ")

        END IF

        ! REMOVE later
        IF (myid == 0) THEN
            WRITE(*, *) ''
            WRITE(*, '("NEW TIMESTEP - PARTICLE TIMEINTEGRATION:")')
        END IF

        ! algorithm for EXPLICIT RK schemes
        DO i = 1, my_particle_list%ifinal

            ! checking activity
            IF (my_particle_list%particles(i)%state < 1) THEN
                CYCLE
            END IF

            ! checking locality (Debug)
            IF (my_particle_list%particles(i)%iproc /= myid) THEN
                WRITE(*, '("ERROR: Particle on wrong proc at start of current timestep")')
                CALL errr(__FILE__, __LINE__)
            END IF

            ! assigning igrid for clarity of the follwing expressions
            igrid = my_particle_list%particles(i)%igrid

            ! checking consistency (Debug)
            gfound = 1
            DO ig = 1, nMyGrids
                IF (my_particle_list%particles(i)%igrid == mygrids(ig)) gfound = 1; EXIT
            END DO

            IF (gfound == 0) THEN
                WRITE(*, '("ERROR: Particle grid ", I0, " is not on this process (", I0, ")")') my_particle_list%particles(i)%igrid, myid
                CALL errr(__FILE__, __LINE__)
            END IF

            ! for debugging
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("Pre Motion - Particle Status:")')
                    CALL print_particle_status(my_particle_list%particles(i))
                    WRITE(*, '()')
            END SELECT

            DO irk = 1, prkscheme%nrk

                CALL prkscheme%get_coeffs(A, B, irk)

                ! for debugging
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*,'("RK Step: ", I0)') irk
                        WRITE(*,'("LS RK Coefficients: ", 2F12.8)') A, B
                        WRITE(*, '()')
                END SELECT

                ! --- ADVECTION VELOCITY ---

                CALL get_mgdims(kk, jj, ii, igrid)

                ! Grid and Field Info
                CALL x_f%get_ptr(x, igrid)
                CALL y_f%get_ptr(y, igrid)
                CALL z_f%get_ptr(z, igrid)

                CALL u_f%get_ptr(u, igrid)
                CALL v_f%get_ptr(v, igrid)
                CALL w_f%get_ptr(w, igrid)

                IF (dinterp_particles) THEN

                    CALL dx_f%get_ptr(dx, igrid)
                    CALL dy_f%get_ptr(dy, igrid)
                    CALL dz_f%get_ptr(dz, igrid)

                    CALL ddx_f%get_ptr(ddx, igrid)
                    CALL ddy_f%get_ptr(ddy, igrid)
                    CALL ddz_f%get_ptr(ddz, igrid)

                    CALL gobert_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                    pu_adv, pv_adv, pw_adv, u, v, w, x, y, z, dx, dy, dz, ddx, ddy, ddz)

                ELSE

                    CALL nearest_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                    pu_adv, pv_adv, pw_adv, u, v, w, x, y, z)

                END IF

                CALL prkstep(pdx(i), pdy(i), pdz(i), pu_adv, pv_adv, pw_adv, dt, A, B)

                ! Particle Boundary Interaction
                CALL move_particle(my_particle_list%particles(i), pdx(i), pdy(i), pdz(i), pdx_eff, pdy_eff, pdz_eff)

                ! TODO: maybe source this out into move_particle (?)
                IF (irk /= prkscheme%nrk) THEN
                    pdx(i) = pdx_eff
                    pdy(i) = pdy_eff
                    pdz(i) = pdz_eff
                END IF

                IF (irk == prkscheme%nrk) THEN

                    ! --- Diffusiion ---
                    CALL get_particle_diffusion(dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff)

                    CALL move_particle(my_particle_list%particles(i), pdx_diff, pdy_diff, pdz_diff, pdx_eff, pdy_eff, pdz_eff)

                    ! for debugging
                    SELECT CASE (TRIM(particle_terminal))
                        CASE ("none")
                            CONTINUE
                        CASE ("normal")
                            CONTINUE
                        CASE ("verbose")
                            WRITE(*,'("Partical Diffusion: ")')
                            WRITE(*, '()')
                    END SELECT

                END IF

            END DO

        END DO

        CALL stop_timer(920)
        CALL stop_timer(900)

    END SUBROUTINE timeintegrate_particles

    SUBROUTINE get_particle_diffusion(dt, pu_diff, pv_diff, pw_diff, dx_diff, dy_diff, dz_diff)

            ! subroutine arguments
            REAL(realk), INTENT(in) :: dt
            REAL(realk), INTENT(out) :: pu_diff, pv_diff, pw_diff, dx_diff, dy_diff, dz_diff

            ! local variables
            REAL(realk) :: rand(3)

            CALL start_timer(922)

            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(rand)

            rand(1) = rand(1) - 0.5_realk
            rand(2) = rand(2) - 0.5_realk
            rand(3) = rand(3) - 0.5_realk

            pu_diff = SQRT(2 * D(1) / dt) * rand(1) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))
            pv_diff = SQRT(2 * D(2) / dt) * rand(2) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))
            pw_diff = SQRT(2 * D(3) / dt) * rand(3) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))

            dx_diff = SQRT(2 * D(1) * dt) * rand(1) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))
            dy_diff = SQRT(2 * D(2) * dt) * rand(2) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))
            dz_diff = SQRT(2 * D(3) * dt) * rand(3) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))

            CALL stop_timer(922)

    END SUBROUTINE get_particle_diffusion

    SUBROUTINE finish_particle_timeintegration()

    END SUBROUTINE finish_particle_timeintegration

END MODULE
