MODULE particle_timeintegration_mod

    USE fields_mod

    USE particle_list_mod
    USE particle_interpolation_mod
    USE particle_diffusion_mod
    USE particle_exchange_mod

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
        TYPE(field_t), POINTER :: diffx_f, diffy_f, diffz_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: diffx, diffy, diffz

        REAL(realk), ALLOCATABLE :: pdx_pot(:), pdy_pot(:), pdz_pot(:)

        INTEGER(intk) :: igrid, i, j, ii, jj, kk, gfound, ig, temp_grid
        REAL(realk) :: pu_adv, pv_adv, pw_adv, p_diffx, p_diffy, p_diffz
        REAL(realk) :: pdx_adv, pdy_adv, pdz_adv, pdx_diff, pdy_diff, pdz_diff, pdx_eff, pdy_eff, pdz_eff

        INTEGER(intk) :: irk
        REAL(realk) :: A, B

        CALL start_timer(900)
        CALL start_timer(920)

        ALLOCATE(pdx_pot(my_particle_list%ifinal))
        ALLOCATE(pdy_pot(my_particle_list%ifinal))
        ALLOCATE(pdz_pot(my_particle_list%ifinal))

        pdx_pot = 0.0
        pdy_pot = 0.0
        pdz_pot = 0.0

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        IF (dinterp_padvection .OR. dinterp_pdiffsuion) THEN
            CALL get_field(dx_f, "DX")
            CALL get_field(dy_f, "DY")
            CALL get_field(dz_f, "DZ")
            CALL get_field(ddx_f, "DDX")
            CALL get_field(ddy_f, "DDY")
            CALL get_field(ddz_f, "DDZ")
        END IF

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        IF (dturb_diff) THEN
            CALL get_field(diffx_f, "P_DIFF_X")
            CALL get_field(diffy_f, "P_DIFF_Y")
            CALL get_field(diffz_f, "P_DIFF_Z")
        END IF

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, *) ''
                WRITE(*, '("=== TIMESTEP ", I0, " - PARTICLE TIMEINTEGRATION ===")') itstep
                WRITE(*, *) ''
            END IF
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
            temp_grid = my_particle_list%particles(i)%igrid

            CALL get_mgdims(kk, jj, ii, igrid)

            ! Grid and Field Info
            CALL x_f%get_ptr(x, igrid)
            CALL y_f%get_ptr(y, igrid)
            CALL z_f%get_ptr(z, igrid)

            IF (dinterp_padvection .OR. dinterp_pdiffsuion) THEN
                CALL dx_f%get_ptr(dx, igrid)
                CALL dy_f%get_ptr(dy, igrid)
                CALL dz_f%get_ptr(dz, igrid)
                CALL ddx_f%get_ptr(ddx, igrid)
                CALL ddy_f%get_ptr(ddy, igrid)
                CALL ddz_f%get_ptr(ddz, igrid)
            END IF

            ! checking consistency (Debug)
            gfound = 1
            DO ig = 1, nMyGrids
                IF (igrid == mygrids(ig)) gfound = 1; EXIT
            END DO

            IF (gfound == 0) THEN
                WRITE(*, '("ERROR: Particle grid ", I0, " is not on this process (", I0, ")")') igrid, myid
                CALL errr(__FILE__, __LINE__)
            END IF

            ! for debugging
            IF (TRIM(particle_terminal) == "verbose") THEN
                WRITE(*,'("Pre Motion - Particle Status:")')
                CALL print_particle_status(my_particle_list%particles(i))
                WRITE(*, '()')
            END IF

            ! --- ADVECTION ---
            CALL start_timer(921)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            DO irk = 1, prkscheme%nrk


                CALL prkscheme%get_coeffs(A, B, irk)

                ! should be obsolete as the effective displacement is also zeroized in move_particle
                pdx_eff = 0.0
                pdy_eff = 0.0
                pdz_eff = 0.0

                ! get particle velocity
                IF (dinterp_padvection) THEN
                    CALL interpolate_lincon(my_particle_list%particles(i), kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                     u, v, w, pu_adv, pv_adv, pw_adv)
                ELSE
                    CALL get_nearest_value(my_particle_list%particles(i), kk, jj, ii, x, y, z, &
                     u, v, w, pu_adv, pv_adv, pw_adv)
                END IF

                CALL prkstep(pdx_pot(i), pdy_pot(i), pdz_pot(i), pu_adv, pv_adv, pw_adv, dt, A, B, pdx_adv, pdy_adv, pdz_adv)

                ! for debugging
                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("---------- Advection RK Step: ", I0, " ----------")') irk
                    WRITE(*,'("Intermediate Velocity ", F12.9, " ", F12.9, " ", F12.9)') pu_adv, pv_adv, pw_adv
                    WRITE(*,'("Intermediate Displacement ", F12.9, " ", F12.9, " ", F12.9)') pdx_adv, pdy_adv, pdz_adv
                    WRITE(*, '()')
                END IF
                CALL stop_timer(921)

                CALL start_timer(922)
                ! Particle Boundary Interaction
                CALL move_particle(my_particle_list%particles(i), pdx_adv, pdy_adv, pdz_adv, pdx_eff, pdy_eff, pdz_eff, temp_grid)
                CALL stop_timer(922)

                CALL start_timer(921)
                pdx_pot(i) = pdx_eff / B
                pdy_pot(i) = pdy_eff / B
                pdz_pot(i) = pdz_eff / B

            END DO

            CALL stop_timer(921)

            ! --- DIFFSUION ---

            IF (ddiffusion) THEN
                CALL start_timer(923)

                IF (dturb_diff) THEN
                    CALL diffx_f%get_ptr(diffx, my_particle_list%particles(i)%igrid)
                    CALL diffy_f%get_ptr(diffy, my_particle_list%particles(i)%igrid)
                    CALL diffz_f%get_ptr(diffz, my_particle_list%particles(i)%igrid)
                END IF

                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("---------- Particle Diffusion ----------")')
                    WRITE(*, '()')
                END IF

                IF (dturb_diff) THEN
                    IF (dinterp_pdiffsuion) THEN
                        CALL interpolate_lincon(my_particle_list%particles(i), kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                        diffx, diffy, diffz, p_diffx, p_diffy, p_diffz)
                    ELSE
                        CALL get_nearest_value(my_particle_list%particles(i), kk, jj, ii, x, y, z, &
                        diffx, diffy, diffz, p_diffx, p_diffy, p_diffz)
                    END IF
                ELSE
                    p_diffx = D(1)
                    p_diffy = D(2)
                    p_diffz = D(3)
                END IF

                CALL stop_timer(923)

                CALL start_timer(924)
                CALL generate_diffusive_displacement(dt, p_diffx, p_diffy, p_diffz, pdx_diff, pdy_diff, pdz_diff)
                CALL stop_timer(924)

                CALL start_timer(925)
                CALL move_particle(my_particle_list%particles(i), pdx_diff, pdy_diff, pdz_diff, pdx_eff, pdy_eff, pdz_eff, temp_grid)
                CALL stop_timer(925)
            END IF
        END DO

        ! TODO: print out particle statistics (terminal)

        DEALLOCATE(pdx_pot)
        DEALLOCATE(pdy_pot)
        DEALLOCATE(pdz_pot)

        CALL stop_timer(920)
        CALL stop_timer(900)

    END SUBROUTINE timeintegrate_particles

    SUBROUTINE finish_particle_timeintegration()

    END SUBROUTINE finish_particle_timeintegration

END MODULE
