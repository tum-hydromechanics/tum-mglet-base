MODULE particle_timeintegration_mod

    USE fields_mod
    USE ib_mod
    USE gc_flowstencils_mod

    USE particle_runtimestat_mod
    USE particle_list_mod
    USE particle_interpolation_mod
    USE particle_diffusion_mod
    USE particle_exchange_mod

    IMPLICIT NONE

    TYPE(rk_2n_t) :: prkscheme

CONTAINS

    SUBROUTINE init_particle_timeintegration()

        USE gc_flowstencils_mod
        USE bound_flow_mod

        CALL start_timer(900)
        CALL start_timer(910)

        ! init rk scheme
        CALL prkscheme%init(prkmethod)

        IF (duse_avg_flow) THEN

            BLOCK
                TYPE(field_t), POINTER :: pwu_avg_f, pwv_avg_f, pww_avg_f
                TYPE(field_t), POINTER :: u_avg_f, v_avg_f, w_avg_f
                INTEGER(intk) :: ilevel

                CALL set_field("PWU_AVG", istag=1, buffers=.TRUE.)
                CALL set_field("PWV_AVG", jstag=1, buffers=.TRUE.)
                CALL set_field("PWW_AVG", kstag=1, buffers=.TRUE.)

                CALL get_field(pwu_avg_f, "PWU_AVG")
                CALL get_field(pwv_avg_f, "PWV_AVG")
                CALL get_field(pww_avg_f, "PWW_AVG")

                CALL get_field(u_avg_f, "U_AVG")
                CALL get_field(v_avg_f, "V_AVG")
                CALL get_field(w_avg_f, "W_AVG")

                pwu_avg_f%arr = u_avg_f%arr
                pwv_avg_f%arr = v_avg_f%arr
                pww_avg_f%arr = w_avg_f%arr

                IF (ib%type == "GHOSTCELL") THEN
                    ! mimic setpointvalues (cannot call setpointvalues because AVG fields have no field buffer) 
                    CALL setpointvalues_all('X', pwu_avg_f, pwv_avg_f, pww_avg_f)
                    CALL setpointvalues_all('Y', pwu_avg_f, pwv_avg_f, pww_avg_f)

                    DO ilevel = minlevel, maxlevel
                        CALL connect(ilevel, 1, v1=pwu_avg_f, v2=pwv_avg_f, v3=pww_avg_f, geom=.TRUE.)
                        CALL bound_flow%bound(ilevel, pwu_avg_f, pwv_avg_f, pww_avg_f)
                    END DO

                    CALL setpointvalues_all('Z', pwu_avg_f, pwv_avg_f, pww_avg_f)
                ELSE
                    DO ilevel = minlevel, maxlevel
                        CALL connect(ilevel, 1, v1=pwu_avg_f, v2=pwv_avg_f, v3=pww_avg_f, geom=.TRUE.)
                        CALL bound_flow%bound(ilevel, pwu_avg_f, pwv_avg_f, pww_avg_f)
                    END DO
                END IF
            END BLOCK

        END IF
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
        TYPE(field_t), POINTER :: pwu_f, pwv_f, pww_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: pwu, pwv, pww

        REAL(realk), ALLOCATABLE :: pdx_pot(:), pdy_pot(:), pdz_pot(:)

        INTEGER(intk) :: igrid, i, ii, jj, kk, gfound, ig, temp_grid
        REAL(realk) :: temp_coord(3)
        REAL(realk) :: pu_adv, pv_adv, pw_adv
        REAL(realk) :: pdx_adv, pdy_adv, pdz_adv, pdx_diff, pdy_diff, pdz_diff
        REAL(realk) :: pdx_eff_tot, pdy_eff_tot, pdz_eff_tot, pdx_eff, pdy_eff, pdz_eff

        INTEGER(intk) :: irk
        REAL(realk) :: A, B

        CALL start_timer(900)
        CALL start_timer(920)

        CALL start_timer(921)

        ALLOCATE(pdx_pot(my_particle_list%ifinal))
        ALLOCATE(pdy_pot(my_particle_list%ifinal))
        ALLOCATE(pdz_pot(my_particle_list%ifinal))

        pdx_pot = 0.0
        pdy_pot = 0.0
        pdz_pot = 0.0

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        IF (dinterp_padvection) THEN
            CALL get_field(dx_f, "DX")
            CALL get_field(dy_f, "DY")
            CALL get_field(dz_f, "DZ")
            CALL get_field(ddx_f, "DDX")
            CALL get_field(ddy_f, "DDY")
            CALL get_field(ddz_f, "DDZ")
        END IF

        IF (duse_avg_flow) THEN
            ! use the point values deduced from the average flow field
            CALL get_field(pwu_f, "PWU_AVG")
            CALL get_field(pwv_f, "PWV_AVG")
            CALL get_field(pww_f, "PWW_AVG")
        ELSE
            IF (ib%type == "GHOSTCELL") THEN
                CALL get_field(pwu_f, "PWU")
                CALL get_field(pwv_f, "PWV")
                CALL get_field(pww_f, "PWW")
            ELSE
                CALL get_field(pwu_f, "U")
                CALL get_field(pwv_f, "V")
                CALL get_field(pww_f, "W")
            END IF
        END IF

        CALL stop_timer(921)

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
            temp_coord(1) = my_particle_list%particles(i)%x
            temp_coord(2) = my_particle_list%particles(i)%y
            temp_coord(3) = my_particle_list%particles(i)%z

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

            ! Grid and Field Info
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL x_f%get_ptr(x, igrid)
            CALL y_f%get_ptr(y, igrid)
            CALL z_f%get_ptr(z, igrid)

            IF (dinterp_padvection) THEN
                CALL dx_f%get_ptr(dx, igrid)
                CALL dy_f%get_ptr(dy, igrid)
                CALL dz_f%get_ptr(dz, igrid)
                CALL ddx_f%get_ptr(ddx, igrid)
                CALL ddy_f%get_ptr(ddy, igrid)
                CALL ddz_f%get_ptr(ddz, igrid)
            END IF

            CALL pwu_f%get_ptr(pwu, igrid)
            CALL pwv_f%get_ptr(pwv, igrid)
            CALL pww_f%get_ptr(pww, igrid)

            CALL stop_timer(921)

            ! for particle runtime statistics (terminal output)
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                pdx_eff_tot = 0.0
                pdy_eff_tot = 0.0
                pdz_eff_tot = 0.0
            END IF


            DO irk = 1, prkscheme%nrk

                CALL start_timer(921)
                CALL prkscheme%get_coeffs(A, B, irk)

                ! should be obsolete as the effective displacement is also zeroized in move_particle
                pdx_eff = 0.0
                pdy_eff = 0.0
                pdz_eff = 0.0

                ! get particle velocity
                IF (dinterp_padvection) THEN
                    CALL interpolate_lincon(my_particle_list%particles(i), kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                     pwu, pwv, pww, pu_adv, pv_adv, pw_adv)
                ELSE
                    CALL get_nearest_value(my_particle_list%particles(i), kk, jj, ii, x, y, z, &
                     pwu, pwv, pww, pu_adv, pv_adv, pw_adv)
                END IF

                CALL prkstep(pdx_pot(i), pdy_pot(i), pdz_pot(i), pu_adv, pv_adv, pw_adv, dt, A, B, pdx_adv, pdy_adv, pdz_adv)

                ! for debugging
                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("---------- Advection RK Step: ", I0, " ----------")') irk
                    WRITE(*,'("Intermediate Velocity ", F12.9, " ", F12.9, " ", F12.9)') pu_adv, pv_adv, pw_adv
                    WRITE(*,'("Intermediate (potential) Displacement ", F12.9, " ", F12.9, " ", F12.9)') pdx_adv, pdy_adv, pdz_adv
                    WRITE(*, '()')
                END IF
                CALL stop_timer(921)

                CALL start_timer(922)
                ! Particle Boundary Interaction
                CALL move_particle(my_particle_list%particles(i), pdx_adv, pdy_adv, pdz_adv, &
                 pdx_eff, pdy_eff, pdz_eff, temp_coord, temp_grid)
                CALL stop_timer(922)

                CALL start_timer(921)
                pdx_pot(i) = pdx_eff / B
                pdy_pot(i) = pdy_eff / B
                pdz_pot(i) = pdz_eff / B

                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    ! for particle runtime statistics (terminal output)
                    pdx_eff_tot = pdx_eff_tot + pdx_eff
                    pdy_eff_tot = pdy_eff_tot + pdy_eff
                    pdz_eff_tot = pdz_eff_tot + pdz_eff
                    psim_max_adv_dx = MAX(psim_max_adv_dx, ABS(pdx_eff_tot))
                    psim_max_adv_dy = MAX(psim_max_adv_dy, ABS(pdy_eff_tot))
                    psim_max_adv_dz = MAX(psim_max_adv_dz, ABS(pdz_eff_tot))
                END IF

                CALL stop_timer(921)

            END DO

            ! --- DIFFSUION ---

            IF (ddiffusion) THEN

                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("---------- Particle Diffusion ----------")')
                    WRITE(*, '()')
                END IF

                CALL start_timer(924)
                CALL generate_diffusive_displacement(dt, D(1), D(2), D(3), pdx_diff, pdy_diff, pdz_diff)
                CALL stop_timer(924)

                CALL start_timer(925)
                CALL move_particle(my_particle_list%particles(i), pdx_diff, pdy_diff, pdz_diff, &
                 pdx_eff, pdy_eff, pdz_eff, temp_coord, temp_grid)
                CALL stop_timer(925)

                ! for particle runtime statistics (terminal output)
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    pdx_eff_tot = pdx_eff_tot + pdx_eff
                    pdy_eff_tot = pdy_eff_tot + pdy_eff
                    pdz_eff_tot = pdz_eff_tot + pdz_eff
                    psim_max_dif_dx = MAX(psim_max_dif_dx, ABS(pdx_eff))
                    psim_max_dif_dy = MAX(psim_max_dif_dy, ABS(pdy_eff))
                    psim_max_dif_dz = MAX(psim_max_dif_dz, ABS(pdz_eff))
                END IF

            END IF

            ! for particle runtime statistics (terminal output)
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                IF (psim_max_disp < SQRT(pdx_eff_tot**2 + pdy_eff_tot**2 + pdz_eff_tot**2)) THEN
                    psim_max_dx = pdx_eff_tot
                    psim_max_dy = pdy_eff_tot
                    psim_max_dz = pdz_eff_tot
                    psim_max_disp = SQRT(pdx_eff_tot**2 + pdy_eff_tot**2 + pdz_eff_tot**2)
                END IF
            END IF

        END DO

        DEALLOCATE(pdx_pot)
        DEALLOCATE(pdy_pot)
        DEALLOCATE(pdz_pot)

        CALL stop_timer(920)
        CALL stop_timer(900)

    END SUBROUTINE timeintegrate_particles

    SUBROUTINE finish_particle_timeintegration()

    END SUBROUTINE finish_particle_timeintegration

END MODULE
