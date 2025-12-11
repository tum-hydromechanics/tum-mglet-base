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

        CALL start_timer(900)
        CALL start_timer(910)

        ! init rk scheme
        CALL prkscheme%init(prkmethod)

        IF (duse_avg_flow) THEN

            BLOCK
                TYPE(field_t), POINTER :: pwu_avg_f, pwv_avg_f, pww_avg_f
                TYPE(field_t), POINTER :: u_avg_f, v_avg_f, w_avg_f

                CALL set_field("PWU_AVG", istag=1, buffers=.TRUE.)
                CALL set_field("PWV_AVG", jstag=1, buffers=.TRUE.)
                CALL set_field("PWW_AVG", kstag=1, buffers=.TRUE.)

                CALL get_field(pwu_avg_f, "PWU_AVG")
                CALL get_field(pwv_avg_f, "PWV_AVG")
                CALL get_field(pww_avg_f, "PWW_AVG")

                CALL get_field(u_avg_f, "U_AVG")
                CALL get_field(v_avg_f, "V_AVG")
                CALL get_field(w_avg_f, "W_AVG")

                CALL setpointvalues(pwu_avg_f, pwv_avg_f, pww_avg_f, u_avg_f, v_avg_f, w_avg_f, .TRUE.)
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
        INTEGER(intk) :: igrid, i, j, ii, jj, kk, gfound, ig, ipart, temp_grid
        REAL(realk) :: pd_eff_tot(3), temp_coord(3)
        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: pwu_f, pwv_f, pww_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: pwu, pwv, pww

        CALL start_timer(900)
        CALL start_timer(920)

        CALL start_timer(921)

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

        CALL count_pog(my_particle_list, grids_np, plist_displ)
        
        DO i = 1, nmy_particle_grids
            igrid = my_particle_grids(i)

            ! checking consistency (Debug)
            gfound = 1
            DO ig = 1, nmygrids
                IF (igrid == mygrids(ig)) gfound = 1; EXIT
            END DO
            
            IF (gfound == 0) THEN
                WRITE(*, '("ERROR: Particle grid ", I0, " is not on this process (", I0, ")")') igrid, myid
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL start_timer(921)

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

            DO j = 1, grids_np(i)

                ipart = plist_displ(i) + j

                ! checking activity
                IF (my_particle_list%particles(ipart)%state < 1) THEN
                    CYCLE
                END IF

                ! checking locality (Debug)
                IF (my_particle_list%particles(ipart)%iproc /= myid) THEN
                    WRITE(*, '("ERROR: Particle on wrong proc at start of current timestep")')
                    CALL errr(__FILE__, __LINE__)
                END IF
                
                IF (.NOT. (igrid == my_particle_list%particles(ipart)%igrid)) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF

                ! for debugging
                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("Pre Motion - Particle Status:")')
                    CALL print_particle_status(my_particle_list%particles(i))
                    WRITE(*, '()')
                END IF

                ! for particle runtime statistics (terminal output)
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    pd_eff_tot = 0.0
                END IF

                temp_grid = my_particle_list%particles(ipart)%igrid
                temp_coord(1) = my_particle_list%particles(ipart)%x
                temp_coord(2) = my_particle_list%particles(ipart)%y
                temp_coord(3) = my_particle_list%particles(ipart)%z

                CALL particle_advection(my_particle_list%particles(ipart), temp_grid, temp_coord, pd_eff_tot, &
                 kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, pwu, pwv, pww, dt)
                
                CALL particle_diffusion(my_particle_list%particles(ipart), temp_grid, temp_coord, pd_eff_tot, dt)

                ! for particle runtime statistics (terminal output)
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    IF (psim_max_disp < SQRT(pd_eff_tot(1)**2 + pd_eff_tot(2)**2 + pd_eff_tot(3)**2)) THEN
                        psim_max_dx = pd_eff_tot(1)
                        psim_max_dy = pd_eff_tot(2)
                        psim_max_dz = pd_eff_tot(3)
                        psim_max_disp = SQRT(pd_eff_tot(1)**2 + pd_eff_tot(2)**2 + pd_eff_tot(3)**2)
                    END IF
                END IF
            END DO 
        END DO

        CALL stop_timer(920)
        CALL stop_timer(900)

    END SUBROUTINE timeintegrate_particles

    SUBROUTINE particle_advection(particle, temp_grid, temp_coord, pd_eff_tot, kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, pwu, pwv, pww, dt)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(inout) :: temp_grid
        REAL(realk), INTENT(inout) :: temp_coord(3)
        REAL(realk), INTENT(inout) :: pd_eff_tot(3)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :), INTENT(in) :: pwu, pwv, pww
        REAL(realk), INTENT(in) :: dt
        
        ! local variables
        INTEGER(intk) :: irk
        REAL(realk) :: A, B
        REAL(realk) :: pu_adv, pv_adv, pw_adv
        REAL(realk) :: pdx_adv, pdy_adv, pdz_adv
        REAL(realk) :: pdx_pot, pdy_pot, pdz_pot
        REAL(realk) :: pdx_eff, pdy_eff, pdz_eff

        CALL start_timer(921)
        DO irk = 1, prkscheme%nrk

            CALL prkscheme%get_coeffs(A, B, irk)

            ! should be obsolete as the effective displacement is also zeroized in move_particle
            pdx_eff = 0.0
            pdy_eff = 0.0
            pdz_eff = 0.0

            ! get particle velocity
            IF (dinterp_padvection) THEN
                CALL interpolate_lincon(particle, kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                 pwu, pwv, pww, pu_adv, pv_adv, pw_adv)
            ELSE
                CALL get_nearest_value(particle, kk, jj, ii, x, y, z, &
                 pwu, pwv, pww, pu_adv, pv_adv, pw_adv)
            END IF

            CALL prkstep(pdx_pot, pdy_pot, pdz_pot, pu_adv, pv_adv, pw_adv, dt, A, B, pdx_adv, pdy_adv, pdz_adv)

            ! for debugging
            IF (TRIM(particle_terminal) == "verbose") THEN
                WRITE(*,'("---------- Advection RK Step: ", I0, " ----------")') irk
                WRITE(*,'("Intermediate Velocity ", F12.9, " ", F12.9, " ", F12.9)') pu_adv, pv_adv, pw_adv
                WRITE(*,'("Intermediate (potential) Displacement ", F12.9, " ", F12.9, " ", F12.9)') pdx_adv, pdy_adv, pdz_adv
                WRITE(*, '()')
            END IF

            ! Particle Boundary Interaction
            CALL stop_timer(921)
            CALL start_timer(922)
            CALL move_particle(particle, pdx_adv, pdy_adv, pdz_adv, &
             pdx_eff, pdy_eff, pdz_eff, temp_coord, temp_grid)
            CALL stop_timer(922)
            CALL start_timer(921)
            pdx_pot = pdx_eff / B
            pdy_pot = pdy_eff / B
            pdz_pot = pdz_eff / B
            
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                ! for particle runtime statistics (terminal output)
                pd_eff_tot(1) = pd_eff_tot(1)+ pdx_eff
                pd_eff_tot(2) = pd_eff_tot(2)+ pdx_eff
                pd_eff_tot(3) = pd_eff_tot(3)+ pdx_eff
                psim_max_adv_dx = MAX(psim_max_adv_dx, ABS(pd_eff_tot(1)))
                psim_max_adv_dy = MAX(psim_max_adv_dy, ABS(pd_eff_tot(2)))
                psim_max_adv_dz = MAX(psim_max_adv_dz, ABS(pd_eff_tot(3)))
            END IF
        END DO
        CALL stop_timer(921)
    
    END SUBROUTINE particle_advection

    SUBROUTINE particle_diffusion(particle, temp_grid, temp_coord, pd_eff_tot, dt)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(inout) :: temp_grid
        REAL(realk), INTENT(inout) :: temp_coord(3)
        REAL(realk), INTENT(inout) :: pd_eff_tot(3)
        REAL(realk), INTENT(in) :: dt

        ! local variables
        REAL(realk) :: pdx_diff, pdy_diff, pdz_diff
        REAL(realk) :: pdx_eff, pdy_eff, pdz_eff

        IF (TRIM(particle_terminal) == "verbose") THEN
            WRITE(*,'("---------- Particle Diffusion ----------")')
            WRITE(*, '()')
        END IF

        CALL start_timer(924)
        CALL generate_diffusive_displacement(dt, D(1), D(2), D(3), pdx_diff, pdy_diff, pdz_diff)
        CALL stop_timer(924)

        CALL start_timer(925)
        CALL move_particle(particle, pdx_diff, pdy_diff, pdz_diff, &
             pdx_eff, pdy_eff, pdz_eff, temp_coord, temp_grid)
        CALL stop_timer(925)

        CALL start_timer(924)
        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
            ! for particle runtime statistics (terminal output)
            pd_eff_tot(1) = pd_eff_tot(1)+ pdx_eff
            pd_eff_tot(2) = pd_eff_tot(2)+ pdx_eff
            pd_eff_tot(3) = pd_eff_tot(3)+ pdx_eff
            psim_max_adv_dx = MAX(psim_max_adv_dx, ABS(pd_eff_tot(1)))
            psim_max_adv_dy = MAX(psim_max_adv_dy, ABS(pd_eff_tot(2)))
            psim_max_adv_dz = MAX(psim_max_adv_dz, ABS(pd_eff_tot(3)))
        END IF
        CALL stop_timer(924)

    END SUBROUTINE particle_diffusion

    SUBROUTINE finish_particle_timeintegration()

        CONTINUE

    END SUBROUTINE finish_particle_timeintegration

END MODULE
