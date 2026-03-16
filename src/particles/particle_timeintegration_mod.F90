MODULE particle_timeintegration_mod

    USE omp_lib

    USE fields_mod
    USE ib_mod
    USE gc_flowstencils_mod

    USE particle_config_mod
    USE particle_ofields_mod
    USE particle_runtimestat_mod
    USE particle_list_mod
    USE particle_interpolation_mod
    USE particle_diffusion_mod
    USE particle_exchange_mod

    IMPLICIT NONE

    TYPE(rk_2n_t) :: prkscheme

    ! rk coefficients as arrays (for offloading) TODO: remove
    INTEGER(intk) :: pnrk
    REAL(realk), ALLOCATABLE :: A_offload(:)
    REAL(realk), ALLOCATABLE :: B_offload(:) 

CONTAINS

    SUBROUTINE init_particle_timeintegration()

        USE gc_flowstencils_mod
        USE bound_flow_mod

        ! local variables
        INTEGER(intk) :: irk 

        CALL start_timer(900)
        CALL start_timer(910)

        ! init rk scheme
        CALL prkscheme%init(prkmethod)

        pnrk = prkscheme%nrk
        
        ALLOCATE(A_offload(pnrk))
        ALLOCATE(B_offload(pnrk))

        DO irk = 1, pnrk
            CALL prkscheme%get_coeffs(A_offload(irk), B_offload(irk), irk)
        END DO

        !$omp target enter data map(to: A_offload(1:pnrk), B_offload(1:pnrk))

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


    SUBROUTINE timeintegrate_particles_target(itstep, dt)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt

        ! local variables
        INTEGER(intk) :: dev_num, num_teams, num_threads
        INTEGER(intk) :: igrid, i, j, k, ipart, temp_grid
        INTEGER(intk) :: ii, jj, kk
        REAL(realk) :: temp_x, temp_y, temp_z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: pwu, pwv, pww
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:) :: obstacles

        INTEGER(intk) :: irk
        REAL(realk) :: pu_adv, pv_adv, pw_adv
        REAL(realk) :: pdx, pdy, pdz
        REAL(realk) :: pdx_pot, pdy_pot, pdz_pot
        REAL(realk) :: pdx_eff, pdy_eff, pdz_eff
        LOGICAL :: dreplace

        CALL start_timer(900)

        IF (dadvection) THEN
            !$omp target update to(u_offload, v_offload, w_offload)
        END IF
        
        !CALL write_particle_list_txt(itstep, 'pre')

#if defined __INTEL_COMPILER
        !$omp target update to(my_particle_list%particles)
#else
        !$omp target update to(my_particle_list)
#endif

        CALL start_timer(920)

        CALL count_pog(my_particle_list, grids_np, plist_displ)

        !$omp target update to(plist_displ(1:nmy_particle_grids))
        !$omp target update to(grids_np(1:nmy_particle_grids))

        dev_num = -99
        num_teams = -99
        num_threads = -99

#if defined __INTEL_COMPILER
        !$omp target map(tofrom: dev_num, num_teams, num_threads) map(mapper(obstacle_t), alloc: obstacles)
#else 
        !$omp target map(tofrom: dev_num, num_teams, num_threads)
#endif
        !$omp teams distribute private(igrid, ipart, temp_grid, ii, jj, kk, temp_x, temp_y, temp_z, &
        !$omp irk, pu_adv, pv_adv, pw_adv, pdx, pdy, pdz, pdx_pot, pdy_pot, pdz_pot, pdx_eff, pdy_eff, pdz_eff, &
        !$omp x, y, z, dx, dy, dz, ddx, ddy, ddz, pwu, pwv, pww, obstacles) reduction(max: num_threads)
        DO i = 1, nmy_particle_grids

#ifdef __GFORTRAN__
            !$omp master
                dev_num = omp_get_device_num()
                num_teams = omp_get_num_teams()
            !$omp end master
#endif
#ifdef __INTEL_COMPILER
            dev_num = omp_get_device_num()
            num_teams = omp_get_num_teams()
#endif

            igrid = my_particle_grids(i)

            CALL get_mgdims_target(kk, jj, ii, igrid)
            
            CALL ptr_to_grid_1(x_offload, igrid, x, 1_intk)
            CALL ptr_to_grid_1(y_offload, igrid, y, 2_intk)
            CALL ptr_to_grid_1(z_offload, igrid, z, 3_intk)

            CALL ptr_to_grid_1(dx_offload, igrid, dx, 1_intk)
            CALL ptr_to_grid_1(dy_offload, igrid, dy, 2_intk)
            CALL ptr_to_grid_1(dz_offload, igrid, dz, 3_intk)
            
            CALL ptr_to_grid_1(ddx_offload, igrid, ddx, 1_intk)
            CALL ptr_to_grid_1(ddy_offload, igrid, ddy, 2_intk)
            CALL ptr_to_grid_1(ddz_offload, igrid, ddz, 3_intk)

            CALL ptr_to_grid_3(u_offload, igrid, pwu)
            CALL ptr_to_grid_3(v_offload, igrid, pwv)
            CALL ptr_to_grid_3(w_offload, igrid, pww)

            obstacles => my_obstacles_offload(obstacle_displ(igrid) + 1: obstacle_displ(igrid) + MAX(1_intk, n_my_obstacles_on_grid(igrid)))
            
            !$omp parallel do private(ipart, temp_grid, temp_x, temp_y, temp_z) firstprivate(igrid, ii, jj, kk, &
            !$omp irk, pu_adv, pv_adv, pw_adv, pdx, pdy, pdz, pdx_pot, pdy_pot, pdz_pot, pdx_eff, pdy_eff, pdz_eff) &
            !$omp shared(x, y, z, dx, dy, dz, ddx, ddy, ddz, pwu, pwv, pww, obstacles)
            DO j = 1, grids_np(i)

                num_threads = omp_get_num_threads()
                
                ipart = plist_displ(i) + j

                temp_grid = my_particle_list%particles(ipart)%igrid
                temp_x = my_particle_list%particles(ipart)%x
                temp_y = my_particle_list%particles(ipart)%y
                temp_z = my_particle_list%particles(ipart)%z

                IF (dadvection) THEN

                    DO irk = 1, pnrk

                        ! get particle velocity
                        CALL interpolate_lincon(my_particle_list%particles(ipart), kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                        pwu, pwv, pww, pu_adv, pv_adv, pw_adv)

                        CALL prkstep(pdx_pot, pdy_pot, pdz_pot, pu_adv, pv_adv, pw_adv, dt, &
                        A_offload(irk), B_offload(irk), pdx, pdy, pdz)

                        ! Particle Boundary Interaction
                        CALL move_particle_target(my_particle_list%particles(ipart), pdx, pdy, pdz, &
                        pdx_eff, pdy_eff, pdz_eff, temp_x, temp_y, temp_z, temp_grid, particle_boundaries, obstacles, dreplace)
                        
                        IF (dreplace) THEN
                            CALL replace_particle_target(my_particle_list%particles(ipart), obstacles, kk, jj, ii, x, y, z, dx, dy, dz)
                        ELSE 
                            CALL update_particle_cell_target(my_particle_list%particles(ipart), kk, jj, ii, x, y, z, dx, dy, dz)
                        END IF

                        pdx_pot = pdx_eff / B_offload(irk)
                        pdy_pot = pdy_eff / B_offload(irk)
                        pdz_pot = pdz_eff / B_offload(irk)

                    ! TODO: reintroduce particle runtime statistics
                    END DO

                END IF

#ifdef _MGLET_OPENMP_
                IF (ddiffusion) THEN
                    CALL generate_diffusive_displacement_target(dt, D(1), D(2), D(3), pdx, pdy, pdz, my_particle_list%particles(ipart)%seed)

                    CALL move_particle_target(my_particle_list%particles(ipart), pdx, pdy, pdz, &
                        pdx_eff, pdy_eff, pdz_eff, temp_x, temp_y, temp_z, temp_grid, particle_boundaries, obstacles, dreplace)

                    IF (dreplace) THEN
                        CALL replace_particle_target(my_particle_list%particles(ipart), obstacles, kk, jj, ii, x, y, z, dx, dy, dz)
                    ELSE 
                        CALL update_particle_cell_target(my_particle_list%particles(ipart), kk, jj, ii, x, y, z, dx, dy, dz)
                    END IF
                END IF
#endif

                ! TODO: reintroduce particle runtime statistics
            END DO
            !$omp end parallel do 
        END DO
        !$omp end teams distribute
        !$omp end target

        CALL stop_timer(920)
        
#if defined __INTEL_COMPILER
        !$omp target update from(my_particle_list%particles)
#else
        !$omp target update from(my_particle_list)
#endif

        !CALL write_particle_list_txt(itstep, 'pos')

        CALL stop_timer(900)

        WRITE(*, '("    Timeintegration on Process:                     ", I9)') myid
        WRITE(*, '("        Device Number:                              ", I9)') dev_num
        WRITE(*, '("        Number of Teams:                            ", I9)') num_teams
        WRITE(*, '("        Max. Number of Threds (per Team):           ", I9)') num_threads
    
    END SUBROUTINE timeintegrate_particles_target

    SUBROUTINE particle_advection_target(particle, temp_grid, temp_x, temp_y, temp_z, kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                                         pwu, pwv, pww, dt, pnrk, A, B, boundaries, obstacles)

        !$omp declare target

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(inout) :: temp_grid
        REAL(realk), INTENT(inout) :: temp_x, temp_y, temp_z
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :), INTENT(in) :: pwu, pwv, pww
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: pnrk
        REAL(realk), INTENT(in) :: A(pnrk), B(pnrk)
        TYPE(particle_boundaries_t), INTENT(in) :: boundaries(ngrid)
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles
        
        
        ! local variables
        INTEGER(intk) :: irk
        REAL(realk) :: pu_adv, pv_adv, pw_adv
        REAL(realk) :: pdx_adv, pdy_adv, pdz_adv
        REAL(realk) :: pdx_pot, pdy_pot, pdz_pot
        REAL(realk) :: pdx_eff, pdy_eff, pdz_eff
        LOGICAL :: dreplace

        DO irk = 1, pnrk

            ! get particle velocity
            CALL interpolate_lincon(particle, kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
             pwu, pwv, pww, pu_adv, pv_adv, pw_adv)

            CALL prkstep(pdx_pot, pdy_pot, pdz_pot, pu_adv, pv_adv, pw_adv, dt, &
             A(irk), B(irk), pdx_adv, pdy_adv, pdz_adv)

            ! Particle Boundary Interaction
            CALL move_particle_target(particle, pdx_adv, pdy_adv, pdz_adv, &
             pdx_eff, pdy_eff, pdz_eff, temp_x, temp_y, temp_z, temp_grid, boundaries, obstacles, dreplace)
            
            IF (dreplace) THEN
                CALL replace_particle_target(particle, obstacles, kk, jj, ii, x, y, z, dx, dy, dz)
            ELSE 
                CALL update_particle_cell_target(particle, kk, jj, ii, x, y, z, dx, dy, dz)
            END IF

            pdx_pot = pdx_eff / B(irk)
            pdy_pot = pdy_eff / B(irk)
            pdz_pot = pdz_eff / B(irk)

        ! TODO: reintroduce particle runtime statistics
        END DO

    END SUBROUTINE particle_advection_target

    SUBROUTINE particle_diffusion_target(particle, temp_grid, temp_x, temp_y, temp_z, kk, jj, ii, x, y, z, dx, dy, dz, dt, seed, boundaries, obstacles)

        !$omp declare target

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(inout) :: temp_grid
        REAL(realk), INTENT(inout) :: temp_x, temp_y, temp_z
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: dx, dy, dz
        REAL(realk), INTENT(in) :: dt
        INTEGER(c_int), INTENT(inout) :: seed
        TYPE(particle_boundaries_t), INTENT(in) :: boundaries(ngrid)
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles

        ! local variables
        REAL(realk) :: pdx_diff, pdy_diff, pdz_diff
        REAL(realk) :: pdx_eff, pdy_eff, pdz_eff
        LOGICAL :: dreplace

        CALL generate_diffusive_displacement_target(dt, D(1), D(2), D(3), pdx_diff, pdy_diff, pdz_diff, seed)

        CALL move_particle_target(particle, pdx_diff, pdy_diff, pdz_diff, &
             pdx_eff, pdy_eff, pdz_eff, temp_x, temp_y, temp_z, temp_grid, boundaries, obstacles, dreplace)

        IF (dreplace) THEN
            CALL replace_particle_target(particle, obstacles, kk, jj, ii, x, y, z, dx, dy, dz)
        ELSE 
            CALL update_particle_cell_target(particle, kk, jj, ii, x, y, z, dx, dy, dz)
        END IF

        ! TODO: reintroduce particle runtime statistics

    END SUBROUTINE particle_diffusion_target


    SUBROUTINE finish_particle_timeintegration()

        !$omp target exit data map(delete: A_offload, B_offload)

        DEALLOCATE(A_offload)
        DEALLOCATE(B_offload)

    END SUBROUTINE finish_particle_timeintegration

END MODULE
