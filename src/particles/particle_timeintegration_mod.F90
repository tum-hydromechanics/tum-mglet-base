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
    USE particle_boundaries_mod
    
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
                    CALL print_particle_status(my_particle_list%particles(ipart))
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


    SUBROUTINE timeintegrate_particles_target3(itstep, dt)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt

        ! local variables (fixed)
        INTEGER(intk) :: dev_num, num_teams, num_threads
        INTEGER(intk) :: igrid, icorn, ipart, i, j, k, irk
        INTEGER(intk) :: ip, ii, jj, kk
        INTEGER(intk) :: pstag(3)
        REAL(realk) :: pdx_pot, pdy_pot, pdz_pot
        REAL(realk) :: dvec(3), cvec(3), n(3)
        REAL(realk) :: bbox(6)

        ! local variables (flexible)
        INTEGER(intk) :: l, int_var_1, int_var_2, int_var_3
        REAL(realk) :: real_var_1, real_var_2, real_var_3, real_var_4, real_var_5
        REAL(realk) :: real_var_6, real_var_7, real_var_8, real_var_9, real_var_10
        REAL(realk) :: real_arr_1(3)

        CALL start_timer(900)

        CALL start_timer(920)

        CALL start_timer(921)

        IF (dadvection) THEN
            !$omp target update to(u_offload, v_offload, w_offload)
        END IF
        
        CALL stop_timer(921)

        !CALL write_particle_list_txt(itstep, 'pre')

        CALL start_timer(922)

#if defined __INTEL_COMPILER
        !$omp target update to(my_particle_list%particles)
#else
        !$omp target update to(my_particle_list)
#endif

        CALL stop_timer(922)

        CALL start_timer(923)

        CALL count_pog(my_particle_list, grids_np, plist_displ)

        !$omp target update to(plist_displ(1:nmy_particle_grids))
        !$omp target update to(grids_np(1:nmy_particle_grids))

        CALL stop_timer(923)

        dev_num = -99
        num_teams = -99
        num_threads = -99

        CALL start_timer(924)

        !$omp target map(tofrom: dev_num, num_teams, num_threads)
        !$omp teams distribute private(igrid, ipart, icorn, pstag, ip, ii, jj, kk, &
        !$omp irk, k, l, dvec, pdx_pot, pdy_pot, pdz_pot, &
        !$omp bbox, cvec, n, int_var_1, int_var_2, int_var_3, real_arr_1, &
        !$omp real_var_1, real_var_2, real_var_3, real_var_4, real_var_5, real_var_6, real_var_7, real_var_8, real_var_9, real_var_10) &
        !$omp reduction(max: num_threads)
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

            ip = ip1d_offload(igrid)
            
            CALL get_bbox_target(bbox(1), bbox(2), bbox(3), bbox(4), bbox(5), bbox(6), igrid)
            
            !$omp parallel do private(ipart, icorn, pstag, irk, cvec, dvec, n, k, l, int_var_1, int_var_2, int_var_3, real_arr_1, &
            !$omp real_var_1, real_var_2, real_var_3, real_var_4, real_var_5, real_var_6, real_var_7, real_var_8, real_var_9, real_var_10) &
            !$omp firstprivate(igrid, ip, ii, jj, kk, pdx_pot, pdy_pot, pdz_pot, bbox)
            DO j = 1, grids_np(i)

                !number of threads in the team (working on j-loop)
                num_threads = omp_get_num_threads()
                
                ipart = plist_displ(i) + j
                
                CALL get_particle_gcorner_target(my_particle_list%particles(ipart), bbox, icorn)

                pstag = 0_intk

                IF (dadvection) THEN
                    
                    pdx_pot = 0.0
                    pdy_pot = 0.0
                    pdz_pot = 0.0

                    DO irk = 1, pnrk

                        ! >>> begin local "conceptual namespace" <<<
                        ! real_var_1 => pu_adv
                        ! real_var_2 => pv_adv
                        ! real_var_3 => pw_adv

                        ! get particle velocity
                        CALL interpolate_lincon_target(my_particle_list%particles(ipart), igrid, kk, jj, ii, real_var_1, real_var_2, real_var_3)

                        CALL prkstep(pdx_pot, pdy_pot, pdz_pot, real_var_1, real_var_2, real_var_3, dt, &
                        A_offload(irk), B_offload(irk), dvec(1), dvec(2), dvec(3))
                        
                        ! >>> end local "conceptual namespace" <<<

                        ! >>> begin local "conceptual namespace" <<<
                        ! int_var_1 => idir
                        ! int_var_2 => iobst_local_new
                        ! int_var_3 => iobst_local_old
                        ! real_var_1 => s
                        ! real_var_2 => temp 
                        ! real_var_3-real_var_10 nested

                        cvec(1) = my_particle_list%particles(ipart)%x
                        cvec(2) = my_particle_list%particles(ipart)%y
                        cvec(3) = my_particle_list%particles(ipart)%z

                        int_var_2 = 0

                        ! to avoid branch divergence here, just iterate to the max. number of iterations that would be a stoping criterion anyways
                        DO k = 1, 10

                            int_var_1 = 0
                            int_var_3 = int_var_2
                            int_var_2 = 0

                            real_var_1 = 1.0

                            ! STEP 1 - OBSTACLES
                            IF (n_my_obstacles_on_grid(igrid) > 0) THEN
                                ! >>> begin local "conceptual namespace" <<<
                                ! real_var_3  => a
                                ! real_var_4  => b0
                                ! real_var_5  => c0
                                ! real_var_6  => sa
                                ! real_var_7  => sb
                                ! real_var_8-real_var_10 nested

                                ! first coefficient
                                real_var_3 = (dvec(1)**2 + dvec(2)**2 + dvec(3)**2)

                                IF (.NOT. (dvec(1)**2 + dvec(2)**2 + dvec(3)**2 == 0.0)) THEN
                                
                                    real_var_4 = 2*cvec(1)*dvec(1) + 2*cvec(2)*dvec(2) + 2*cvec(3)*dvec(3)
                                    real_var_5 = cvec(1)**2 + cvec(2)**2 + cvec(3)**2

                                    ! iterate over all obstacles of the grid
                                    DO l = 1, n_my_obstacles_on_grid(igrid)

                                        ! check if a particle interacts with the obstacle it has been deflected from in the previous timestep
                                        IF (l == int_var_3 .OR. my_obstacles_offload(obstacle_displ(igrid) + l)%iobst < 0) THEN
                                            CYCLE
                                        END IF

                                        ! >>> begin local "conceptual namespace" <<<
                                        ! real_var_8  => b
                                        ! real_var_9  => c
                                        ! real_var_10 => d

                                        real_var_8 = real_var_4 - &
                                            2*my_obstacles_offload(obstacle_displ(igrid) + l)%x*dvec(1) - &
                                            2*my_obstacles_offload(obstacle_displ(igrid) + l)%y*dvec(2) - &
                                            2*my_obstacles_offload(obstacle_displ(igrid) + l)%z*dvec(3)
                                        real_var_9 = real_var_5 + &
                                            my_obstacles_offload(obstacle_displ(igrid) + l)%x**2 + &
                                            my_obstacles_offload(obstacle_displ(igrid) + l)%y**2 + &
                                            my_obstacles_offload(obstacle_displ(igrid) + l)%z**2 - &
                                            2*cvec(obstacle_displ(igrid) + l)*my_obstacles_offload(obstacle_displ(igrid) + l)%x - &
                                            2*cvec(obstacle_displ(igrid) + l)*my_obstacles_offload(obstacle_displ(igrid) + l)%y - &
                                            2*cvec(obstacle_displ(igrid) + l)*my_obstacles_offload(obstacle_displ(igrid) + l)%z - &
                                            my_obstacles_offload(i)%radius**2
                                        real_var_10 = real_var_8**2 - 4*real_var_3*real_var_9

                                        IF (real_var_10 < EPSILON(0.0_realk)) THEN
                                            CYCLE
                                        END IF

                                        real_var_6 = (-real_var_8 + SQRT(real_var_10)) / (2 * real_var_3)
                                        real_var_7 = (-real_var_8 - SQRT(real_var_10)) / (2 * real_var_3)
                                        ! >>> end local "conceptual namespace" <<<

                                        ! >>> begin local "conceptual namespace" <<<
                                        ! real_var_8  => sc
                                        ! real_var_9  => sd

                                        ! if a particle moves towards an obstacle, limit its motion to the closest intersection yet
                                        IF (real_var_6 >= 0.0 .AND. real_var_7 >= 0.0) THEN
                                            real_var_8 = MIN(real_var_6, real_var_7)
                                            IF (real_var_8 < real_var_1) THEN
                                                real_var_1 = real_var_8
                                                int_var_2 = i
                                            END IF
                                        ELSEIF (real_var_6 <= 0.0 .AND. real_var_7 <= 0.0) THEN
                                            CYCLE
                                        ELSE
                                            real_var_8 = MIN(real_var_6, real_var_7)
                                            real_var_9 = MAX(real_var_6, real_var_7)

                                            IF (ABS(real_var_8) < ABS(real_var_9)) THEN
                                                real_var_1 = 0.0
                                                int_var_2 = i
                                                EXIT
                                            ELSEIF (ABS(real_var_8) >= ABS(real_var_9)) THEN
                                                CYCLE
                                            END IF

                                        END IF
                                        ! >>> end local "conceptual namespace" <<<
                                    END DO
                                    ! >>> end local "conceptual namespace" <<<
                                END IF
                            END IF

                            ! STEP 2 - GRID BOUNDARIES
                            IF (real_var_1 > 0.0_realk) THEN
                                ! >>> begin local "conceptual namespace" <<<
                                ! real_arr_1(1)  => l_a(1)
                                ! real_arr_1(2)  => l_a(2)
                                ! real_arr_1(3)  => l_a(3)
                                ! real_var_6  => rx
                                ! real_var_7  => ry
                                ! real_var_8  => rz
                                ! real_var_9  => rmax
                                ! real_var_10 => ratio

                                real_var_9 = 1.0

                                ! abs distance of particle to grid boundaries
                                real_arr_1(1) = ABS(cvec(1) - particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%face_coord(1))
                                IF (real_arr_1(1) < EPSILON(real_arr_1(1))) THEN
                                    real_var_6 = SIGN(HUGE(real_var_6), particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(1) * ((-1.0) ** pstag(1)) * dvec(1)) * ABS(dvec(1))
                                ELSE
                                    real_var_6 = particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(1) * ((-1_intk) ** pstag(1)) * ((real_var_1 * dvec(1)) / (real_arr_1(1)))
                                END IF

                                IF (real_var_9 <= real_var_6) THEN 
                                    real_var_9 = real_var_6
                                    int_var_1 = 1
                                END IF

                                real_arr_1(2) = ABS(cvec(2) - particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%face_coord(2))
                                IF (real_arr_1(2) < EPSILON(real_arr_1(2))) THEN
                                    real_var_7 = SIGN(HUGE(real_var_7), particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(2) * ((-1.0) ** pstag(2)) * dvec(2)) * ABS(dvec(2))
                                ELSE
                                    real_var_7 = particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(2) * ((-1_intk) ** pstag(2)) * ((real_var_1 * dvec(2)) / (real_arr_1(2)))
                                END IF

                                IF (real_var_9 <= real_var_7) THEN 
                                    real_var_9 = real_var_7
                                    int_var_1 = 2
                                END IF

                                real_arr_1(3) = ABS(cvec(3) - particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%face_coord(3))
                                IF (real_arr_1(3) < EPSILON(real_arr_1(3))) THEN
                                    real_var_8 = SIGN(HUGE(real_var_8), particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(3) * ((-1.0) ** pstag(3)) * dvec(3)) * ABS(dvec(3))
                                ELSE
                                    real_var_8 = particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(3) * ((-1_intk) ** pstag(3)) * ((real_var_1 * dvec(3)) / (real_arr_1(3)))
                                END IF

                                IF (real_var_9 <= real_var_8) THEN 
                                    real_var_9 = real_var_8
                                    int_var_1 = 3
                                END IF
                                
                                IF (int_var_1 == 0) THEN
                                    cvec(1) = cvec(1) + dvec(1) * real_var_1
                                    dvec(1) = dvec(1) - dvec(1) * real_var_1
                                    cvec(2) = cvec(2) + dvec(2) * real_var_1
                                    dvec(2) = dvec(2) - dvec(2) * real_var_1
                                    cvec(3) = cvec(3) + dvec(3) * real_var_1
                                    dvec(3) = dvec(3) - dvec(3) * real_var_1
                                ELSE
                                    int_var_2 = 0_intk

                                    real_var_10 = real_arr_1(int_var_1) / ABS(dvec(int_var_1))

                                    cvec(int_var_1) = particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%face_coord(int_var_1) 
                                    dvec(int_var_1) = dvec(int_var_1) - particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(int_var_1) * real_arr_1(int_var_1)
                                    l = MOD(int_var_1, 3) + 1
                                    cvec(l) = cvec(l) + (real_var_10 * dvec(l)) - EPSILON(cvec(l)) * particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(l)
                                    dvec(l) = dvec(l) - (real_var_10 * dvec(l))
                                    l = MOD(int_var_1 + 1, 3) + 1
                                    cvec(l) = cvec(l) + (real_var_10 * dvec(l)) - EPSILON(cvec(l)) * particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(l)
                                    dvec(l) = dvec(l) - (real_var_10 * dvec(l))
                                END IF
                                ! >>> end local "conceptual namespace" <<<
                            END IF

                            IF (0 < int_var_2) THEN

                                ! reflect at obstacle
                                ! compute normal vector
                                n(1) = cvec(1) - my_obstacles_offload(obstacle_displ(igrid) + int_var_2)%x
                                n(2) = cvec(2) - my_obstacles_offload(obstacle_displ(igrid) + int_var_2)%y
                                n(3) = cvec(3) - my_obstacles_offload(obstacle_displ(igrid) + int_var_2)%z

                                ! magnitude
                                real_var_2 = SQRT(n(1)**2 + n(2)**2 + n(3)**2)

                                n(1) = n(1) / real_var_2
                                n(2) = n(2) / real_var_2
                                n(3) = n(3) / real_var_2

                                ! alter displacement verctor
                                ! dot product
                                real_var_2 = MIN((n(1) * dvec(1) + n(2) * dvec(2) + n(3) * dvec(3)), 0.0)

                                dvec(1) = dvec(1) - 2 * real_var_2 * n(1)
                                dvec(2) = dvec(2) - 2 * real_var_2 * n(2)
                                dvec(3) = dvec(3) - 2 * real_var_2 * n(3)

                            ELSEIF (0 < int_var_1) THEN
                                
                                CALL get_gcorner_normal(particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn), pstag, int_var_1, n)
                                
                                ! reflect at grid boundary
                                ! dot product
                                real_var_2 = MIN((n(1) * dvec(1) + n(2) * dvec(2) + n(3) * dvec(3)), 0.0)

                                dvec(1) = dvec(1) - 2 * real_var_2 * n(1)
                                dvec(2) = dvec(2) - 2 * real_var_2 * n(2)
                                dvec(3) = dvec(3) - 2 * real_var_2 * n(3)

                                !update pstag (normal vector idir (int_var_1) component must be zero or point inwards for this method to work)
                                pstag(int_var_1) = MAX(pstag(int_var_1) + (1 - pstag(int_var_1)) + &
                                 NINT(n(int_var_1) * particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(int_var_1)), 0)
                            END IF

                        END DO

                        pdx_pot = (cvec(1) - my_particle_list%particles(ipart)%x) / B_offload(irk)
                        pdy_pot = (cvec(2) - my_particle_list%particles(ipart)%y) / B_offload(irk)
                        pdz_pot = (cvec(3) - my_particle_list%particles(ipart)%z) / B_offload(irk)

                        !my_particle_list%particles(ipart)%xyz_abs(1) = my_particle_list%particles(ipart)%xyz_abs(1) + (cvec(1) - my_particle_list%particles(ipart)%x)
                        !my_particle_list%particles(ipart)%xyz_abs(2) = my_particle_list%particles(ipart)%xyz_abs(2) + (cvec(2) - my_particle_list%particles(ipart)%y)
                        !my_particle_list%particles(ipart)%xyz_abs(3) = my_particle_list%particles(ipart)%xyz_abs(3) + (cvec(3) - my_particle_list%particles(ipart)%z)

                        my_particle_list%particles(ipart)%x = cvec(1)
                        my_particle_list%particles(ipart)%y = cvec(2)
                        my_particle_list%particles(ipart)%z = cvec(3)

                        CALL update_particle_cell_target(my_particle_list%particles(ipart), kk, jj, ii, &
                        x_offload(ip:ip+ii-1), y_offload(ip:ip+jj-1), z_offload(ip:ip+kk-1), &
                        dx_offload(ip:ip+ii-1), dy_offload(ip:ip+jj-1), dz_offload(ip:ip+kk-1))

                        ! >>> end local "conceptual namespace" <<<
                    
                    ! TODO: reintroduce particle runtime statistics
                    END DO

                END IF

#ifdef _MGLET_OPENMP_
                IF (ddiffusion) THEN
                    CALL generate_diffusive_displacement_target(dt, D(1), D(2), D(3), dvec(1), dvec(2), dvec(3), my_particle_list%particles(ipart)%seed)

                    ! >>> end local "conceptual namespace" <<<

                    ! >>> begin local "conceptual namespace" <<<
                    ! int_var_1 => idir
                    ! int_var_2 => iobst_local_new
                    ! int_var_3 => iobst_local_old
                    ! real_var_1 => s
                    ! real_var_2 => temp 

                    cvec(1) = my_particle_list%particles(ipart)%x
                    cvec(2) = my_particle_list%particles(ipart)%y
                    cvec(3) = my_particle_list%particles(ipart)%z

                    int_var_2 = 0

                    ! to avoid branch divergence here, just iterate to the max. number of iterations that would be a stoping criterion anyways
                    DO k = 1, 10

                        int_var_1 = 0
                        int_var_3 = int_var_2
                        int_var_2 = 0

                        real_var_1 = 1.0

                        ! STEP 1 - OBSTACLES
                        IF (n_my_obstacles_on_grid(igrid) > 0) THEN
                            ! >>> begin local "conceptual namespace" <<<
                            ! real_var_3  => a
                            ! real_var_4  => b0
                            ! real_var_5  => c0
                            ! real_var_6  => sa
                            ! real_var_7  => sb
                            ! real_var_8-real_var_10 nested

                            ! first coefficient
                            real_var_3 = (dvec(1)**2 + dvec(2)**2 + dvec(3)**2)

                            IF (.NOT. (dvec(1)**2 + dvec(2)**2 + dvec(3)**2 == 0.0)) THEN
                            
                                real_var_4 = 2*cvec(1)*dvec(1) + 2*cvec(2)*dvec(2) + 2*cvec(3)*dvec(3)
                                real_var_5 = cvec(1)**2 + cvec(2)**2 + cvec(3)**2

                                ! iterate over all obstacles of the grid
                                DO l = 1, n_my_obstacles_on_grid(igrid)

                                    ! check if a particle interacts with the obstacle it has been deflected from in the previous timestep
                                    IF (l == int_var_3 .OR. my_obstacles_offload(obstacle_displ(igrid) + l)%iobst < 0) THEN
                                        CYCLE
                                    END IF

                                    ! >>> begin local "conceptual namespace" <<<
                                    ! real_var_8  => b
                                    ! real_var_9  => c
                                    ! real_var_10 => d

                                    real_var_8 = real_var_4 - &
                                        2*my_obstacles_offload(obstacle_displ(igrid) + l)%x*dvec(1) - &
                                        2*my_obstacles_offload(obstacle_displ(igrid) + l)%y*dvec(2) - &
                                        2*my_obstacles_offload(obstacle_displ(igrid) + l)%z*dvec(3)
                                    real_var_9 = real_var_5 + &
                                        my_obstacles_offload(obstacle_displ(igrid) + l)%x**2 + &
                                        my_obstacles_offload(obstacle_displ(igrid) + l)%y**2 + &
                                        my_obstacles_offload(obstacle_displ(igrid) + l)%z**2 - &
                                        2*cvec(1)*my_obstacles_offload(obstacle_displ(igrid) + l)%x - &
                                        2*cvec(2)*my_obstacles_offload(obstacle_displ(igrid) + l)%y - &
                                        2*cvec(3)*my_obstacles_offload(obstacle_displ(igrid) + l)%z - &
                                        my_obstacles_offload(obstacle_displ(igrid) + l)%radius**2
                                    real_var_10 = real_var_8**2 - 4*real_var_3*real_var_9

                                    IF (real_var_10 < EPSILON(0.0_realk)) THEN
                                        CYCLE
                                    END IF

                                    real_var_6 = (-real_var_8 + SQRT(real_var_10)) / (2 * real_var_3)
                                    real_var_7 = (-real_var_8 - SQRT(real_var_10)) / (2 * real_var_3)
                                    ! >>> end local "conceptual namespace" <<<

                                    ! >>> begin local "conceptual namespace" <<<
                                    ! real_var_8  => sc
                                    ! real_var_9  => sd

                                    ! if a particle moves towards an obstacle, limit its motion to the closest intersection yet
                                    IF (real_var_6 >= 0.0 .AND. real_var_7 >= 0.0) THEN
                                        real_var_8 = MIN(real_var_6, real_var_7)
                                        IF (real_var_8 < real_var_1) THEN
                                            real_var_1 = real_var_8
                                            int_var_2 = i
                                        END IF
                                    ELSEIF (real_var_6 <= 0.0 .AND. real_var_7 <= 0.0) THEN
                                        CYCLE
                                    ELSE
                                        real_var_8 = MIN(real_var_6, real_var_7)
                                        real_var_9 = MAX(real_var_6, real_var_7)

                                        IF (ABS(real_var_8) < ABS(real_var_9)) THEN
                                            real_var_1 = 0.0
                                            int_var_2 = i
                                            EXIT
                                        ELSEIF (ABS(real_var_8) >= ABS(real_var_9)) THEN
                                            CYCLE
                                        END IF

                                    END IF
                                    ! >>> end local "conceptual namespace" <<<
                                END DO
                                ! >>> end local "conceptual namespace" <<<
                            END IF
                        END IF

                        ! STEP 2 - GRID BOUNDARIES
                        IF (real_var_1 > 0.0_realk) THEN
                            ! >>> begin local "conceptual namespace" <<<
                            ! real_arr_1(1)  => l_a(1)
                            ! real_arr_1(2)  => l_a(2)
                            ! real_arr_1(3)  => l_a(3)
                            ! real_var_6  => rx
                            ! real_var_7  => ry
                            ! real_var_8  => rz
                            ! real_var_9  => rmax
                            ! real_var_10 => ratio

                            int_var_1 = 0
                            real_var_9 = 1.0

                            ! abs distance of particle to grid boundaries
                            real_arr_1(1) = ABS(cvec(1) - particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%face_coord(1))
                            IF (real_arr_1(1) < EPSILON(real_arr_1(1))) THEN
                                real_var_6 = SIGN(HUGE(real_var_6), particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(1) * ((-1.0) ** pstag(1)) * dvec(1)) * ABS(dvec(1))
                            ELSE
                                real_var_6 = particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(1) * ((-1_intk) ** pstag(1)) * ((real_var_1 * dvec(1)) / (real_arr_1(1)))
                            END IF

                            IF (real_var_9 <= real_var_6) THEN 
                                real_var_9 = real_var_6
                                int_var_1 = 1
                            END IF

                            real_arr_1(2) = ABS(cvec(2) - particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%face_coord(2))
                            IF (real_arr_1(2) < EPSILON(real_arr_1(2))) THEN
                                real_var_7 = SIGN(HUGE(real_var_7), particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(2) * ((-1.0) ** pstag(2)) * dvec(2)) * ABS(dvec(2))
                            ELSE
                                real_var_7 = particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(2) * ((-1_intk) ** pstag(2)) * ((real_var_1 * dvec(2)) / (real_arr_1(2)))
                            END IF

                            IF (real_var_9 <= real_var_7) THEN 
                                real_var_9 = real_var_7
                                int_var_1 = 2
                            END IF

                            real_arr_1(3) = ABS(cvec(3) - particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%face_coord(3))
                            IF (real_arr_1(3) < EPSILON(real_arr_1(3))) THEN
                                real_var_8 = SIGN(HUGE(real_var_8), particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(3) * ((-1.0) ** pstag(3)) * dvec(3)) * ABS(dvec(3))
                            ELSE
                                real_var_8 = particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(3) * ((-1_intk) ** pstag(3)) * ((real_var_1 * dvec(3)) / (real_arr_1(3)))
                            END IF

                            IF (real_var_9 <= real_var_8) THEN 
                                real_var_9 = real_var_8
                                int_var_1 = 3
                            END IF

                            IF (int_var_1 == 0) THEN
                                cvec(1) = cvec(1) + dvec(1) * real_var_1
                                dvec(1) = dvec(1) - dvec(1) * real_var_1
                                cvec(2) = cvec(2) + dvec(2) * real_var_1
                                dvec(2) = dvec(2) - dvec(2) * real_var_1
                                cvec(3) = cvec(3) + dvec(3) * real_var_1
                                dvec(3) = dvec(3) - dvec(3) * real_var_1
                            ELSE
                                int_var_2 = 0_intk

                                real_var_10 = real_arr_1(int_var_1) / ABS(dvec(int_var_1))

                                cvec(int_var_1) = particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%face_coord(int_var_1) 
                                dvec(int_var_1) = dvec(int_var_1) - particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(int_var_1) * real_arr_1(int_var_1)
                                l = MOD(int_var_1, 3) + 1
                                cvec(l) = cvec(l) + (real_var_10 * dvec(l)) - EPSILON(cvec(l)) * particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(l)
                                dvec(l) = dvec(l) - (real_var_10 * dvec(l))
                                l = MOD(int_var_1 + 1, 3) + 1
                                cvec(l) = cvec(l) + (real_var_10 * dvec(l)) - EPSILON(cvec(l)) * particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(l)
                                dvec(l) = dvec(l) - (real_var_10 * dvec(l))
                            END IF
                            ! >>> end local "conceptual namespace" <<<
                        END IF

                        IF (0 < int_var_2) THEN

                            ! reflect at obstacle
                            ! compute normal vector
                            n(1) = cvec(1) - my_obstacles_offload(obstacle_displ(igrid) + int_var_2)%x
                            n(2) = cvec(2) - my_obstacles_offload(obstacle_displ(igrid) + int_var_2)%y
                            n(3) = cvec(3) - my_obstacles_offload(obstacle_displ(igrid) + int_var_2)%z

                            ! magnitude
                            real_var_2 = SQRT(n(1)**2 + n(2)**2 + n(3)**2)

                            n(1) = n(1) / real_var_2
                            n(2) = n(2) / real_var_2
                            n(3) = n(3) / real_var_2

                            ! alter displacement verctor
                            ! dot product
                            real_var_2 = MIN((n(1) * dvec(1) + n(2) * dvec(2) + n(3) * dvec(3)), 0.0)

                            dvec(1) = dvec(1) - 2 * real_var_2 * n(1)
                            dvec(2) = dvec(2) - 2 * real_var_2 * n(2)
                            dvec(3) = dvec(3) - 2 * real_var_2 * n(3)

                        ELSEIF (0 < int_var_1) THEN
                            
                            CALL get_gcorner_normal(particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn), pstag, int_var_1, n)
                            
                            ! reflect at grid boundary
                            ! dot product
                            real_var_2 = MIN((n(1) * dvec(1) + n(2) * dvec(2) + n(3) * dvec(3)), 0.0)

                            dvec(1) = dvec(1) - 2 * real_var_2 * n(1)
                            dvec(2) = dvec(2) - 2 * real_var_2 * n(2)
                            dvec(3) = dvec(3) - 2 * real_var_2 * n(3)

                            !update pstag (normal vector idir (int_var_1) component must be zero or point inwards for this method to work)
                            pstag(int_var_1) = MAX(pstag(int_var_1) + (1 - pstag(int_var_1)) + &
                             NINT(n(int_var_1) * particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn)%location(int_var_1)), 0)
                        END IF

                    END DO

                    pdx_pot = (cvec(1) - my_particle_list%particles(ipart)%x) / B_offload(irk)
                    pdy_pot = (cvec(2) - my_particle_list%particles(ipart)%y) / B_offload(irk)
                    pdz_pot = (cvec(3) - my_particle_list%particles(ipart)%z) / B_offload(irk)

                    !my_particle_list%particles(ipart)%xyz_abs(1) = my_particle_list%particles(ipart)%xyz_abs(1) + (cvec(1) - my_particle_list%particles(ipart)%x)
                    !my_particle_list%particles(ipart)%xyz_abs(2) = my_particle_list%particles(ipart)%xyz_abs(2) + (cvec(2) - my_particle_list%particles(ipart)%y)
                    !my_particle_list%particles(ipart)%xyz_abs(3) = my_particle_list%particles(ipart)%xyz_abs(3) + (cvec(3) - my_particle_list%particles(ipart)%z)

                    my_particle_list%particles(ipart)%x = cvec(1)
                    my_particle_list%particles(ipart)%y = cvec(2)
                    my_particle_list%particles(ipart)%z = cvec(3)

                    ! >>> end local "conceptual namespace" <<<

                    !CALL update_particle_cell_target(my_particle_list%particles(ipart), kk, jj, ii, x, y, z, dx, dy, dz)
                END IF
#endif
                ! >>>>>>>>>>>> UPDATE OF PARTICLE COORDINATES (neccesary for PER boundaries) AND GRID <<<<<<<<<<<<
               
                ! >>> begin local "conceptual namespace" <<<
                ! int_var_1 => destgrid

                CALL get_gcorner_neighbour(particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn), pstag, int_var_1)

                IF (int_var_1 == 0) int_var_1 = my_particle_list%particles(ipart)%igrid

!IF(int_var_1 < 1 .OR. int_var_1 > ngrid) WRITE(*,*) ">>>>>>>>>>>>>>>>>> ipart: ", my_particle_list%particles(ipart)%ipart, " igrid: ", &
!   my_particle_list%particles(ipart)%igrid, &
!  "destgrid: ", int_var_1, " ip: ", ip, "Pstag: ", pstag
                
                CALL update_coordinates_target(my_particle_list%particles(ipart), int_var_1, 99, bbox)
                
                my_particle_list%particles(ipart)%igrid = int_var_1

                ip = ip1d_offload(my_particle_list%particles(ipart)%igrid)

                IF (ip > 0) THEN
                    CALL get_mgdims_target(kk, jj, ii, my_particle_list%particles(ipart)%igrid)

                    CALL set_particle_cell_target(my_particle_list%particles(ipart), kk, jj, ii, &
                    x_offload(ip:ip+ii-1), y_offload(ip:ip+jj-1), z_offload(ip:ip+kk-1), &
                    dx_offload(ip:ip+ii-1), dy_offload(ip:ip+jj-1), dz_offload(ip:ip+kk-1))
                END IF

                !CALL update_coordinates_target3(particle_gcorner_boundaries((igrid - 1) * 8_intk + icorn), pstag, my_particle_list%particles(ipart))

                ! >>> end local "conceptual namespace" <<<

                ! TODO: reintroduce particle runtime statistics
            END DO
            !$omp end parallel do 
        END DO
        !$omp end teams distribute
        !$omp end target

        CALL stop_timer(924)

        CALL start_timer(925)
        
#if defined __INTEL_COMPILER
        !$omp target update from(my_particle_list%particles)
#else
        !$omp target update from(my_particle_list)
#endif

        CALL stop_timer(925)

        !CALL write_particle_list_txt(itstep, 'pos')

        WRITE(*, '("    Timeintegration on Process:                     ", I9)') myid
        WRITE(*, '("        Device Number:                              ", I9)') dev_num
        WRITE(*, '("        Number of Teams:                            ", I9)') num_teams
        WRITE(*, '("        Max. Number of Threds (per Team):           ", I9)') num_threads

        CALL stop_timer(920)

        CALL stop_timer(900)
    
    END SUBROUTINE timeintegrate_particles_target3


    SUBROUTINE finish_particle_timeintegration()

        !$omp target exit data map(delete: A_offload, B_offload)

        DEALLOCATE(A_offload)
        DEALLOCATE(B_offload)

    END SUBROUTINE finish_particle_timeintegration

END MODULE
