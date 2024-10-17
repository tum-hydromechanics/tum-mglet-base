MODULE particle_timeintegration_mod

    USE particle_interpolation_mod
    USE particle_exchange_mod
    USE particle_list_mod
    USE particle_boundaries_mod

    IMPLICIT NONE

    REAL(realk), ALLOCATABLE :: px_bup(:), py_bup(:), pz_bup(:)
    REAL(realk), ALLOCATABLE :: pu_itm(:,:), pv_itm(:,:), pw_itm(:,:)

    TYPE(particle_rk_t) :: prkscheme
    REAL(realk), ALLOCATABLE :: c(:), b(:), a(:,:)

CONTAINS

    SUBROUTINE init_particle_timeintegration()

        ! local variables
        INTEGER(intk), PARAMETER :: units_v(7) = [0, 1, -1, 0, 0, 0, 0]

        TYPE(field_t), POINTER :: u, v, w
        TYPE(field_t), POINTER :: buu, buv, buw

        CALL start_timer(900)
        CALL start_timer(910)

        ! set backup fields (U_BackUp  => UBU; V and W analogously), used to store flow field at the beginning of each timestep
        ! until the beginning of the next timestep; for particle timeintegration
        CALL set_field("UBU", istag=1, units=units_v, buffers=.TRUE.)
        CALL set_field("VBU", jstag=1, units=units_v, buffers=.TRUE.)
        CALL set_field("WBU", kstag=1, units=units_v, buffers=.TRUE.)

        ! set velocity field for Interpolation in rk steps (U Linear Interpolation => ULI; V and W analogously)
        CALL set_field("ULI", istag=1, units=units_v, buffers=.TRUE.)
        CALL set_field("VLI", jstag=1, units=units_v, buffers=.TRUE.)
        CALL set_field("WLI", kstag=1, units=units_v, buffers=.TRUE.)

        CALL prkscheme%init(prkmethod)

        ALLOCATE(px_bup(my_particle_list%max_np))
        ALLOCATE(py_bup(my_particle_list%max_np))
        ALLOCATE(pz_bup(my_particle_list%max_np))

        ALLOCATE(pu_itm(my_particle_list%max_np, prkscheme%nrk - 1))
        ALLOCATE(pv_itm(my_particle_list%max_np, prkscheme%nrk - 1))
        ALLOCATE(pw_itm(my_particle_list%max_np, prkscheme%nrk - 1))

        ! init rk scheme
        CALL prkscheme%get_coeffs(c, b, a)

        ! for debugging
        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                WRITE(*,'("RK-Coefficients:")')
                WRITE(*,'("c(1) = ", F12.8)') c(1)
                WRITE(*,'("c(2) = ", F12.8)') c(2)
                WRITE(*,'("c(3) = ", F12.8)') c(3)
                WRITE(*,'("b(1) = ", F12.8)') b(1)
                WRITE(*,'("b(2) = ", F12.8)') b(2)
                WRITE(*,'("b(3) = ", F12.8)') b(3)
                WRITE(*,'("a(2,1) = ", F12.8)') a(2,1)
                WRITE(*,'("a(3,1) = ", F12.8)') a(3,1)
                WRITE(*,'("a(3,2) = ", F12.8)') a(3,2)
                WRITE(*, '()')
        END SELECT

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_timeintegration

    SUBROUTINE prepare_particle_timeintegration()

        CALL update_backup_fields()

        IF (SIZE(px_bup) /= my_particle_list%max_np) THEN

            DEALLOCATE(px_bup)
            DEALLOCATE(py_bup)
            DEALLOCATE(pz_bup)
            DEALLOCATE(pu_itm)
            DEALLOCATE(pv_itm)
            DEALLOCATE(pw_itm)

            ALLOCATE(px_bup(my_particle_list%max_np))
            ALLOCATE(py_bup(my_particle_list%max_np))
            ALLOCATE(pz_bup(my_particle_list%max_np))
            ALLOCATE(pu_itm(my_particle_list%max_np, prkscheme%nrk - 1))
            ALLOCATE(pv_itm(my_particle_list%max_np, prkscheme%nrk - 1))
            ALLOCATE(pw_itm(my_particle_list%max_np, prkscheme%nrk - 1))

        END IF

    END SUBROUTINE prepare_particle_timeintegration

    SUBROUTINE timeintegrate_particles(itstep, dt)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt

        ! local variables
        !TYPE(baseparticle_t) :: particle_clone
        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: u_f, v_f, w_f
        TYPE(field_t), POINTER :: ubu_f, vbu_f, wbu_f
        TYPE(field_t), POINTER :: uli_f, vli_f, wli_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uli, vli, wli

        INTEGER(intk) :: igrid, i, j, ii, jj, kk, gfound, ig, destgrid
        REAL(realk) :: pu_adv, pv_adv, pw_adv, pu_diff, pv_diff, pw_diff
        REAL(realk) :: pdx, pdy, pdz, pdx_diff, pdy_diff, pdz_diff, pdx_eff, pdy_eff, pdz_eff

        INTEGER(intk) :: irk

        CALL start_timer(900)
        CALL start_timer(920)

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        CALL get_field(ubu_f, "UBU")
        CALL get_field(vbu_f, "VBU")
        CALL get_field(wbu_f, "WBU")

        CALL get_field(uli_f, "ULI")
        CALL get_field(vli_f, "VLI")
        CALL get_field(wli_f, "WLI")

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

        ! store the initial particle coordinates in backup array
        DO i = 1, my_particle_list%ifinal

            px_bup(i) = my_particle_list%particles(i)%x
            py_bup(i) = my_particle_list%particles(i)%y
            pz_bup(i) = my_particle_list%particles(i)%z

        END DO

        ! algorithm for EXPLICIT RK schemes
        ! iterate over all particles withing the rk iteration so the field interpolation can be used für all particles
        DO irk = 1, prkscheme%nrk

            DO i = 1, my_particle_list%ifinal

                CALL time_interpolate_field(c(irk), u_f, v_f, w_f, ubu_f, vbu_f, wbu_f, uli_f, vli_f, wli_f)

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

                ! --- ADVECTION VELOCITY ---

                CALL get_mgdims(kk, jj, ii, igrid)

                ! Grid and Field Info
                CALL x_f%get_ptr(x, igrid)
                CALL y_f%get_ptr(y, igrid)
                CALL z_f%get_ptr(z, igrid)

                CALL uli_f%get_ptr(uli, igrid)
                CALL vli_f%get_ptr(vli, igrid)
                CALL wli_f%get_ptr(wli, igrid)

                IF (dinterp_particles) THEN

                    CALL dx_f%get_ptr(dx, igrid)
                    CALL dy_f%get_ptr(dy, igrid)
                    CALL dz_f%get_ptr(dz, igrid)

                    CALL ddx_f%get_ptr(ddx, igrid)
                    CALL ddy_f%get_ptr(ddy, igrid)
                    CALL ddz_f%get_ptr(ddz, igrid)

                    CALL gobert_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                    pu_adv, pv_adv, pw_adv, uli, vli, wli, x, y, z, dx, dy, dz, ddx, ddy, ddz)

                ELSE

                    CALL nearest_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                    pu_adv, pv_adv, pw_adv, uli, vli, wli, x, y, z)

                END IF

                ! store intermediate velocity
                pu_itm(i, irk) = pu_adv
                pv_itm(i, irk) = pv_adv
                pw_itm(i, irk) = pw_adv

                IF (irk /= prkscheme%nrk) THEN
                    ! zeroize current intermediate particle displacement
                    pdx = 0
                    pdy = 0
                    pdz = 0

                    ! Advection
                    DO j = 1, prkscheme%nrk
                        pdx = pdx + dt * pu_itm(i, j) * a(irk + 1, j)
                        pdy = pdy + dt * pv_itm(i, j) * a(irk + 1, j)
                        pdz = pdz + dt * pw_itm(i, j) * a(irk + 1, j)
                    END DO

                    ! Particle Boundary Interaction
                    CALL move_particle(my_particle_list%particles(i), pdx, pdy, pdz, pdx_eff, pdy_eff, pdz_eff)

                ELSEIF (irk == prkscheme%nrk) THEN

                    ! set particle coordinates to backup value
                    my_particle_list%particles(i)%x = px_bup(i)
                    my_particle_list%particles(i)%y = py_bup(i)
                    my_particle_list%particles(i)%z = pz_bup(i)

                    CALL update_particle_cell(my_particle_list%particles(i))

                    ! --- Diffusiion ---
                    CALL get_particle_diffusion(dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff)

                    ! --- Advection ---
                    pdx = dt * b(irk) * pu_adv
                    pdy = dt * b(irk) * pv_adv
                    pdz = dt * b(irk) * pw_adv

                    DO j = 1, prkscheme%nrk - 1
                        pdx = pdx + dt * b(j) * pu_itm(i, j)
                        pdy = pdy + dt * b(j) * pv_itm(i, j)
                        pdz = pdz + dt * b(j) * pw_itm(i, j)
                    END DO

                    ! for debugging
                    SELECT CASE (TRIM(particle_terminal))
                        CASE ("none")
                            CONTINUE
                        CASE ("normal")
                            CONTINUE
                        CASE ("verbose")
                            WRITE(*,'("Pre Physical Particle Motion:")')
                            CALL print_particle_status(my_particle_list%particles(i))
                            WRITE(*,'("Unhindered Motion", I0, "(dt = ", F12.8,"):")') my_particle_list%particles(i)%ipart, dt
                            WRITE(*,'("Particle ADVECTION [m] = ", 3F12.8)') pdx, pdy, pdz
                            WRITE(*,'("Particle DIFFUSION [m] = ", 3F12.8)') pdx_diff, pdy_diff, pdz_diff
                            WRITE(*, '()')
                    END SELECT

                    ! --- Superposition ---
                    pdx = pdx + pdx_diff
                    pdy = pdy + pdy_diff
                    pdz = pdz + pdz_diff

                    ! THIS IS THE ACTUAL "PHYSICAL" DISPLACEMENT OF THE PARTICLE
                    CALL move_particle(my_particle_list%particles(i), pdx, pdy, pdz, pdx_eff, pdy_eff, pdz_eff)

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

    SUBROUTINE update_backup_fields()

        ! local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f
        TYPE(field_t), POINTER :: ubu_f, vbu_f, wbu_f

        CALL start_timer(900)
        CALL start_timer(920)
        CALL start_timer(921)

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        CALL get_field(ubu_f, "UBU")
        CALL get_field(vbu_f, "VBU")
        CALL get_field(wbu_f, "WBU")

        ! Copy the velocity field u, v, w into ubu, vbu, wbu
        ubu_f%arr = u_f%arr
        vbu_f%arr = v_f%arr
        wbu_f%arr = w_f%arr

        ! Copy buffers from u, v, w to ubu, vbu, wbu
        ubu_f%buffers = u_f%buffers
        vbu_f%buffers = v_f%buffers
        wbu_f%buffers = w_f%buffers

        CALL stop_timer(921)
        CALL stop_timer(920)
        CALL stop_timer(900)

    END SUBROUTINE update_backup_fields

    SUBROUTINE finish_particle_timeintegration()

        DEALLOCATE(px_bup)
        DEALLOCATE(py_bup)
        DEALLOCATE(pz_bup)
        DEALLOCATE(pu_itm)
        DEALLOCATE(pv_itm)
        DEALLOCATE(pw_itm)

    END SUBROUTINE finish_particle_timeintegration

END MODULE
