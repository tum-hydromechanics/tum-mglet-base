MODULE particle_timeintegration_mod

    USE particle_interpolation_mod
    USE particle_exchange_mod
    USE particle_list_mod
    USE particle_boundaries_mod

    IMPLICIT NONE

CONTAINS

    SUBROUTINE timeintegrate_particles(itstep, dt)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt

        ! local variables
        TYPE(baseparticle_t) :: particle_clone
        INTEGER(intk) :: igrid, i, ii, jj, kk, gfound, ig, destgrid
        REAL(realk) :: p_u, p_v, p_w

        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: xstag_f, ystag_f, zstag_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: u_f, v_f, w_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: xstag, ystag, zstag
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz, ddx, ddy, ddz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w

        REAL(realk) :: rand(3), diffusion_dx, diffusion_dy, diffusion_dz, advection_dx, advection_dy, advection_dz, pdx, pdy, pdz

        CALL start_timer(910)

        ! is calling the geometric fields obsolete when using core_mode -> corefields_mod ?
        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(xstag_f, "XSTAG")
        CALL get_field(ystag_f, "YSTAG")
        CALL get_field(zstag_f, "ZSTAG")

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
            WRITE(*, *) '' ! REMOVE later
            WRITE(*, '("NEW TIMESTEP - PARTICLE TIMEINTEGRATION:")')
        END IF

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

            ! Grid and Field Info
            CALL x_f%get_ptr(x, igrid)
            CALL y_f%get_ptr(y, igrid)
            CALL z_f%get_ptr(z, igrid)

            CALL xstag_f%get_ptr(xstag, igrid)
            CALL ystag_f%get_ptr(ystag, igrid)
            CALL zstag_f%get_ptr(zstag, igrid)

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

            END IF

            CALL get_mgdims(kk, jj, ii, igrid)

            ! Advection
            IF (dinterp_particles) THEN

                ! gobert
                CALL get_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                 p_u, p_v, p_w, u, v, w, x, y, z, dx, dy, dz, ddx, ddy, ddz)

            ELSE

                ! nearest cell
                CALL get_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                 p_u, p_v, p_w, u, v, w, x, y, z)

            END IF

            advection_dx = p_u * dt
            advection_dy = p_v * dt
            advection_dz = p_w * dt

            ! Diffusion
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(rand)

            rand(1) = rand(1) - 0.5_realk
            rand(2) = rand(2) - 0.5_realk
            rand(3) = rand(3) - 0.5_realk

            diffusion_dx = SQRT(2 * D(1) * dt) * rand(1) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))
            diffusion_dy = SQRT(2 * D(2) * dt) * rand(2) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))
            diffusion_dz = SQRT(2 * D(3) * dt) * rand(3) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))

            ! Advection + Diffusion
            pdx = advection_dx + diffusion_dx
            pdy = advection_dy + diffusion_dy
            pdz = advection_dz + diffusion_dz

            ! for debugging
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("Pre Particle Motion:")')
                    CALL print_particle_status(my_particle_list%particles(i))
                    WRITE(*,'("Unhindered Motion", I0, "(dt = ", F12.8,"):")') my_particle_list%particles(i)%ipart, dt
                    WRITE(*,'("Particle ADVECTION [m] = ", 3F12.8)') advection_dx, advection_dy, advection_dz
                    WRITE(*,'("Particle DIFFUSION [m] = ", 3F12.8)') diffusion_dx, diffusion_dy, diffusion_dz
                    WRITE(*, '()')
            END SELECT

            ! particle_clone = my_particle_list%particles(i)

            ! Particle Boundary Interaction
            CALL move_particle(my_particle_list%particles(i), pdx, pdy, pdz)

            ! my_particle_list%particles(i) = particle_clone

            ! count particles to be sent (per process)
            ! CALL prepare_particle_exchange(my_particle_list%particles(i))

        END DO

        ! migration of particles between grids
        ! and also across MPI ranks (processes)
        CALL exchange_particles(my_particle_list, itstep)

        CALL stop_timer(910)

    END SUBROUTINE timeintegrate_particles

END MODULE
