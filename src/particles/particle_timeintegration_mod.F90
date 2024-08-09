! file for particle motion and timeintegration

MODULE particle_timeintegration_mod

    USE particle_list_mod
    USE particle_interpolation_mod
    USE particle_exchange_mod
    USE core_mod

    IMPLICIT NONE

CONTAINS

    SUBROUTINE timeintegrate_particles(dt)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: dt

        ! local variables
        INTEGER(intk) :: igrid, nbrgrid, nbrproc, i, ii, jj, kk, gfound, ig
        REAL(realk) :: p_u, p_v, p_w
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, nbrminx, nbrmaxx, nbrminy, nbrmaxy, nbrminz, nbrmaxz

        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: xstag_f, ystag_f, zstag_f
        TYPE(field_t), POINTER :: u_f, v_f, w_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: xstag, ystag, zstag
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w

        REAL(realk) :: rand(3), diffusion_dx, diffusion_dy, diffusion_dz, advection_dx, advection_dy, advection_dz, pdx, pdy, pdz

        !is calling the geometric fields obsolete when using core_mode -> corefields_mod ?
        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(xstag_f, "XSTAG")
        CALL get_field(ystag_f, "YSTAG")
        CALL get_field(zstag_f, "ZSTAG")

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        IF (myid == 0) THEN
            WRITE(*, *) '' ! REMOVE later
            WRITE(*, *) 'New Timestep ...' ! REMOVE later
        END IF

        DO i = 1, my_particle_list%ifinal

            ! checking activity
            IF ( my_particle_list%particles(i)%is_active /= 1 ) THEN
                CYCLE
            END IF

            ! checking locality (Debug)
            IF ( my_particle_list%particles(i)%iproc /= myid ) THEN
                WRITE(*,*) 'Particle on wrong proc at start of timestep'
                CALL errr(__FILE__, __LINE__)
            END IF

            ! getting the grid
            igrid = my_particle_list%particles(i)%igrid !just for clarity of the follwing expressions

            ! checking consistency (Debug)
            gfound = 1
            DO ig = 1, nMyGrids
                IF ( my_particle_list%particles(i)%igrid == mygrids(ig) ) gfound = 1; EXIT
            END DO

            IF ( gfound == 0 ) THEN
                WRITE(*, '("Particle grid ", I0, " is not on this process (", I0, ")")') my_particle_list%particles(i)%igrid, myid
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

            CALL get_mgdims(kk, jj, ii, igrid)

            ! EXPLICIT EULER FOR NOW

            ! Advection
            CALL get_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                 p_u, p_v, p_w, xstag, ystag, zstag, u, v, w, x, y, z)

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

            pdx = advection_dx + diffusion_dx
            pdy = advection_dy + diffusion_dy
            pdz = advection_dz + diffusion_dz

            ! Advection + DIffusion
            my_particle_list%particles(i)%x = my_particle_list%particles(i)%x + pdx
            my_particle_list%particles(i)%y = my_particle_list%particles(i)%y + pdy
            my_particle_list%particles(i)%z = my_particle_list%particles(i)%z + pdz

            ! for debugging
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("Timeintegrating Particle ", I0, "(dt = ", F12.8,"):")') my_particle_list%particles(i)%ipart, dt
                    WRITE(*,'("Particle ADVECTION [m] = ", 3F12.8)') advection_dx, advection_dy, advection_dz
                    WRITE(*,'("Particle DIFFUSION [m] = ", 3F12.8)') diffusion_dx, diffusion_dy, diffusion_dz
                    WRITE(*, *) ' '
            END SELECT

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("Pre Migration:")')
                    CALL print_particle_status(my_particle_list%particles(i))
            END SELECT

        END DO

        ! migration of particles between grids
        ! and also across MPI ranks (processes)
        CALL exchange_particles(my_particle_list)

    END SUBROUTINE timeintegrate_particles

END MODULE
