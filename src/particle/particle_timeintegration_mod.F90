! file for particle motion and timeintegration

MODULE particle_timeintegration_mod

    USE particle_list_mod
    USE particle_interpolation_mod
    USE particle_connect_mod
    !USE flow_mod (NOT NEEDED)

    IMPLICIT NONE

CONTAINS

    SUBROUTINE timeintegrate_particles(dt)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: dt

        ! local variables
        INTEGER(intk) :: igrid, nbrgrid, i, ii, jj, kk
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

        DO i = 1, my_particle_list%ifinal

            IF (.NOT. my_particle_list%particles(i)%is_active) THEN

                CYCLE

            END IF

            igrid = my_particle_list%particles(i)%igrid !just for clarity of the follwing expressions

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

            ! Connect
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            IF (my_particle_list%particles(i)%x < minx .OR. &
                my_particle_list%particles(i)%x > maxx .OR. &
                my_particle_list%particles(i)%y < miny .OR. &
                my_particle_list%particles(i)%y > maxy .OR. &
                my_particle_list%particles(i)%z < minz .OR. &
                my_particle_list%particles(i)%z > maxz) THEN

                CALL get_target_grid(my_particle_list%particles(i), nbrgrid)

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*,'("Particle ", I0, " left grid ", I0, " at ", 3F12.6)') my_particle_list%particles(i)%ipart, my_particle_list%particles(i)%igrid, my_particle_list%particles(i)%x, &
                        my_particle_list%particles(i)%y, my_particle_list%particles(i)%z
                END SELECT


                ! coordinate adaption if particle passes periodic boundary
                CALL get_bbox(nbrminx, nbrmaxx, nbrminy, nbrmaxy, nbrminz, nbrmaxz, nbrgrid)

                ! new grid for particle
                my_particle_list%particles(i)%igrid = nbrgrid

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*,'("New Grid: ", I0)') my_particle_list%particles(i)%igrid
                END SELECT

                ! PERIODIC BOUNDARY

                IF (my_particle_list%particles(i)%x < nbrminx) THEN
                    my_particle_list%particles(i)%x = nbrmaxx - ABS(my_particle_list%particles(i)%x - minx)
                END IF

                IF (nbrmaxx < my_particle_list%particles(i)%x) THEN
                    my_particle_list%particles(i)%x = nbrminx + ABS(my_particle_list%particles(i)%x - maxx)
                END IF

                IF (my_particle_list%particles(i)%y < nbrminy) THEN
                    my_particle_list%particles(i)%y = nbrmaxy - ABS(my_particle_list%particles(i)%y - miny)
                END IF

                IF (nbrmaxy < my_particle_list%particles(i)%y) THEN
                    my_particle_list%particles(i)%y = nbrminy + ABS(my_particle_list%particles(i)%y - maxy)
                END IF

                IF (my_particle_list%particles(i)%z < nbrminz) THEN
                    my_particle_list%particles(i)%z = nbrmaxz - ABS(my_particle_list%particles(i)%z - minz)
                END IF

                IF (nbrmaxz < my_particle_list%particles(i)%z) THEN
                    my_particle_list%particles(i)%z = nbrminz + ABS(my_particle_list%particles(i)%z - maxz)
                END IF

                ! NON PERIODIC BOUNDARY
                !IF (my_particle_list%particles(i)%x < nbrminx .OR. &
                !    my_particle_list%particles(i)%x > nbrmaxx .OR. &
                !    my_particle_list%particles(i)%y < nbrminy .OR. &
                !    my_particle_list%particles(i)%y > nbrmaxy .OR. &
                !    my_particle_list%particles(i)%z < nbrminz .OR. &
                !    my_particle_list%particles(i)%z > nbrmaxz) THEN

                    !SELECT CASE (TRIM(particle_terminal))
                    !    CASE ("none")
                    !        CONTINUE
                    !    CASE ("normal")
                    !        CONTINUE
                    !    CASE ("verbose")
                    !        WRITE(*,'("Particle ", I0, " crossed Periodic Boundary and was deactivated.")') my_particle_list%particles(i)%ipart
                    !END SELECT

                    !my_particle_list%particles(i)%is_active = .FALSE.

                    !my_particle_list%active_np = my_particle_list%active_np - 1_intk

                !END IF

                CALL my_particle_list%particles(i)%get_p_ijkcell()

            ELSE

                CALL my_particle_list%particles(i)%update_p_ijkcell(pdx, pdy, pdz)

            END IF

            CALL my_particle_list%particles(i)%print_status()

        END DO

    END SUBROUTINE timeintegrate_particles

END MODULE
