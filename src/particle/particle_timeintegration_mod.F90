! file for particle motion and timeintegration

MODULE particle_timeintegration_mod

    USE particle_list_mod
    !USE flow_mod

    IMPLICIT NONE

CONTAINS

    SUBROUTINE timeintegrate_particles(dt)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: dt

        ! local variables
        INTEGER(intk) :: igrid, i, ii, jj, kk
        REAL(realk) :: p_u, p_v, p_w
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: xstag_f, ystag_f, zstag_f
        TYPE(field_t), POINTER :: u_f, v_f, w_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: xstag, ystag, zstag
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w

        REAL(realk) :: rand(3), diffusion_dx, diffusion_dy, diffusion_dz

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

            igrid = my_particle_list%particles(i)%igrid !just for clarity of the coming expressions

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

            IF (dinterp_particles) THEN

                CALL interp_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                 p_u, p_v, p_w, xstag, ystag, zstag, u, v, w, x, y, z)

            ELSE

                CALL get_particle_uvw(kk, jj, ii, my_particle_list%particles(i), &
                 p_u, p_v, p_w, u, v, w, x, y, z)

            END IF

            ! TODO: proper timeintegration
            ! dummy method for now (explicit euler)

            ! how should 2 dimensional cases be handeled ?
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(rand)

            rand(1) = rand(1) - 0.5_realk
            rand(2) = rand(2) - 0.5_realk
            rand(3) = rand(3) - 0.5_realk

            diffusion_dx = SQRT(2 * D * dt) * rand(1) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))
            diffusion_dy = SQRT(2 * D * dt) * rand(2) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))
            diffusion_dz = SQRT(2 * D * dt) * rand(3) / SQRT(rand(1)**(2) + rand(2)**(2) + rand(3)**(2))

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("Particle ", I0, " moved from ", 3F12.6)', advance = "no") my_particle_list%particles(i)%ipart, my_particle_list%particles(i)%x, &
                     my_particle_list%particles(i)%y, my_particle_list%particles(i)%z
            END SELECT

            my_particle_list%particles(i)%x = my_particle_list%particles(i)%x + diffusion_dx + p_u * dt
            my_particle_list%particles(i)%y = my_particle_list%particles(i)%y + diffusion_dy + p_v * dt
            my_particle_list%particles(i)%z = my_particle_list%particles(i)%z + diffusion_dz + p_w * dt

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("   to ", 3F12.6)') my_particle_list%particles(i)%x, &
                     my_particle_list%particles(i)%y, my_particle_list%particles(i)%z
                    WRITE(*,'("Particle velocity = ", 3F12.6, " | dt = ", F12.6)') p_u, p_v, p_w, dt
                    WRITE(*, *) ' '
            END SELECT

            ! deactivate particles that leave the domain (only for simple tests with one grid for now) ...

            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            IF (my_particle_list%particles(i)%x < minx .OR. &
                my_particle_list%particles(i)%x > maxx .OR. &
                my_particle_list%particles(i)%y < miny .OR. &
                my_particle_list%particles(i)%y > maxy .OR. &
                my_particle_list%particles(i)%z < minz .OR. &
                my_particle_list%particles(i)%z > maxz) THEN

                my_particle_list%particles(i)%is_active = .FALSE.

                my_particle_list%active_np = my_particle_list%active_np - 1_intk

                SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*,'("Particle ", I0, " left domian at ", 3F12.6)') my_particle_list%particles(i)%ipart, my_particle_list%particles(i)%x, &
                     my_particle_list%particles(i)%y, my_particle_list%particles(i)%z
                END SELECT

            END IF

        END DO

    END SUBROUTINE timeintegrate_particles

    !------------------------------

    SUBROUTINE get_particle_uvw(kk, jj, ii, particle, p_u, p_v, p_w, &
         u, v, w, x, y, z)
    ! no interpolation of velocity, just the velocity components of the nearest staggered cells respectively

                !subroutine_arguments
                INTEGER(intk), INTENT(in) :: kk, jj, ii
                CLASS(baseparticle_t), INTENT(in) :: particle
                REAL(realk), INTENT(out) :: p_u, p_v, p_w
                REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
                REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)

                !local variables
                INTEGER(intk) :: p_istag, p_jstag, p_kstag

                p_istag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%x - x(particle%ijkcell(1)) ), intk ))
                p_jstag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%y - x(particle%ijkcell(2)) ), intk ))
                p_kstag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%z - x(particle%ijkcell(3)) ), intk ))

                p_u = u(particle%ijkcell(3), particle%ijkcell(2), particle%ijkcell(1) + p_istag - 1)
                p_v = v(particle%ijkcell(3), particle%ijkcell(2) + p_jstag - 1, particle%ijkcell(1))
                p_w = w(particle%ijkcell(3) + p_kstag - 1, particle%ijkcell(2), particle%ijkcell(1))

    END SUBROUTINE get_particle_uvw

    !------------------------------

    SUBROUTINE interp_particle_uvw(kk, jj, ii, particle, p_u, p_v, p_w, &
         xstag, ystag, zstag, u, v, w, x, y, z)

            ! subroutine arguments
            INTEGER(intk), INTENT(in) :: kk, jj, ii
            CLASS(baseparticle_t), INTENT(in) :: particle
            REAL(realk), INTENT(out) :: p_u, p_v, p_w
            REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
            REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
            REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)

            ! local variables
            INTEGER(intk) :: p_i, p_j, p_k, p_istag, p_jstag, p_kstag
            REAL(realk) :: p_x, p_y, p_z
            REAL(realk), DIMENSION(6) :: p_bbox_u, p_bbox_v, p_bbox_w

            !just for readability of the coming expressions
            p_i = particle%ijkcell(1)
            p_j = particle%ijkcell(2)
            p_k = particle%ijkcell(3)

            !just for readability
            p_x = particle%x
            p_y = particle%y
            p_z = particle%z

            p_istag = MAX( 0_intk, NINT( SIGN(1.0_realk, p_x - x(p_i)), intk ) )
            p_jstag = MAX( 0_intk, NINT( SIGN(1.0_realk, p_y - y(p_j)), intk ) )
            p_kstag = MAX( 0_intk, NINT( SIGN(1.0_realk, p_z - z(p_k)), intk ) )

            p_bbox_u(1) = xstag(p_i - 1)
            p_bbox_u(2) = xstag(p_i)
            p_bbox_u(3) = y(p_j + p_jstag - 1)
            p_bbox_u(4) = y(p_j + p_jstag)
            p_bbox_u(5) = z(p_k + p_kstag - 1)
            p_bbox_u(6) = z(p_k + p_kstag)

            p_bbox_v(1) = x(p_i + p_istag - 1)
            p_bbox_v(2) = x(p_i + p_istag)
            p_bbox_v(3) = ystag(p_j - 1)
            p_bbox_v(4) = ystag(p_j)
            p_bbox_v(5) = z(p_k + p_kstag - 1)
            p_bbox_v(6) = z(p_k + p_kstag)

            p_bbox_w(1) = x(p_i + p_istag - 1)
            p_bbox_w(2) = x(p_i + p_istag)
            p_bbox_w(3) = y(p_j + p_jstag - 1)
            p_bbox_w(4) = y(p_j + p_jstag)
            p_bbox_w(5) = zstag(p_k - 1)
            p_bbox_w(6) = zstag(p_k)

            CALL interp_trilinear(p_x, p_y, p_z, p_bbox_u, p_u, &
                            u(p_k + p_kstag - 1, p_j + p_jstag -1, p_i - 1),&
                            u(p_k + p_kstag    , p_j + p_jstag -1, p_i - 1),&
                            u(p_k + p_kstag - 1, p_j + p_jstag   , p_i - 1),&
                            u(p_k + p_kstag    , p_j + p_jstag   , p_i - 1),&
                            u(p_k + p_kstag - 1, p_j + p_jstag -1, p_i    ),&
                            u(p_k + p_kstag    , p_j + p_jstag -1, p_i    ),&
                            u(p_k + p_kstag - 1, p_j + p_jstag   , p_i    ),&
                            u(p_k + p_kstag    , p_j + p_jstag   , p_i    ))

            CALL interp_trilinear(p_x, p_y, p_z, p_bbox_v, p_v,&
                            v(p_k + p_kstag - 1, p_j -1, p_i + p_istag - 1),&
                            v(p_k + p_kstag    , p_j -1, p_i + p_istag - 1),&
                            v(p_k + p_kstag - 1, p_j   , p_i + p_istag - 1),&
                            v(p_k + p_kstag    , p_j   , p_i + p_istag - 1),&
                            v(p_k + p_kstag - 1, p_j -1, p_i + p_istag    ),&
                            v(p_k + p_kstag    , p_j -1, p_i + p_istag    ),&
                            v(p_k + p_kstag - 1, p_j   , p_i + p_istag    ),&
                            v(p_k + p_kstag    , p_j   , p_i + p_istag    ))

            CALL interp_trilinear(p_x, p_y, p_z, p_bbox_w, p_w,&
                            w(p_k - 1, p_j + p_jstag -1, p_i + p_istag - 1),&
                            w(p_k    , p_j + p_jstag -1, p_i + p_istag - 1),&
                            w(p_k - 1, p_j + p_jstag   , p_i + p_istag - 1),&
                            w(p_k    , p_j + p_jstag   , p_i + p_istag - 1),&
                            w(p_k - 1, p_j + p_jstag -1, p_i + p_istag    ),&
                            w(p_k    , p_j + p_jstag -1, p_i + p_istag    ),&
                            w(p_k - 1, p_j + p_jstag   , p_i + p_istag    ),&
                            w(p_k    , p_j + p_jstag   , p_i + p_istag    ))

    END SUBROUTINE interp_particle_uvw

    !------------------------------

    SUBROUTINE interp_trilinear(x, y, z, bbox, &
        f, f19, f20, f21, f22, f23, f24, f25, f26) ! Trilinear interpolation for a rectilinear grid !

        ! subroutine arguments
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), DIMENSION(6), INTENT(in) :: bbox
        REAL(realk), INTENT(in) :: f19, f20, f21, f22, f23, f24, f25, f26 ! numbers indicate the corners (see convention on the numbering of faces/edges/corners)
        REAL(realk), INTENT(out) :: f

        ! local variables
        REAL(realk) :: xd, yd, zd

        xd = (x - bbox(1)) / (bbox(2) - bbox(1))
        yd = (y - bbox(3)) / (bbox(4) - bbox(3))
        zd = (z - bbox(5)) / (bbox(6) - bbox(5))

        f = (f19 * (1 - xd) + f23 * xd)  * (1 - yd) + (f21 * (1 - xd) + f25 * xd) * yd * (1 - zd) + &
            (f20 * (1 - xd) + f24 * xd) * (1 - yd) + (f22 * (1 - xd) + f26 * xd) * yd * zd

    END SUBROUTINE interp_trilinear

END MODULE
