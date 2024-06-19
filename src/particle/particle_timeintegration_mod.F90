! file for particle motion and timeintegration

MODULE particle_timeintegration_mod

    USE particle_list_mod
    USE flow_mod

    IMPLICIT NONE

    REAL(realk), PARAMETER :: D = 0.1 ! Diffusion constant in m²/s (homogeneous and isotropic diffusion for now)

CONTAINS

    SUBROUTINE timeintegrate_particles(dt)

        ! subroutine arguments
        REAL(realk) :: dt

        ! local variables
        INTEGER(intk) :: igrid, i, ii, jj, kk
        REAL(realk) :: p_u, p_v, p_w
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        TYPE(field_t), POINTER :: x_f, y_f, z_f
        TYPE(field_t), POINTER :: xstag_f, ystag_f, zstag_f
        TYPE(field_t), POINTER :: u_f, v_f, w_f

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: x, y, z
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: xstag, ystag, zstag
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w

        REAL(realk) :: rand_dx, rand_dy, rand_dz, diffusion_dx, diffusion_dy, diffusion_dz

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

            IF (.NOT. my_particle_list%particles(i)%is_init) THEN

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

            CALL interp_particle_uvw(kk, jj, ii, &
                my_particle_list%particles(i), &
                p_u, p_v, p_w, xstag, ystag, zstag, u, v, w, x, y, z)

            ! TODO: proper timeintegration
            ! dummy method for now (explicit euler)

            ! how should 2 dimensional cases be handeled ?
            CALL RANDOM_NUMBER(rand_dx)
            CALL RANDOM_NUMBER(rand_dy)
            CALL RANDOM_NUMBER(rand_dz)

            diffusion_dx = SQRT(2 * D * dt) * rand_dx / SQRT(rand_dx**(2) + rand_dy**(2) + rand_dz**(2))
            diffusion_dy = SQRT(2 * D * dt) * rand_dy / SQRT(rand_dx**(2) + rand_dy**(2) + rand_dz**(2))
            diffusion_dz = SQRT(2 * D * dt) * rand_dz / SQRT(rand_dx**(2) + rand_dy**(2) + rand_dz**(2))

            my_particle_list%particles(i)%x = my_particle_list%particles(i)%x + p_u * dt + diffusion_dx
            my_particle_list%particles(i)%y = my_particle_list%particles(i)%y + p_v * dt + diffusion_dy
            my_particle_list%particles(i)%z = my_particle_list%particles(i)%z + p_w * dt + diffusion_dz

            ! deactivate particles that leave the domain (only for simple tests with one grid for now) ...

            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            IF (my_particle_list%particles(i)%x < minx .OR. &
                my_particle_list%particles(i)%x > maxx .OR. &
                my_particle_list%particles(i)%y < miny .OR. &
                my_particle_list%particles(i)%y > maxy .OR. &
                my_particle_list%particles(i)%z < minz .OR. &
                my_particle_list%particles(i)%z > maxz) THEN

                my_particle_list%particles(i)%is_init = .FALSE.

                my_particle_list%particle_stored(i) = .FALSE.

                WRITE(*,'("Particle ", I0, " left domian!")') my_particle_list%particles(i)%ipart
                WRITE(*,'("Final Particle Coordinates: ", 3F6.2)') my_particle_list%particles(i)%x, &
                 my_particle_list%particles(i)%y, my_particle_list%particles(i)%z

            END IF

        END DO

    END SUBROUTINE timeintegrate_particles

    !------------------------------

    SUBROUTINE interp_particle_uvw(kk, jj, ii, particle, p_u, p_v, p_w, &
        xstag, ystag, zstag, u, v, w, x, y, z)

            ! subroutine arguments
            INTEGER(intk), INTENT(in) :: kk, jj, ii
            CLASS(baseparticle_t), INTENT(in) :: particle
            REAL(realk), INTENT(out) :: p_u, p_v, p_w
            REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
            REAL(realk), INTENT(in) :: u(kk,jj,ii), v(kk,jj,ii), w(kk,jj,ii)
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

            p_istag = NINT( SIGN(1.0_realk, p_x - x(p_i)), intk )
            p_jstag = NINT( SIGN(1.0_realk, p_y - x(p_j)), intk )
            p_kstag = NINT( SIGN(1.0_realk, p_z - x(p_k)), intk )

            p_bbox_u(1) = xstag(p_i - 1)
            p_bbox_u(2) = xstag(p_i)
            p_bbox_u(3) = ystag(p_j + p_jstag - 1)
            p_bbox_u(4) = ystag(p_j + p_jstag)
            p_bbox_u(5) = zstag(p_k + p_kstag - 1)
            p_bbox_u(6) = zstag(p_k + p_kstag)

            p_bbox_v(1) = xstag(p_i + p_istag - 1)
            p_bbox_v(2) = xstag(p_i + p_istag)
            p_bbox_v(3) = ystag(p_j - 1)
            p_bbox_v(4) = ystag(p_j)
            p_bbox_v(5) = zstag(p_k + p_kstag - 1)
            p_bbox_v(6) = zstag(p_k + p_kstag)

            p_bbox_w(1) = xstag(p_i + p_istag - 1)
            p_bbox_w(2) = xstag(p_i + p_istag)
            p_bbox_w(3) = ystag(p_j + p_jstag - 1)
            p_bbox_w(4) = ystag(p_j + p_jstag)
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
        REAL(realk), INTENT(in) :: f19, f20, f21, f22, f23, f24, f25, f26 ! number indicate the corners (see convention on the numbering of faces/edges/corners)
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
