MODULE particle_interpolation_mod

    ! This module is responsible for:
    ! Interpolation of the flow field to deduce the advective velocity of particles.

    USE particle_core_mod

    IMPLICIT NONE

CONTAINS    !===================================

    SUBROUTINE get_particle_uvw(kk, jj, ii, particle, &
                 p_u, p_v, p_w, u, v, w, x, y, z, dx, dy, dz, ddx, ddy, ddz)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        TYPE(baseparticle_t), INTENT(in) :: particle
        REAL(realk), INTENT(out) :: p_u, p_v, p_w
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in), OPTIONAL :: dx(ii), dy(jj), dz(kk), ddx(ii), ddy(jj), ddz(kk)

        CALL start_timer(923)

        IF (dinterp_particles .AND. PRESENT(dx) .AND. PRESENT(dy) .AND. PRESENT(dz) &
         .AND. PRESENT(ddx) .AND. PRESENT(ddy) .AND. PRESENT(ddz)) THEN

            CALL gobert_particle_uvw(kk, jj, ii, particle, &
             p_u, p_v, p_w, u, v, w, x, y, z, dx, dy, dz, ddx, ddy, ddz)

        ELSE

            CALL nearest_particle_uvw(kk, jj, ii, particle, &
             p_u, p_v, p_w, u, v, w, x, y, z)

        END IF

        CALL stop_timer(923)

    END SUBROUTINE get_particle_uvw

    !------------------------------

    ! no interpolation of velocity, just the velocity components of the nearest staggered cells respectively
    SUBROUTINE nearest_particle_uvw(kk, jj, ii, particle, p_u, p_v, p_w, &
         u, v, w, x, y, z)

        !subroutine_arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        TYPE(baseparticle_t), INTENT(in) :: particle
        REAL(realk), INTENT(out) :: p_u, p_v, p_w
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)

        !local variables
        INTEGER(intk) :: p_istag, p_jstag, p_kstag

        CALL start_timer(923)

        p_istag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%x - x(particle%ijkcell(1)) ), intk ))
        p_jstag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%y - x(particle%ijkcell(2)) ), intk ))
        p_kstag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%z - x(particle%ijkcell(3)) ), intk ))

        p_u = u(particle%ijkcell(3), particle%ijkcell(2), particle%ijkcell(1) + p_istag - 1)
        p_v = v(particle%ijkcell(3), particle%ijkcell(2) + p_jstag - 1, particle%ijkcell(1))
        p_w = w(particle%ijkcell(3) + p_kstag - 1, particle%ijkcell(2), particle%ijkcell(1))

        CALL stop_timer(923)

    END SUBROUTINE nearest_particle_uvw

    !------------------------------

    ! from Gobert et. al, LAGRANGIAN SCALAR TRACKING FOR LAMINAR MICROMIXING AT HIGH SCHMIDT NUMBERS, 2006
    SUBROUTINE gobert_particle_uvw(kk, jj, ii, particle, p_u, p_v, p_w, &
         u, v, w, x, y, z, dx, dy, dz, ddx, ddy, ddz)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        TYPE(baseparticle_t), INTENT(in) :: particle
        REAL(realk), INTENT(out) :: p_u, p_v, p_w
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk), dx(ii), dy(jj), dz(kk), ddx(ii), ddy(jj), ddz(kk)

        ! local variables
        INTEGER(intk) :: p_i, p_j, p_k
        REAL(realk) :: p_x, p_y, p_z, alpha, beta, gamma, delta

        CALL start_timer(923)

        !just for readability of the following expressions
        p_i = particle%ijkcell(1)
        p_j = particle%ijkcell(2)
        p_k = particle%ijkcell(3)

        !just for readability of the following expressions
        p_x = particle%x
        p_y = particle%y
        p_z = particle%z

        ! u interpolation
        alpha = (u(p_k, p_j, p_i) - u(p_k, p_j, p_i - 1)) / ddx(p_i)

        beta = 0.25 * ((u(p_k, p_j + 1, p_i) + u(p_k, p_j + 1, p_i - 1) - u(p_k, p_j, p_i) - u(p_k, p_j, p_i - 1)) / dy(p_j) &
         + (u(p_k, p_j, p_i) + u(p_k, p_j, p_i - 1) - u(p_k, p_j - 1, p_i) - u(p_k, p_j - 1, p_i - 1)) / dy(p_j -1))

        gamma = 0.25 * ((u(p_k + 1, p_j, p_i) + u(p_k + 1, p_j, p_i - 1) - u(p_k, p_j, p_i) - u(p_k, p_j, p_i - 1)) / dz(p_k) &
         + (u(p_k, p_j, p_i) + u(p_k, p_j, p_i - 1) - u(p_k - 1, p_j, p_i) - u(p_k - 1, p_j, p_i - 1)) / dz(p_k - 1))

        delta = 0.5 * (u(p_k, p_j, p_i) + u(p_k, p_j, p_i - 1) &
         - alpha * (ddx(p_i) - dx(p_i - 1)) - beta * (ddy(p_j) - dy(p_j - 1)) - gamma * (ddz(p_k) - dz(p_k - 1)))

        p_u = alpha * (p_x - x(p_i)) + beta * (p_y - y(p_j)) + gamma * (p_z - z(p_k)) + delta

        ! v interpolation
        alpha = (v(p_k, p_j, p_i) - v(p_k, p_j - 1, p_i)) / ddy(p_j)

        beta = 0.25 * ((v(p_k, p_j, p_i + 1) + v(p_k, p_j - 1, p_i + 1) - v(p_k, p_j, p_i) - v(p_k, p_j - 1, p_i)) / dx(p_i) &
         + (v(p_k, p_j, p_i) + v(p_k, p_j - 1, p_i) - v(p_k, p_j, p_i - 1) - v(p_k, p_j - 1, p_i - 1)) /dx(p_i - 1))

        gamma = 0.25 * ((v(p_k + 1, p_j, p_i) + v(p_k + 1, p_j - 1, p_i) - v(p_k, p_j, p_i) - v(p_k, p_j - 1, p_i)) / dz(p_k) &
         + (v(p_k, p_j, p_i) + v(p_k, p_j - 1, p_i) - v(p_k - 1, p_j, p_i) - v(p_k - 1, p_j - 1, p_i)) / dz(p_k - 1))

        delta = 0.5 * (v(p_k, p_j, p_i) + v(p_k, p_j - 1, p_i) &
         - alpha * (ddy(p_j) - dy(p_j - 1)) - beta * (ddx(p_i) - dx(p_i - 1)) - gamma * (ddz(p_k) - dz(p_k - 1)))

        p_v = alpha * (p_y - y(p_j)) + beta * (p_x - x(p_i)) + gamma * (p_z - z(p_k)) + delta

        ! w interpolation
        alpha = (w(p_k, p_j, p_i) - w(p_k - 1 ,p_j ,p_i)) / ddz(p_k)

        beta = 0.25 * ((w(p_k, p_j, p_i + 1) + w(p_k - 1, p_j, p_i + 1) - w(p_k, p_j, p_i) - w(p_k - 1, p_j, p_i)) / dx(p_i) &
         + (w(p_k, p_j, p_i) + w(p_k - 1, p_j, p_i) - w(p_k, p_j, p_i - 1) - w(p_k - 1, p_j, p_i - 1)) / dx(p_i - 1))

        gamma = 0.25 * ((w(p_k, p_j + 1, p_i) + w(p_k - 1, p_j + 1, p_i) - w(p_k, p_j, p_i) - w(p_k - 1, p_j, p_i)) / dy(p_j) &
         + (w(p_k, p_j, p_i) + w(p_k - 1, p_j, p_i) - w(p_k, p_j - 1, p_i) - w(p_k - 1, p_j - 1, p_i)) / dy(p_j - 1))

        delta = 0.5 * (w(p_k, p_j, p_i) + w(p_k - 1, p_j, p_i) &
         - alpha * (ddz(p_k) - dz(p_k - 1)) - beta * (ddx(p_i) - dx(p_i - 1)) - gamma * (ddy(p_j) - dy(p_j - 1)))

        p_w = alpha * (p_z - z(p_k)) + beta * (p_x - x(p_i)) + gamma * (p_y - y(p_j)) + delta

        CALL stop_timer(923)

    END SUBROUTINE gobert_particle_uvw

    !------------------------------

    ! THIS ROUTINE IS NOT SUITABLE AS IT IS NOT MASS CONSERVATIVE (divergence of the interpolated velocity is not zero)
    SUBROUTINE trilinear_particle_uvw(kk, jj, ii, particle, p_u, p_v, p_w, &
         xstag, ystag, zstag, u, v, w, x, y, z)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        TYPE(baseparticle_t), INTENT(in) :: particle
        REAL(realk), INTENT(out) :: p_u, p_v, p_w
        REAL(realk), INTENT(in) :: xstag(ii), ystag(jj), zstag(kk)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)

        ! local variables
        INTEGER(intk) :: p_i, p_j, p_k, p_istag, p_jstag, p_kstag
        REAL(realk) :: p_x, p_y, p_z
        REAL(realk), DIMENSION(6) :: p_bbox_u, p_bbox_v, p_bbox_w

        !just for readability of the following expressions
        p_i = particle%ijkcell(1)
        p_j = particle%ijkcell(2)
        p_k = particle%ijkcell(3)

        !just for readability of the following expressions
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

    END SUBROUTINE trilinear_particle_uvw

    !------------------------------

    !Trilinear interpolation for a rectilinear grid !
    SUBROUTINE interp_trilinear(x, y, z, bbox, &
        f, f19, f20, f21, f22, f23, f24, f25, f26)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), DIMENSION(6), INTENT(in) :: bbox
        ! numbers indicate the corners (see convention on the numbering of faces/edges/corners)
        REAL(realk), INTENT(in) :: f19, f20, f21, f22, f23, f24, f25, f26
        REAL(realk), INTENT(out) :: f

        ! local variables
        REAL(realk) :: xd, yd, zd

        xd = (x - bbox(1)) / (bbox(2) - bbox(1))
        yd = (y - bbox(3)) / (bbox(4) - bbox(3))
        zd = (z - bbox(5)) / (bbox(6) - bbox(5))

        f = ((f19 * (1 - xd) + f23 * xd)  * (1 - yd) + (f21 * (1 - xd) + f25 * xd) * yd) * (1 - zd) + &
            ((f20 * (1 - xd) + f24 * xd) * (1 - yd) + (f22 * (1 - xd) + f26 * xd) * yd) * zd

    END SUBROUTINE interp_trilinear

    SUBROUTINE time_interpolate_field(brk, u_f, v_f, w_f, ubu_f, vbu_f, wbu_f, uli_f, vli_f, wli_f)

        ! subroutine arguments
        REAL(realk), INTENT(in) :: brk
        TYPE(field_t), INTENT(in) :: u_f, v_f, w_f
        TYPE(field_t), INTENT(in) :: ubu_f, vbu_f, wbu_f
        TYPE(field_t), INTENT(inout) :: uli_f, vli_f, wli_f

        ! local variables
        INTEGER(intk) :: g, i, j, k, ii, jj, kk, igrid
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ubu, vbu, wbu
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uli, vli, wli

        CALL start_timer(921)

        DO g = 1, nmygrids
            igrid = mygrids(g)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            CALL ubu_f%get_ptr(ubu, igrid)
            CALL vbu_f%get_ptr(vbu, igrid)
            CALL wbu_f%get_ptr(wbu, igrid)

            CALL uli_f%get_ptr(uli, igrid)
            CALL vli_f%get_ptr(vli, igrid)
            CALL wli_f%get_ptr(wli, igrid)

            DO i = 1, ii
                DO j = 1, jj
                    DO k = 1, kk
                        uli(k, j, i) = ubu(k, j, i) * (1 - brk) + u(k, j, i) * brk
                        vli(k, j, i) = vbu(k, j, i) * (1 - brk) + v(k, j, i) * brk
                        wli(k, j, i) = wbu(k, j, i) * (1 - brk) + w(k, j, i) * brk
                    END DO
                END DO
            END DO

        END DO

        CALL stop_timer(921)

    END SUBROUTINE time_interpolate_field

END MODULE particle_interpolation_mod