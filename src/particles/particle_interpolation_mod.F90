MODULE particle_interpolation_mod

    ! This module is responsible for:
    ! Interpolation of (staggered) fields at given point coordinates (particle location)

    USE particle_core_mod

    IMPLICIT NONE

CONTAINS    !===================================

    ! get nearest vector components of staggered component
    ! no interpolation of vector components, just the values of the nearest staggered cells respectively
    SUBROUTINE get_nearest_value(particle, kk, jj, ii, x, y, z, &
         v1, v2, v3, p_v1, p_v2, p_v3)

        !subroutine_arguments
        TYPE(baseparticle_t), INTENT(in) :: particle
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: v1(kk, jj, ii), v2(kk, jj, ii), v3(kk, jj, ii)
        REAL(realk), INTENT(out) :: p_v1, p_v2, p_v3

        !local variables
        INTEGER(intk) :: p_istag, p_jstag, p_kstag

        p_istag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%x - x(particle%ijkcell(1)) ), intk ))
        p_jstag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%y - x(particle%ijkcell(2)) ), intk ))
        p_kstag = MAX( 0_intk, NINT( SIGN( 1.0_realk, particle%z - x(particle%ijkcell(3)) ), intk ))

        p_v1 = v1(particle%ijkcell(3), particle%ijkcell(2), particle%ijkcell(1) + p_istag - 1)
        p_v2 = v2(particle%ijkcell(3), particle%ijkcell(2) + p_jstag - 1, particle%ijkcell(1))
        p_v3 = v3(particle%ijkcell(3) + p_kstag - 1, particle%ijkcell(2), particle%ijkcell(1))

    END SUBROUTINE get_nearest_value

    ! linear and differntially conservative interpolation of staggered vector field
    ! from Gobert et. al, LAGRANGIAN SCALAR TRACKING FOR LAMINAR MICROMIXING AT HIGH SCHMIDT NUMBERS, 2006
    SUBROUTINE interpolate_lincon(particle, kk, jj, ii, x, y, z, dx, dy, dz, ddx, ddy, ddz, &
         v1, v2, v3, p_v1, p_v2, p_v3)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(in) :: particle
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk), dx(ii), dy(jj), dz(kk), ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: v1(kk, jj, ii), v2(kk, jj, ii), v3(kk, jj, ii)
        REAL(realk), INTENT(out) :: p_v1, p_v2, p_v3 ! particle values in x/y/z direction (velocity or diffusion constant)

        ! local variables
        INTEGER(intk) :: p_i, p_j, p_k
        INTEGER(intk) :: p_ip, p_im, p_jp, p_jm, p_kp, p_km
        REAL(realk) :: p_x, p_y, p_z, alpha, beta, gamma, delta

        !just for readability of the following expressions
        p_i = particle%ijkcell(1)
        p_j = particle%ijkcell(2)
        p_k = particle%ijkcell(3)

        ! limit the shifted indices to the grid limits!
        ! if a shifted index would be outside the grid index range of [1, ii] (or [1, jj], [1, kk])
        ! it is set to be the respective grid index limit.
        ! that should be equivalent to just giving any cell that is referred but is "located outside the grid" (i.e. does not exist)
        ! the velocity value of the "nearest" cell on the grid (i.e. cell that acutally exists)
        p_im = MAX(MIN(p_i - 1, ii), 1)
        p_ip = MAX(MIN(p_i + 1, ii), 1)
        p_jm = MAX(MIN(p_j - 1, jj), 1)
        p_jp = MAX(MIN(p_j + 1, jj), 1)
        p_km = MAX(MIN(p_k - 1, kk), 1)
        p_kp = MAX(MIN(p_k + 1, kk), 1)

        !just for readability of the following expressions
        p_x = particle%x
        p_y = particle%y
        p_z = particle%z

        ! u interpolation
        alpha = (v1(p_k, p_j, p_i) - v1(p_k, p_j, p_im)) / ddx(p_i)

        beta = 0.25 * ((v1(p_k, p_jp, p_i) + v1(p_k, p_jp, p_im) - v1(p_k, p_j, p_i) - v1(p_k, p_j, p_im)) / dy(p_j) &
         + (v1(p_k, p_j, p_i) + v1(p_k, p_j, p_im) - v1(p_k, p_jm, p_i) - v1(p_k, p_jm, p_im)) / dy(p_j -1))

        gamma = 0.25 * ((v1(p_kp, p_j, p_i) + v1(p_kp, p_j, p_im) - v1(p_k, p_j, p_i) - v1(p_k, p_j, p_im)) / dz(p_k) &
         + (v1(p_k, p_j, p_i) + v1(p_k, p_j, p_im) - v1(p_km, p_j, p_i) - v1(p_km, p_j, p_im)) / dz(p_km))

        delta = 0.5 * (v1(p_k, p_j, p_i) + v1(p_k, p_j, p_im) &
         - alpha * (ddx(p_i) - dx(p_im)) - beta * (ddy(p_j) - dy(p_jm)) - gamma * (ddz(p_k) - dz(p_km)))

        p_v1 = alpha * (p_x - x(p_i)) + beta * (p_y - y(p_j)) + gamma * (p_z - z(p_k)) + delta

        ! v interpolation
        alpha = (v2(p_k, p_j, p_i) - v2(p_k, p_jm, p_i)) / ddy(p_j)

        beta = 0.25 * ((v2(p_k, p_j, p_ip) + v2(p_k, p_jm, p_ip) - v2(p_k, p_j, p_i) - v2(p_k, p_jm, p_i)) / dx(p_i) &
         + (v2(p_k, p_j, p_i) + v2(p_k, p_jm, p_i) - v2(p_k, p_j, p_im) - v2(p_k, p_jm, p_im)) /dx(p_im))

        gamma = 0.25 * ((v2(p_kp, p_j, p_i) + v2(p_kp, p_jm, p_i) - v2(p_k, p_j, p_i) - v2(p_k, p_jm, p_i)) / dz(p_k) &
         + (v2(p_k, p_j, p_i) + v2(p_k, p_jm, p_i) - v2(p_km, p_j, p_i) - v2(p_km, p_jm, p_i)) / dz(p_km))

        delta = 0.5 * (v2(p_k, p_j, p_i) + v2(p_k, p_jm, p_i) &
         - alpha * (ddy(p_j) - dy(p_jm)) - beta * (ddx(p_i) - dx(p_im)) - gamma * (ddz(p_k) - dz(p_km)))

        p_v2 = alpha * (p_y - y(p_j)) + beta * (p_x - x(p_i)) + gamma * (p_z - z(p_k)) + delta

        ! w interpolation
        alpha = (v3(p_k, p_j, p_i) - v3(p_km ,p_j ,p_i)) / ddz(p_k)

        beta = 0.25 * ((v3(p_k, p_j, p_ip) + v3(p_km, p_j, p_ip) - v3(p_k, p_j, p_i) - v3(p_km, p_j, p_i)) / dx(p_i) &
         + (v3(p_k, p_j, p_i) + v3(p_km, p_j, p_i) - v3(p_k, p_j, p_im) - v3(p_km, p_j, p_im)) / dx(p_im))

        gamma = 0.25 * ((v3(p_k, p_jp, p_i) + v3(p_km, p_jp, p_i) - v3(p_k, p_j, p_i) - v3(p_km, p_j, p_i)) / dy(p_j) &
         + (v3(p_k, p_j, p_i) + v3(p_km, p_j, p_i) - v3(p_k, p_jm, p_i) - v3(p_km, p_jm, p_i)) / dy(p_jm))

        delta = 0.5 * (v3(p_k, p_j, p_i) + v3(p_km, p_j, p_i) &
         - alpha * (ddz(p_k) - dz(p_km)) - beta * (ddx(p_i) - dx(p_im)) - gamma * (ddy(p_j) - dy(p_jm)))

        p_v3 = alpha * (p_z - z(p_k)) + beta * (p_x - x(p_i)) + gamma * (p_y - y(p_j)) + delta

    END SUBROUTINE interpolate_lincon


END MODULE particle_interpolation_mod