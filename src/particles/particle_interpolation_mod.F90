MODULE particle_interpolation_mod

    ! This module is responsible for:
    ! Interpolation of (staggered) fields at given point coordinates (particle location)

    USE particle_basetype_mod
    USE particle_ofields_mod

    IMPLICIT NONE

CONTAINS    !===================================

    ! get nearest vector components of staggered component
    ! no interpolation of vector components, just the values of the nearest staggered cells respectively
    SUBROUTINE get_nearest_value(particle, kk, jj, ii, x, y, z, &
         v1, v2, v3, p_v1, p_v2, p_v3)

        !$omp declare target

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

        !$omp declare target

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

        ! TODO: FIX GRADIENTS AT NO FLUX GRID BOUNDARIES !!!

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
         + (v1(p_k, p_j, p_i) + v1(p_k, p_j, p_im) - v1(p_k, p_jm, p_i) - v1(p_k, p_jm, p_im)) / dy(p_jm))

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


    ! linear and differntially conservative interpolation of staggered vector field
    ! from Gobert et. al, LAGRANGIAN SCALAR TRACKING FOR LAMINAR MICROMIXING AT HIGH SCHMIDT NUMBERS, 2006
    SUBROUTINE interpolate_lincon_target(particle, igrid, kk, jj, ii, p_v1, p_v2, p_v3)

        !$omp declare target

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(in) :: particle
        INTEGER(intk), INTENT(in) :: igrid, kk, jj, ii
        REAL(realk), INTENT(out) :: p_v1, p_v2, p_v3 ! particle values in x/y/z direction (velocity or diffusion constant)

        ! local variables
        INTEGER(intk) :: ip3d, ip1d(3), i, j, k
        REAL(realk) :: x, ddx, y, ddy, z, ddz
        REAL(realk) :: dx(-1:0), dy(-1:0), dz(-1:0)
        REAL(realk) :: v1(-1:1,-1:1,-1:0), v2(-1:1,-1:0,-1:1), v3(-1:0,-1:1,-1:1)
        INTEGER(intk) :: p_i(-1:1), p_j(-1:1), p_k(-1:1)
        REAL(realk) :: alpha, beta, gamma, delta

        ! TODO: FIX GRADIENTS AT NO FLUX GRID BOUNDARIES !!!

        p_i(-1) = MAX(MIN(particle%ijkcell(1) - 1, ii), 1)
        p_i( 0) = particle%ijkcell(1)
        p_i( 1) = MAX(MIN(particle%ijkcell(1) + 1, ii), 1)

        p_j(-1) = MAX(MIN(particle%ijkcell(2) - 1, jj), 1)
        p_j( 0) = particle%ijkcell(2)
        p_j( 1) = MAX(MIN(particle%ijkcell(2) + 1, jj), 1)
        
        p_k(-1) = MAX(MIN(particle%ijkcell(3) - 1, kk), 1)
        p_k( 0) = particle%ijkcell(3)
        p_k( 1) = MAX(MIN(particle%ijkcell(3) + 1, kk), 1)

        ip1d(:) = ip1d_offload(:,igrid)

        x = x_offload(ip1d(1) + particle%ijkcell(1) - 1) 
        dx(-1) = dx_offload(ip1d(1) + p_i(-1) - 1)
        dx( 0) = dx_offload(ip1d(1) + particle%ijkcell(1) -1)
        ddx = ddx_offload(ip1d(1) + particle%ijkcell(1) - 1)

        y = y_offload(ip1d(2) + particle%ijkcell(2) - 1)
        dy(-1) = dy_offload(ip1d(2) + p_j(-1) - 1)
        dy( 0) = dy_offload(ip1d(2) + particle%ijkcell(2) -1)
        ddy = ddy_offload(ip1d(2) + particle%ijkcell(2) - 1)

        z = z_offload(ip1d(3) + particle%ijkcell(3) - 1) 
        dz(-1) = dz_offload(ip1d(3) + p_k(-1) - 1)
        dz( 0) = dz_offload(ip1d(3) + particle%ijkcell(3) -1)
        ddz = ddz_offload(ip1d(3) + particle%ijkcell(3) - 1)

        ip3d = ip3d_offload(igrid)

        DO i=-1, 0
            DO j=-1, 1
                DO k=-1, 1
                    v1(k, j, i) = u_offload(ip3d + &
                    (p_k(k) - 1) + &
                    (p_j(j) - 1) * kk + &
                    (p_i(i) - 1) * kk * jj)
                END DO
            END DO
        END DO 
        
        DO i=-1, 1
            DO j=-1, 0
                DO k=-1, 1
                    v2(k, j, i) = v_offload(ip3d + &
                    (p_k(k) - 1) + &
                    (p_j(j) - 1) * kk + &
                    (p_i(i) - 1) * kk * jj)
                END DO
            END DO
        END DO

        DO i=-1, 1
            DO j=-1, 1
                DO k=-1, 0
                    v3(k, j, i) = w_offload(ip3d + &
                    (p_k(k) - 1) + &
                    (p_j(j) - 1) * kk + &
                    (p_i(i) - 1) * kk * jj)
                END DO
            END DO
        END DO 

        ! u interpolation
        alpha = (v1(0, 0, 0) - v1(0, 0, -1)) / ddx

        beta = 0.25 * ((v1(0, 1, 0) + v1(0, 1, -1) - v1(0, 0, 0) - v1(0, 0, -1)) / dy(0) &
         + (v1(0, 0, 0) + v1(0, 0, -1) - v1(0, -1, 0) - v1(0, -1, -1)) / dy(-1))

        gamma = 0.25 * ((v1(1, 0, 0) + v1(1, 0, -1) - v1(0, 0, 0) - v1(0, 0, -1)) / dz(0) &
         + (v1(0, 0, 0) + v1(0, 0, -1) - v1(-1, 0, 0) - v1(-1, 0, -1)) / dz(-1))

        delta = 0.5 * (v1(0, 0, 0) + v1(0, 0, -1) &
         - alpha * (ddx - dx(-1)) - beta * (ddy - dy(-1)) - gamma * (ddz - dz(-1)))

        p_v1 = alpha * (particle%x - x) + beta * (particle%y - y) + gamma * (particle%z - z) + delta

        ! v interpolation
        alpha = (v2(0, 0, 0) - v2(0, -1, 0)) / ddy

        beta = 0.25 * ((v2(0, 0, 1) + v2(0, -1, 1) - v2(0, 0, 0) - v2(0, -1, 0)) / dx(0) &
         + (v2(0, 0, 0) + v2(0, -1, 0) - v2(0, 0, -1) - v2(0, -1, -1)) /dx(-1))

        gamma = 0.25 * ((v2(1, 0, 0) + v2(1, -1, 0) - v2(0, 0, 0) - v2(0, -1, 0)) / dz(0) &
         + (v2(0, 0, 0) + v2(0, -1, 0) - v2(-1, 0, 0) - v2(-1, -1, 0)) / dz(-1))

        delta = 0.5 * (v2(0, 0, 0) + v2(0, -1, 0) &
         - alpha * (ddy - dy(-1)) - beta * (ddx - dx(-1)) - gamma * (ddz - dz(-1)))

        p_v2 = alpha * (particle%y - y) + beta * (particle%x - x) + gamma * (particle%z - z) + delta

        ! w interpolation
        alpha = (v3(0, 0, 0) - v3(-1 ,0 ,0)) / ddz

        beta = 0.25 * ((v3(0, 0, 1) + v3(-1, 0, 1) - v3(0, 0, 0) - v3(-1, 0, 0)) / dx(0) &
         + (v3(0, 0, 0) + v3(-1, 0, 0) - v3(0, 0, -1) - v3(-1, 0, -1)) / dx(-1))

        gamma = 0.25 * ((v3(0, -1, 0) + v3(-1, -1, 0) - v3(0, 0, 0) - v3(-1, 0, 0)) / dy(0) &
         + (v3(0, 0, 0) + v3(-1, 0, 0) - v3(0, -1, 0) - v3(-1, -1, 0)) / dy(-1))

        delta = 0.5 * (v3(0, 0, 0) + v3(-1, 0, 0) &
         - alpha * (ddz - dz(-1)) - beta * (ddx - dx(-1)) - gamma * (ddy - dy(-1)))

        p_v3 = alpha * (particle%z - z) + beta * (particle%x - x) + gamma * (particle%y - y) + delta

    END SUBROUTINE interpolate_lincon_target



END MODULE particle_interpolation_mod