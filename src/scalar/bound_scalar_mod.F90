MODULE bound_scalar_mod
    USE core_mod
    USE scacore_mod
    USE flow_mod, ONLY: ilesmodel, gmol, rho, qwallfix
    USE offload_helper_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: bound_sca

CONTAINS
    SUBROUTINE bound_sca(ilevel, t_f, sca_prmol, timeph)
        INTEGER(intk), INTENT(in) :: ilevel
        REAL(realk), INTENT(in) :: sca_prmol
        TYPE(field_t), INTENT(in) :: t_f
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        INTEGER(intk) :: igrid, iface

        ! Simplification: Only BCs on the same level are handled
        DO igrid = 1, nmygrids

            DO iface = 1, 6
                SELECT CASE(iface)
                CASE(1)
                    CALL bfront(igrid, iface, t_f, sca_prmol, timeph)
                CASE(2)
                    CALL bfront(igrid, iface, t_f, sca_prmol, timeph)
                CASE(3)
                    CALL bright(igrid, iface, t_f, sca_prmol, timeph)
                CASE(4)
                    CALL bright(igrid, iface, t_f, sca_prmol, timeph)
                CASE(5)
                    CALL bbottom(igrid, iface, t_f, sca_prmol, timeph)
                CASE(6)
                    CALL bbottom(igrid, iface, t_f, sca_prmol, timeph)
                END SELECT
            END DO
        END DO
    END SUBROUTINE bound_sca

    SUBROUTINE bfront(igrid, iface, t_f, sca_prmol, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface
        TYPE(field_t), INTENT(in), OPTIONAL :: t_f
        REAL(realk), INTENT(in), OPTIONAL :: timeph
        REAL(realk), INTENT(in) :: sca_prmol

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i3, istag2, dir
        REAL(realk) :: area, adv, diff, gamma, tout

        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtu, t, bt, u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, ddx, ddy, ddz

        CALL ptr_to_grid3(qtu_offload, igrid, qtu)

        CALL t_f%get_ptr(t, igrid)

        CALL ptr_to_grid3(u_offload, igrid, u)
        CALL ptr_to_grid3(v_offload, igrid, v)
        CALL ptr_to_grid3(w_offload, igrid, w)
        CALL ptr_to_grid3(bt_offload, igrid, bt)

        CALL ptr_to_grid_x(ddx_offload, igrid, ddx)
        CALL ptr_to_grid_y(ddy_offload, igrid, ddy)
        CALL ptr_to_grid_z(ddz_offload, igrid, ddz)

        CALL ptr_to_grid_x(dx_offload, igrid, dx)

        CALL get_mgdims_target(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (1)
            ! Front
            i3 = 3
            istag2 = 2
            dir = -1
        CASE (2)
            ! Back
            i3 = ii - 2
            istag2 = ii - 2
            dir = 1
        END SELECT

        gamma = gmol / rho / sca_prmol

        !$omp target teams distribute parallel do collapse(2)
        DO j = 1, jj
            DO k = 1, kk
                IF (-dir*u(k, j, istag2) >= 0.0) THEN
                    ! flow into the domain (requires specified value)
                    ! OLD: tout = 2.0*tbuf(k, j, 1) - t(k, j, i3)
                    tout = -t(k, j, i3)
                ELSE
                    ! flow out of the domain (zero-gradient)
                    tout = t(k, j, i3)
                END IF

                ! Adding to the flux (diffusive and convective)
                ! adv and diff are fluxes /into/ domain!
                adv = -dir*u(k, j, istag2)*0.5*(tout + t(k, j, i3))
                diff = gamma*(tout - t(k, j, i3))/dx(istag2)
                area = ddy(j)*ddz(k)

                ! flux in /coordinate direction/ (ref. fluxbalance)
                qtu(k, j, istag2) = -dir*(adv + diff)*area
            END DO
        END DO
        !$omp end target teams distribute parallel do
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, t_f, sca_prmol, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface
        TYPE(field_t), INTENT(in) :: t_f
        REAL(realk), INTENT(in), OPTIONAL :: timeph
        REAL(realk), INTENT(in) :: sca_prmol

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, i, jstag2, dir
        REAL(realk) :: area
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtv, t, bt, u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dy, ddx, ddy, ddz

        ! Fetch pointers
        CALL t_f%get_ptr(t, igrid)

        CALL ptr_to_grid3(qtv_offload, igrid, qtv)

        CALL ptr_to_grid3(u_offload, igrid, u)
        CALL ptr_to_grid3(v_offload, igrid, v)
        CALL ptr_to_grid3(w_offload, igrid, w)
        CALL ptr_to_grid3(bt_offload, igrid, bt)

        CALL ptr_to_grid_x(ddx_offload, igrid, ddx)
        CALL ptr_to_grid_y(ddy_offload, igrid, ddy)
        CALL ptr_to_grid_z(ddz_offload, igrid, ddz)

        CALL ptr_to_grid_y(dy_offload, igrid, dy)

        CALL get_mgdims_target(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (3)
            ! Right
            jstag2 = 2
            dir = -1
        CASE (4)
            ! Left
            jstag2 = jj - 2
            dir = 1
        END SELECT

        ! SWA ctyp with sbctype(idx)=1
        DO i = 1, ii
            DO k = 1, kk
                ! Wall buffer tbuf contains flux at this boundary
                area = ddx(i)*ddz(k)
                ! OLD: qtv(k, jstag2, i) = -dir*tbuf(k, i, 1)*area
                qtv(k, jstag2, i) = 0
            END DO
        END DO
    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, t_f, sca_prmol, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface
        TYPE(field_t), INTENT(in) :: t_f
        REAL(realk), INTENT(in), OPTIONAL :: timeph
        REAL(realk), INTENT(in) :: sca_prmol

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: j, i, k3, kstag2, dir
        REAL(realk) :: area, diff, gamma2dx, uquer
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtw, t, bt, u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dz, ddx, ddy, ddz

        ! Fetch pointers
        CALL ptr_to_grid3(qtw_offload, igrid, qtw)

        CALL t_f%get_ptr(t, igrid)

        CALL ptr_to_grid3(u_offload, igrid, u)
        CALL ptr_to_grid3(v_offload, igrid, v)
        CALL ptr_to_grid3(w_offload, igrid, w)
        CALL ptr_to_grid3(bt_offload, igrid, bt)

        CALL ptr_to_grid_x(ddx_offload, igrid, ddx)
        CALL ptr_to_grid_y(ddy_offload, igrid, ddy)
        CALL ptr_to_grid_z(ddz_offload, igrid, ddz)

        CALL ptr_to_grid_z(dz_offload, igrid, dz)

        CALL get_mgdims_target(kk, jj, ii, igrid)

        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (5)
            ! Bottom
            k3 = 3
            kstag2 = 2
            dir = -1
        CASE (6)
            ! Top
            k3 = kk - 2
            kstag2 = kk - 2
            dir = 1
        END SELECT


        ! SWA ctyp with sbctype(idx)=0
        ! tbuf=1.0
        IF (ilesmodel == 0) THEN
            gamma2dx = 2.0 * gmol / rho / sca_prmol / dz(kstag2)
            DO i = 1, ii
                DO j = 1, jj
                    ! Setting the scalar diffusive flux from the wall
                    ! Wall buffer tbuf contains set scalar value
                    ! OLD: diff = gamma2dx*(tbuf(j, i, 1) - t(k3, j, i))
                    diff = gamma2dx*(1.0 - t(k3, j, i))
                    qtw(kstag2, j, i) = -dir*diff*ddx(i)*ddy(j)
                END DO
            END DO
        ELSE
            DO i = 2, ii
                DO j = 2, jj
                    ! Setting the scalar flux with a wall model
                    ! Wall buffer tbuf contains set scalar value
                    area = ddx(i)*ddy(j)
                    uquer = SQRT( &
                        (u(k3, j, i-1) + (u(k3, j, i)-u(k3, j, i-1)) &
                            /ddx(i)*ddx(i-1)*0.5)**2.0 + &
                        (v(k3, j-1, i) + (v(k3, j, i)-v(k3, j-1, i)) &
                            /ddy(j)*ddy(j-1)*0.5 )**2.0)
                    ! OLD: qtw(kstag2, j, i) = -dir*qwallfix(tbuf(j, i, 1), &
                    ! OLD:    t(k3, j, i), uquer, ddz(k3), prmol)*area
                    qtw(kstag2, j, i) = -dir*qwallfix(1.0, &
                        t(k3, j, i), uquer, ddz(k3), sca_prmol)*area
                END DO
            END DO
        END IF
    END SUBROUTINE bbottom

END MODULE bound_scalar_mod
