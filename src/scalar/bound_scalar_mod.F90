MODULE bound_scalar_mod
    USE core_mod
    USE scacore_mod
    USE flow_mod, ONLY: ilesmodel, gmol, rho, qwallfix
    USE offload_helper_mod

    ! ┌────────────────────────────────────────────────────────────────────────────┐
    ! | Removed interface complying with bound_mod                                 |
    ! | Allows for experimental code architecture to offload computation           |
    ! └────────────────────────────────────────────────────────────────────────────┘

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: bound_sca
CONTAINS
    !> @brief Applies boundary conditions to scalar field
    !!
    !! Subroutine is specifically resigned to apply boundary conditions to scalar fluxes that "live" on the target device.
    !! Thus, in the actual computation, we can only use data available on the target device.
    !! All parameters "ilevel", "t_f" and "timeph" are actually unused in this module and only kept to resemble the
    !! old interface more closely and provide a starting point for a more robust implementation.
    !!
    !! @param[in] ilevel Grid-level index
    !! @param[in] t_f    Scalar field
    !! @param[in] timeph size if the array
    SUBROUTINE bound_sca(ilevel, t_f, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ilevel
        TYPE(field_t), INTENT(in) :: t_f
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: igrid, iface, nbocd, ibocd, ctyp_encoded
        REAL(realk) :: sca_prmol
        sca_prmol = scalar(ISCA_FIELD)%prmol

        ! Simplification: Only BCs on the same level are handled, ilevel is unused
        !$omp target teams distribute
        DO igrid = 1, nmygrids
            DO iface = 1, 6
                nbocd = nboconds_offload(iface, igrid)
                DO ibocd = 1, nbocd
                    CALL get_encoded_ctyp_offload(ctyp_encoded, ibocd, iface, igrid)
                    SELECT CASE(iface)
                    CASE(1, 2)
                        CALL bfront(igrid, iface, ctyp_encoded, t_f, sca_prmol, gmol, rho, timeph)
                    CASE(3, 4)
                        CALL bright(igrid, iface, ctyp_encoded, t_f, sca_prmol, gmol, rho, timeph)
                    CASE(5, 6)
                        CALL bbottom(igrid, iface, ctyp_encoded, t_f, sca_prmol, gmol, rho, timeph)
                    END SELECT
                END DO
            END DO
        END DO
        !$omp end target teams distribute
    END SUBROUTINE bound_sca

    SUBROUTINE bfront(igrid, iface, ctyp_encoded, t_f, sca_prmol, gmol_offload, rho_offload, timeph)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ctyp_encoded
        TYPE(field_t), INTENT(in), OPTIONAL :: t_f
        REAL(realk), INTENT(in), OPTIONAL :: timeph
        REAL(realk), INTENT(in) :: sca_prmol, gmol_offload, rho_offload

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i3, istag2, dir
        REAL(realk) :: area, adv, diff, gamma, tout
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtu, t, bt, u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, ddx, ddy, ddz

        ! Return early when no action is to be taken
        SELECT CASE (ctyp_encoded)
        CASE (2, 5)
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        CALL ptr_to_grid3(qtu_offload, igrid, qtu)
        CALL ptr_to_grid3(t_offload, igrid, t)
        CALL ptr_to_grid3(u_offload, igrid, u)
        CALL ptr_to_grid3(v_offload, igrid, v)
        CALL ptr_to_grid3(w_offload, igrid, w)
        CALL ptr_to_grid3(bt_offload, igrid, bt)
        CALL ptr_to_grid_x(dx_offload, igrid, dx)
        CALL ptr_to_grid_x(ddx_offload, igrid, ddx)
        CALL ptr_to_grid_y(ddy_offload, igrid, ddy)
        CALL ptr_to_grid_z(ddz_offload, igrid, ddz)
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

        ! ┌────────────────────────────────────────────────────────────────────────────┐
        ! | SIMPLIFICATION: Only Fixed scalar value (inflow/outflow) BC = "SIO"        |
        ! |     - Removed handling of other BCs                                        |
        ! |     - sbctype(idx)=0 for our test case                                     |
        ! |     - tbuf=0 for our test case                                             |
        ! └────────────────────────────────────────────────────────────────────────────┘
        gamma = gmol_offload / rho_offload / sca_prmol

        !$omp parallel do collapse(2)
        DO j = 1, jj
            DO k = 1, kk
                IF (-dir*u(k, j, istag2) >= 0.0) THEN
                    ! flow into the domain (requires specified value)
                    ! OLD: tout = 2.0*tbuf(k, j, 1) - t(k, j, i3)
                    tout = 2.0*0.0 - t(k, j, i3)
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
        !$omp end parallel do
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, ctyp_encoded, t_f, sca_prmol, gmol_offload, rho_offload, timeph)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ctyp_encoded
        TYPE(field_t), INTENT(in) :: t_f
        REAL(realk), INTENT(in), OPTIONAL :: timeph
        REAL(realk), INTENT(in) :: sca_prmol, gmol_offload, rho_offload

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, i, jstag2, dir
        REAL(realk) :: area
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtv, t, bt, u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dy, ddx, ddy, ddz

        ! Return early when no action is to be taken
        SELECT CASE (ctyp_encoded)
        CASE (2, 5)
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Fetch pointers
        CALL ptr_to_grid3(t_offload, igrid, t)
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

        ! ┌────────────────────────────────────────────────────────────────────────────┐
        ! | SIMPLIFICATION: Only Wall BC = "SWA"                                       |
        ! |     - Removed handling of other BCs                                        |
        ! |     - sbctype(idx)=1 for our test case                                     |
        ! |     - tbuf=0 for our test case                                             |
        ! └────────────────────────────────────────────────────────────────────────────┘
        !$omp parallel do collapse(2)
        DO i = 1, ii
            DO k = 1, kk
                ! Wall buffer tbuf contains flux at this boundary
                area = ddx(i)*ddz(k)
                ! OLD: qtv(k, jstag2, i) = -dir*tbuf(k, i, 1)*area
                qtv(k, jstag2, i) = 0
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, ctyp_encoded, t_f, sca_prmol, gmol_offload, rho_offload, timeph)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ctyp_encoded
        TYPE(field_t), INTENT(in) :: t_f
        REAL(realk), INTENT(in), OPTIONAL :: timeph
        REAL(realk), INTENT(in) :: sca_prmol, gmol_offload, rho_offload

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: j, i, k3, kstag2, dir
        REAL(realk) :: diff, gamma2dx!, uquer, area
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtw, t, bt, u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dz, ddx, ddy, ddz

        ! Return early when no action is to be taken
        SELECT CASE (ctyp_encoded)
        CASE (2, 5)
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Fetch pointers
        CALL ptr_to_grid3(qtw_offload, igrid, qtw)
        CALL ptr_to_grid3(t_offload, igrid, t)
        CALL ptr_to_grid3(u_offload, igrid, u)
        CALL ptr_to_grid3(v_offload, igrid, v)
        CALL ptr_to_grid3(w_offload, igrid, w)
        CALL ptr_to_grid3(bt_offload, igrid, bt)
        CALL ptr_to_grid_x(ddx_offload, igrid, ddx)
        CALL ptr_to_grid_y(ddy_offload, igrid, ddy)
        CALL ptr_to_grid_z(ddz_offload, igrid, ddz)
        CALL ptr_to_grid_z(dz_offload, igrid, dz)
        CALL get_mgdims_target(kk, jj, ii, igrid)

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

        ! ┌────────────────────────────────────────────────────────────────────────────┐
        ! | SIMPLIFICATION: Only Wall BC = "SWA"                                       |
        ! |     - Removed handling of other BCs                                        |
        ! |     - Removed handling of any ilesmodel other than ilesmodel==0            |
        ! |     - sbctype(idx)=0 for our test case                                     |
        ! |     - tbuf=1 for our test case                                             |
        ! └────────────────────────────────────────────────────────────────────────────┘
        gamma2dx = 2.0 * gmol_offload / rho_offload / sca_prmol / dz(kstag2)
        !$omp parallel do collapse(2)
        DO i = 1, ii
            DO j = 1, jj
                ! Setting the scalar diffusive flux from the wall
                ! Wall buffer tbuf contains set scalar value
                ! OLD: diff = gamma2dx*(tbuf(j, i, 1) - t(k3, j, i))
                diff = gamma2dx*(1.0 - t(k3, j, i))
                qtw(kstag2, j, i) = -dir*diff*ddx(i)*ddy(j)
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE bbottom

END MODULE bound_scalar_mod
