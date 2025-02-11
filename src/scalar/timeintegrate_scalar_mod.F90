MODULE timeintegrate_scalar_mod
    USE core_mod
    USE ib_mod, ONLY: parent, ftoc, ib
    USE scacore_mod
    USE flow_mod, ONLY: ilesmodel, gmol, rho
    USE bound_scalar_mod
    USE itinfo_scalar_mod
    USE gc_scastencils_mod
    USE offload_helper_mod
    USE pointers_mod
    USE boussinesqterm_mod, ONLY: has_buoyancy

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: timeintegrate_scalar, itinfo_scalar

CONTAINS
    SUBROUTINE timeintegrate_scalar(itstep, ittot, timeph, dt, irk, rkscheme)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: irk
        TYPE(rk_2n_t), INTENT(in) :: rkscheme

        ! Local variables
        INTEGER(intk) :: ilevel, l
        REAL(realk) :: frhs, fu, dtrk, dtrki
        TYPE(field_t), POINTER :: t, told, dt_f
        TYPE(field_t) :: qtt

        IF (.NOT. solve_scalar) RETURN
        CALL start_timer(400)
        
        ! Scalar buffer allocation
        CALL start_timer(401)
        ! Local temporary storage ("scrap")
        CALL qtt%init("QTT")
        CALL stop_timer(401)

        ! Update u, v and w fields on the device for next scalar timestep
        ! These were solved for in the flow simulation on the host
        !$omp target update to(u_offload, v_offload, w_offload)

        ! Perform actual time integration
        ! ┌────────────────────────────────────────────────────────────────────────────┐
        ! | SIMPLIFICATION: Only one scalar field                                      |
        ! | *Removed loop over multiple scalar fields*                                 |
        ! | Allows us to simply keep a mapped version of the scalar field              |
        ! | Saves and extra mapping of the t field                                     |
        ! └────────────────────────────────────────────────────────────────────────────┘
        CALL start_timer(402)
        DO l = 1, ISCA_FIELD
            ! Fetch scalar field
            CALL get_field(t, scalar(l)%name)
            CALL get_field(dt_f, "D"//TRIM(scalar(l)%name))
            ! Copy to "T_OLD" 
            ! ONLY REQUIRED FOR Boussinesq-Approximation
            IF (has_buoyancy) THEN
                CALL get_field(told, TRIM(scalar(l)%name)//"_OLD")
                told%arr = t%arr
            ENDIF

            ! Map scalar field to device for new timestep, which is required by tstsca4 and bound_sca
            ! Use t_offload as it is a pointer to the t%arr
            !$omp target update to(t_offload)
    
            ! TSTSCA4 zeroize qtu, qtv, qtw before use internally
            CALL tstsca4(scalar(l))

            ! ┌────────────────────────────────────────────────────────────────────────────┐
            ! | SIMPLIFICATION: Only grids on the same level                               |
            ! | *Removed minlevel and maxlevel*                                            |
            ! └────────────────────────────────────────────────────────────────────────────┘
            ! This operation apply boundary conditions to qtu_offload, qtv_offload, qtw_offload ONLY!
            ! Does not modify t-field at all!
            DO ilevel = ISCA_LEVEL, ISCA_LEVEL
                ! Simplification only same grid level: CALL parent(ilevel, qtu, qtv, qtw)
                CALL bound_sca(ilevel, t)
            END DO

            ! Remove comments for verifications of values
            !!$omp target update from(qtu_offload, qtv_offload, qtw_offload)
            !print *, "---------------"
            !print *, MAXVAL(qtu_offload)
            !print *, MAXVAL(qtv_offload)
            !print *, MAXVAL(qtw_offload)
    
            ! fluxbalance zeroize qtt before use internally
            CALL fluxbalance(qtt)

            ! Ghost cell "flux" boundary condition applied to qtt field
            IF (ib%type == "GHOSTCELL") THEN
                CALL set_scastencils("P", scalar(l), qtt=qtt)
            END IF
    
            ! In IRK 1, FRHS is zero, therefore we do not need to zeroize
            ! the dt field before each step
            CALL rkscheme%get_coeffs(frhs, fu, dtrk, dtrki, irk)
    
            ! dT_j = A_j*dT_(j-1) + QTT
            ! T_j = T_(j-1) + B_j*dT_j
            CALL rkstep(t%arr, dt_f%arr, qtt%arr, frhs, dt*fu)
    
            ! Mask blocked cells
            CALL maskbt(t)
    
            ! Ghost cell "value" boundary condition applied to t field
            IF (ib%type == "GHOSTCELL") THEN
                CALL connect(layers=2, s1=t, corners=.TRUE.)
                CALL set_scastencils("P", scalar(l), t=t)
            END IF
    
            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, t%arr, t%arr, 'T')
            END DO
    
            CALL connect(layers=2, s1=t, corners=.TRUE.)

            ! TODO: Fill ghost layers of T (maybe only at last IRK?)
        END DO
        
        CALL stop_timer(402)

        CALL qtt%finish()
        
        CALL stop_timer(400)
    END SUBROUTINE timeintegrate_scalar


    SUBROUTINE itinfo_scalar(itstep, ittot, timeph, dt, exploded)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(inout) :: exploded

        ! Local variables
        INTEGER(intk) :: i, l, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f

        REAL(realk), POINTER, CONTIGUOUS :: t(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        REAL(realk) :: tmean(nsca), tmeansqr(nsca)

        IF (.NOT. solve_scalar) RETURN
        CALL start_timer(400)
        CALL start_timer(420)

        ! Get fields
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        ! Compute CFL, divergence and kinetic energy
        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            DO l = 1, nsca
                CALL get_fieldptr(t, scalar(l)%name, igrid)
                CALL comp_tmean(tmean(l), tmeansqr(l), kk, jj, ii, t, ddx, &
                    ddy, ddz)
            END DO

            CALL itinfo_scalar_sample(igrid, tmean, tmeansqr)
        END DO

        CALL itinfo_scalar_print(itstep, ittot, timeph, exploded)

        CALL stop_timer(420)
        CALL stop_timer(400)
    END SUBROUTINE itinfo_scalar


    SUBROUTINE tstsca4(sca)
        ! Subroutine arguments
        TYPE(scalar_t), INTENT(in) :: sca

        ! Local variables
        INTEGER(intk) :: igrid

        REAL(realk) :: prmol_offload, prturb_offload
        INTEGER(intk) :: kayscrawford_offload
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Expensive data to offload

        CALL start_timer(410)

        !$omp target update to(g_offload)

        prmol_offload = sca%prmol
        prturb_offload = prturb
        kayscrawford_offload = sca%kayscrawford
        
        DO igrid = 1, nmygrids

            CALL get_mgbasb_target(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL tstsca4_grid(kk, jj, ii, &
                prmol_offload, kayscrawford_offload, prturb_offload, &
                nfro, nbac, nrgt, nlft, nbot, ntop, ilesmodel, &
                gmol, rho, igrid)
        END DO

        CALL stop_timer(410)
    END SUBROUTINE tstsca4

    !> @brief Kernel function to timeintegrate a scalar field for a specific grid
    !! 
    !! This whole function shall be executed on the target device.
    !! It is executed on a target-device per-grid basis, meaning one grid is processed at a time.
    !!
    !! @param[in] kk, jj, ii    Grid-specific dimensions in the z, y, and x directions
    !! @param[in] sca_prmol      Scalar field specific value of prmol
    !! @param[in] sca_kayscrawford Scalar field specific value of kayscrawford
    !! @param[in] sca_prturb    Scalar field-specific turbulent Prandtl number
    !! @param[in] nfro, nbac, nrgt, nlft, nbot, ntop Number of boundary cells in each direction
    !! @param[in] ilesmodel_offlad Encoded index of the Large Eddy Simulation (LES) model
    !! @param[in] gmol_offload  Flow parameter gmol on target device
    !! @param[in] rho_offload   Flow parameter rho on target device
    !! @param[in] igrid         Grid index being processed
    SUBROUTINE tstsca4_grid(kk, jj, ii, &
            sca_prmol, sca_kayscrawford, sca_prturb, nfro, nbac, nrgt, nlft, &
            nbot, ntop, ilesmodel_offlad, gmol_offload, rho_offload, igrid)  
        
            ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii, igrid
        REAL(realk), INTENT(in) :: sca_prmol, sca_prturb
        INTEGER(intk), INTENT(IN) :: sca_kayscrawford
        INTEGER(intk), INTENT(IN) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk), INTENT(IN) :: ilesmodel_offlad
        REAL(realk), INTENT(IN) :: gmol_offload, rho_offload

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: nfu, nbu, nrv, nlv, nbw, ntw
        INTEGER(intk) :: iles
        REAL(realk) :: gsca(ii)
        REAL(realk) :: adv, diff, area
        REAL(realk) :: gscamol, gtgmolp, gtgmoln

        ! Usually, the fluxes across the grid boundaries are already set
        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

        ! Only for CON boundaries, fluxes are computed for one more layer
        ! (this avoids a connect on qtu, qtv, qtw)
        IF (nfro == 7) nfu = 1
        IF (nbac == 7) nbu = 1
        IF (nrgt == 7) nrv = 1
        IF (nlft == 7) nlv = 1
        IF (nbot == 7) nbw = 1
        IF (ntop == 7) ntw = 1

        iles = 1
        IF (ilesmodel_offlad == 0) iles = 0

        !! Create the parallel region for per-grid parallelization
        !! All GPU threads work on a single grid simultaneously
        !! The execution moves sequentially from one grid to another.
        !! Repeat this for the X, Y and Z direction
        
        ! X direction
        !$omp target data map(to: gsca)
        !$omp target
        BLOCK
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: g
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: t, u, bt, qtu
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddy, ddz, rdx

            ! Get grid-specific data required for the computation
            CALL ptr_to_grid3(qtu_offload, igrid, qtu)
            CALL ptr_to_grid3(g_offload, igrid, g)
            CALL ptr_to_grid3(t_offload, igrid, t)
            CALL ptr_to_grid3(u_offload, igrid, u)
            CALL ptr_to_grid3(bt_offload, igrid, bt)
            CALL ptr_to_grid_y(ddy_offload, igrid, ddy)
            CALL ptr_to_grid_z(ddz_offload, igrid, ddz)
            CALL ptr_to_grid_x(rdx_offload, igrid, rdx)

            !$omp teams distribute parallel do collapse(2)
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    ! Scalar diffusivity LES/DNS computation
                    IF (iles == 1) THEN
                        DO k = 3, kk-2
                            gscamol = gmol_offload/rho_offload/sca_prmol
                            gtgmolp = (g(k, j, i) - gmol_offload)/gmol_offload
                            gtgmoln = (g(k, j, i+1) - gmol_offload)/gmol_offload

                            ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                            gsca(k) = gscamol &
                            + (g(k, j, i+1) + g(k, j, i) - 2.0*gmol_offload) / rho_offload &
                            / (sca_prt(sca_prmol, sca_kayscrawford, sca_prturb, gtgmoln) + sca_prt(sca_prmol, sca_kayscrawford, sca_prturb, gtgmolp))

                            ! Limit gsca here MAX(...,0): no negative diffusion!
                            gsca(k) = MAX(gscamol, gsca(k))
                        END DO
                    ELSE
                        !$omp simd
                        DO k = 3, kk-2
                            gsca(k) = gmol_offload/rho_offload/sca_prmol
                        END DO
                        !$omp end simd
                    END IF

                    ! Final asembly
                    !$omp simd
                    DO k = 3, kk-2
                        ! Convective fluxes
                        ! It is assumed that the velocity field is already masked
                        ! with BU, BV, BW = no new masking necessary (!)
                        adv = (ddy(j)*ddz(k)) * u(k, j, i) &
                            * 0.5 * (t(k, j, i) + t(k, j, i+1))
                        ! Depending on the knowledge about the the cell and its
                        ! neighbours it is determined if faces are blocked (=0)
                        ! or open (=1)
                        area = bt(k, j, i)*bt(k, j, i+1)*(ddy(j)*ddz(k))
                        diff = -gsca(k)*rdx(i)*(t(k, j, i+1) - t(k, j, i))*area                  
                        ! Final result
                        !print *, adv
                        qtu(k, j, i) = adv + diff
                    END DO
                    !$omp end simd
                END DO
            END DO
            !$omp end teams distribute parallel do
        END BLOCK
        !$omp end target
        !$omp end target data

        ! Y direction
        !$omp target data map(to: gsca)
        !$omp target
        BLOCK
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: g
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: t, v, bt, qtv
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx, ddz, rdy

            CALL ptr_to_grid3(qtv_offload, igrid, qtv)
            CALL ptr_to_grid3(g_offload, igrid, g)
            CALL ptr_to_grid3(t_offload, igrid, t)
            CALL ptr_to_grid3(v_offload, igrid, v)
            CALL ptr_to_grid3(bt_offload, igrid, bt)
            CALL ptr_to_grid_x(ddx_offload, igrid, ddx)
            CALL ptr_to_grid_z(ddz_offload, igrid, ddz)
            CALL ptr_to_grid_y(rdy_offload, igrid, rdy)

            !$omp teams distribute parallel do collapse(2)        
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    ! Scalar diffusivity LES/DNS computation
                    IF (iles == 1) THEN
                        DO k = 3, kk-2
                            gscamol = gmol_offload/rho_offload/sca_prmol
                            gtgmolp = (g(k, j, i) - gmol_offload)/gmol_offload
                            gtgmoln = (g(k, j+1, i) - gmol_offload)/gmol_offload

                            ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                            gsca(k) = gscamol &
                                + (g(k, j+1, i) + g(k, j, i) - 2.0*gmol_offload) / rho_offload &
                                / (sca_prt(sca_prmol, sca_kayscrawford, sca_prturb, gtgmoln) + sca_prt(sca_prmol, sca_kayscrawford, sca_prturb, gtgmolp))

                            ! Limit gsca here MAX(...,0): no negative diffusion!
                            gsca(k) = MAX(gscamol, gsca(k))
                        END DO
                    ELSE
                        !$omp simd
                        DO k = 3, kk-2
                            gsca(k) = gmol_offload/rho_offload/sca_prmol
                        END DO
                        !$omp end simd
                    END IF

                    ! Final asembly
                    !$omp simd
                    DO k = 3, kk-2
                        ! Convective fluxes
                        ! It is assumed that the velocity field is already masked
                        ! with BU, BV, BW = no new masking necessary (!)
                        adv = (ddx(i)*ddz(k)) * v(k, j, i) &
                            * 0.5 * (t(k, j, i) + t(k, j+1, i))

                        ! Depending on the knowledge about the the cell and its
                        ! neighbours it is determined if faces are blocked (=0)
                        ! or open (=1)
                        area = bt(k, j, i)*bt(k, j+1, i)*(ddx(i)*ddz(k))
                        diff = -gsca(k)*rdy(j)*(t(k, j+1, i) - t(k, j, i))*area

                        ! Final result
                        qtv(k, j, i) = adv + diff
                    END DO
                    !$omp end simd
                END DO
            END DO
            !$omp end teams distribute parallel do
        END BLOCK
        !$omp end target
        !$omp end target data

        ! Z direction
        !$omp target data map(to: gsca)
        !$omp target
        BLOCK
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: g
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: t, w, bt, qtw
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx, ddy, rdy

            CALL ptr_to_grid3(qtw_offload, igrid, qtw)
            CALL ptr_to_grid3(g_offload, igrid, g)
            CALL ptr_to_grid3(t_offload, igrid, t)
            CALL ptr_to_grid3(w_offload, igrid, w)
            CALL ptr_to_grid3(bt_offload, igrid, bt)
            CALL ptr_to_grid_x(ddx_offload, igrid, ddx)
            CALL ptr_to_grid_y(ddy_offload, igrid, ddy)
            CALL ptr_to_grid_y(rdy_offload, igrid, rdy)
            
            !$omp teams distribute parallel do collapse(2)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    ! Scalar diffusivity LES/DNS computation
                    IF (iles == 1) THEN
                        DO k = 3-nbw, kk-3+ntw
                            gscamol = gmol_offload/rho_offload/sca_prmol
                            gtgmolp = (g(k, j, i) - gmol_offload)/gmol_offload
                            gtgmoln = (g(k+1, j, i) - gmol_offload)/gmol_offload

                            ! 1/Re * 1/Pr + 1/Re_t * 1/Pr_t:
                            gsca(k) = gscamol &
                                + (g(k+1, j, i) + g(k, j, i) - 2.0*gmol_offload) / rho_offload &
                                / (sca_prt(sca_prmol, sca_kayscrawford, sca_prturb, gtgmoln) + sca_prt(sca_prmol, sca_kayscrawford, sca_prturb, gtgmolp))

                            ! Limit gsca here MAX(...,0): no negative diffusion!
                            gsca(k) = MAX(gscamol, gsca(k))
                        END DO
                    ELSE
                        !$omp simd
                        DO k = 3-nbw, kk-3+ntw
                            gsca(k) = gmol_offload/rho_offload/sca_prmol
                        END DO
                        !$omp end simd
                    END IF

                    ! Final asembly
                    !$omp simd
                    DO k = 3-nbw, kk-3+ntw
                        ! Convective fluxes
                        ! It is assumed that the velocity field is already masked
                        ! with BU, BV, BW = no new masking necessary (!)
                        adv = (ddx(i)*ddy(j)) * w(k, j, i) &
                            * 0.5 * (t(k, j, i) + t(k+1, j, i))

                        ! Depending on the knowledge about the the cell and its
                        ! neighbours it is determined if faces are blocked (=0)
                        ! or open (=1)
                        area = bt(k, j, i)*bt(k+1, j, i)*(ddx(i)*ddy(j))
                        diff = -gsca(k)*rdy(j)*(t(k+1, j, i) - t(k, j, i))*area

                        ! Final result
                        qtw(k, j, i) = adv + diff
                    END DO
                    !$omp end simd
                END DO
            END DO
            !$omp end teams distribute parallel do
        END BLOCK
        !$omp end target
        !$omp end target data

    END SUBROUTINE tstsca4_grid

    !> @brief Calculation of local turbulent Prandtl number
    !!
    !! Target device ready implementation mirroring the prt Function in scacore_mod.
    !!
    !! @param[in] sca_prmol        Scalar-field specific prmol
    !! @param[in] sca_kayscrawford Scalar-field specific kayscrawford
    !! @param[in] sca_prturb       Scalar-field specific prturb
    !! @return     Local turbulent Prandtl number
    FUNCTION sca_prt(sca_prmol, sca_kayscrawford, sca_prturb, sca_gtgmol) RESULT(res)
        !$omp declare target
        REAL(realk), INTENT(IN) :: sca_prmol, sca_prturb, sca_gtgmol
        INTEGER(intk), INTENT(IN) :: sca_kayscrawford

        REAL(realk) :: res
        REAL(realk) :: kayscrawford

        IF (sca_kayscrawford == 0) THEN
            res = sca_prturb
        ELSE
            IF (sca_gtgmol > 0.0) THEN
                kayscrawford = 0.5882 + 0.228*sca_gtgmol &
                    - 0.0441*sca_gtgmol**2*(1.0 - exp(-5.165/sca_gtgmol))
            ELSE
                kayscrawford = sca_prmol
            ENDIF
            res = kayscrawford
        END IF
    END FUNCTION sca_prt

    SUBROUTINE fluxbalance(qtt_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: qtt_f

        ! Local variables
        INTEGER(intk) :: igrid
        INTEGER(intk) :: kk, jj, ii

        CALL start_timer(411)

        DO igrid = 1, nmygrids
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL fluxbalance_grid(kk, jj, ii, igrid)
        END DO

        !$omp target update from(qtt_offload)
        qtt_f%arr = qtt_offload

        CALL stop_timer(411)
    END SUBROUTINE fluxbalance

    !> @brief Kernel function for fluxbalance on a specific grid.
    !!
    !! This subroutine computes the net flux resulting from neighboring grid cells.
    !!  This whole function shall be executed on the target device.
    !!
    !! @param[in] kk, jj, ii  Grid-specific dimensions in the z, y, and x directions
    !! @param[in] igrid        Grid index being processed
    SUBROUTINE fluxbalance_grid(kk, jj, ii, igrid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: kk, jj, ii, igrid

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: netflux

        ! Set INTENT(out) to zero
        !qtt = 0.0

        !$omp target
        BLOCK
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: qtu, qtv, qtw, qtt
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: rddx, rddy, rddz
            CALL ptr_to_grid3(qtu_offload, igrid, qtu)
            CALL ptr_to_grid3(qtv_offload, igrid, qtv)
            CALL ptr_to_grid3(qtw_offload, igrid, qtw)
            CALL ptr_to_grid3(qtt_offload, igrid, qtt)

            CALL ptr_to_grid_x(rddx_offload, igrid, rddx)
            CALL ptr_to_grid_y(rddy_offload, igrid, rddy)
            CALL ptr_to_grid_z(rddz_offload, igrid, rddz)
            !$omp teams distribute parallel do collapse(2)
            DO i = 3, ii-2
                DO j = 3, jj-2
                    !$omp simd
                    DO k = 3, kk-2
                        ! Computing netflux resulting from exchange with neighbors
                        netflux = qtu(k, j, i-1) - qtu(k, j, i) + qtv(k, j-1, i) &
                            - qtv(k, j, i) + qtw(k-1, j, i) - qtw(k, j, i)

                        qtt(k, j, i) = rddz(k)*rddy(j)*rddx(i)*netflux
                    END DO
                    !$omp end simd
                END DO
            END DO
            !$omp end teams distribute parallel do
        END BLOCK
        !$omp end target
    END SUBROUTINE fluxbalance_grid

    SUBROUTINE comp_tmean(tmean, tmeansqr, kk, jj, ii, t, ddx, ddy, ddz)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: tmean, tmeansqr
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: t(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        REAL(realk) :: vsum, vol
        INTEGER(intk) :: k, j, i

        tmean = 0.0
        tmeansqr = 0.0
        vsum = 0.0

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    vol = ddx(i)*ddy(j)*ddz(k)
                    tmean = tmean + vol*t(k, j, i)
                    tmeansqr = tmeansqr + vol*t(k, j, i)**2
                    vsum = vsum + vol
                END DO
            END DO
        END DO

        tmean = tmean/vsum
        tmeansqr = tmeansqr/vsum
    END SUBROUTINE comp_tmean
END MODULE timeintegrate_scalar_mod
