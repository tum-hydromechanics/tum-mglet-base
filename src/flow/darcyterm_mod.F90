MODULE darcyterm_mod
    USE core_mod
    USE flowcore_mod, ONLY: gmol
    USE region_mod, ONLY: has_regions, regionid, get_region_id

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE :: darcyregion_t
        INTEGER(intk) :: id
        REAL(realk) :: cdarcy
    END TYPE darcyregion_t

    LOGICAL, PROTECTED :: has_darcy = .FALSE.
    INTEGER(intk), PROTECTED :: ndarcyregions = 0
    TYPE(darcyregion_t), ALLOCATABLE, PROTECTED :: darcyregions(:)

    PUBLIC :: init_darcyterm, finish_darcyterm, darcyterm

CONTAINS

    SUBROUTINE init_darcyterm()
        ! Subroutine arguments
        ! none...

        ! Local variables
        TYPE(config_t) :: darcyconf, rc
        CHARACTER(len=64) :: jsonptr
        CHARACTER(len=nchar_name) :: regionname
        INTEGER(intk) :: n, id

        has_darcy = .FALSE.
        IF (.NOT. fort7%exists("/flow/darcy")) THEN
            RETURN
        END IF

        IF (.NOT. has_regions) THEN
            WRITE(*, '(A)') "flow/darcy given, but no /regions are defined"
            CALL errr(__FILE__, __LINE__)
        END IF
        has_darcy = .TRUE.

        CALL fort7%get(darcyconf, "/flow/darcy")
        CALL darcyconf%get_size("/regions", ndarcyregions)
        ALLOCATE(darcyregions(ndarcyregions))

        DO n = 1, ndarcyregions
            WRITE(jsonptr, '("/regions/", I0)') n-1
            CALL darcyconf%get(rc, jsonptr)

            CALL rc%get_value("/region", regionname)
            id = get_region_id(regionname)
            IF (id < 0) THEN
                WRITE(*, '(A)') "flow/darcy: unknown region: "// &
                    TRIM(regionname)
                CALL errr(__FILE__, __LINE__)
            END IF
            darcyregions(n)%id = id

            CALL rc%get_value("/cdarcy", darcyregions(n)%cdarcy)

            CALL rc%finish()
        END DO
        CALL darcyconf%finish()

        IF (myid == 0) THEN
            WRITE(*, '("DARCY TERM:")')
            DO n = 1, ndarcyregions
                WRITE(*, '(2X, "Region id ", I0, ": cdarcy = ", G0)') &
                    darcyregions(n)%id, darcyregions(n)%cdarcy
            END DO
            WRITE(*, '()')
        END IF
    END SUBROUTINE init_darcyterm

    
    SUBROUTINE finish_darcyterm()
        ! placeholder 
        CONTINUE
    END SUBROUTINE finish_darcyterm


    ! Applies the isotropic Darcy resistance as an operator-split
    ! correction to the already RK-updated velocity - NOT an addition to
    ! the explicit RHS accumulator (contrast with boussinesqterm/
    ! coriolisterm, which modify uo/vo/wo before rkstep). Caller is
    ! responsible for calling this after the per-substep halo-exchange/BC
    ! block and before the pressure projection - see timeintegration_mod.
    SUBROUTINE darcyterm(u_f, v_f, w_f, dtrki)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f, v_f, w_f
        REAL(realk), INTENT(in) :: dtrki

        ! Local variables
        INTEGER(intk) :: i, n, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        INTEGER(ifk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: rid

        IF (.NOT. has_darcy) RETURN
        CALL start_timer(380)

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL regionid%get_ptr(rid, igrid)

            ! Legacy applies one closed-form correction per configured
            ! Darcy-active region, sequentially - mirrored here rather
            ! than trying to blend multiple regions' coefficients into a
            ! single pass, which would be ambiguous at a boundary between
            ! two differently-parameterized regions
            DO n = 1, ndarcyregions
                CALL darcyterm_grid(kk, jj, ii, u, v, w, rid, &
                    darcyregions(n)%id, darcyregions(n)%cdarcy, dtrki, &
                    nfro, nbac, nrgt, nlft, nbot, ntop)
            END DO
        END DO

        CALL stop_timer(380)
    END SUBROUTINE darcyterm


    SUBROUTINE darcyterm_grid(kk, jj, ii, u, v, w, rid, regid, cdarcy, &
            dtrki, nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii), &
            w(kk, jj, ii)
        INTEGER(ifk), INTENT(in) :: rid(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: regid
        REAL(realk), INTENT(in) :: cdarcy, dtrki
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        INTEGER(ifk) :: regid_ifk
        REAL(realk) :: fak

        regid_ifk = INT(regid, ifk)

        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

        ! CON = 7
        IF (nbac == 7) nbu = 1
        IF (nlft == 7) nlv = 1
        IF (ntop == 7) ntw = 1

        ! OP1 = 3
        IF (nfro == 3) nfu = 1
        IF (nbac == 3) nbu = 1
        IF (nrgt == 3) nrv = 1
        IF (nlft == 3) nlv = 1
        IF (nbot == 3) nbw = 1
        IF (ntop == 3) ntw = 1

        ! U-velocity - blend across the face-normal (I-direction)
        ! neighbors only. Safe under the per-substep halo exchange
        ! (timeintegration_mod.F90's CALL connect(..., normal=.true.)
        ! only refreshes the face-normal component - a transverse-neighbor
        ! blend would not be, see Darcy-Forchheimer/NOTES.md).
        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    fak = 0.5*(ind(rid(k, j, i)) + ind(rid(k, j, i+1)))
                    IF (fak > 0.0) THEN
                        u(k, j, i) = u(k, j, i) &
                            /(1.0 + fak*cdarcy*gmol*dtrki)
                    END IF
                END DO
            END DO
        END DO

        ! V-velocity - blend across the face-normal (J-direction)
        ! neighbors only
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    fak = 0.5*(ind(rid(k, j, i)) + ind(rid(k, j+1, i)))
                    IF (fak > 0.0) THEN
                        v(k, j, i) = v(k, j, i) &
                            /(1.0 + fak*cdarcy*gmol*dtrki)
                    END IF
                END DO
            END DO
        END DO

        ! W-velocity - blend across the face-normal (K-direction)
        ! neighbors only
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    fak = 0.5*(ind(rid(k, j, i)) + ind(rid(k+1, j, i)))
                    IF (fak > 0.0) THEN
                        w(k, j, i) = w(k, j, i) &
                            /(1.0 + fak*cdarcy*gmol*dtrki)
                    END IF
                END DO
            END DO
        END DO

    CONTAINS
        PURE REAL(realk) FUNCTION ind(id)
            INTEGER(ifk), INTENT(in) :: id
            IF (id == regid_ifk) THEN
                ind = 1.0
            ELSE
                ind = 0.0
            END IF
        END FUNCTION ind
    END SUBROUTINE darcyterm_grid
END MODULE darcyterm_mod