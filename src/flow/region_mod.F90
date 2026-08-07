MODULE region_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: region_type_box = 1

    TYPE :: region_t
        CHARACTER(len=nchar_name) :: name
        INTEGER(intk) :: id
        INTEGER(intk) :: itype
        REAL(realk) :: box_min(3)
        REAL(realk) :: box_max(3)
    END TYPE region_t

    LOGICAL, PROTECTED :: has_regions = .FALSE.
    INTEGER(intk), PROTECTED :: nregions = 0
    TYPE(region_t), ALLOCATABLE, PROTECTED :: regions(:)

    ! Per-cell region-id field, 0 = not part of any region. Built once at
    ! init from static grid geometry - not exchanged/re-built per timestep,
    ! and not registered in the field_t "get_field" registry (that registry
    ! is field_t-only); consumers 'USE region_mod, ONLY: regionid' directly.
    TYPE(intfield_t), PROTECTED, TARGET :: regionid

    PUBLIC :: has_regions, regionid, init_region, finish_region, get_region_id

CONTAINS

    SUBROUTINE init_region()
        ! Subroutine arguments
        ! none...

        ! Local variables
        TYPE(config_t) :: rc
        CHARACTER(len=64) :: jsonptr
        CHARACTER(len=nchar_name) :: ctype
        INTEGER(intk) :: n, m, i, igrid, kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: x(:), y(:), z(:)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        INTEGER(ifk), POINTER, CONTIGUOUS :: rid(:, :, :)

        has_regions = .FALSE.
        nregions = 0
        IF (.NOT. fort7%exists("/regions")) THEN
            RETURN
        END IF
        has_regions = .TRUE.

        CALL fort7%get_size("/regions", nregions)
        ALLOCATE(regions(nregions))

        DO n = 1, nregions
            WRITE(jsonptr, '("/regions/", I0)') n-1
            CALL fort7%get(rc, jsonptr)

            CALL rc%get_value("/name", regions(n)%name)
            CALL rc%get_value("/id", regions(n)%id)
            CALL rc%get_value("/type", ctype)

            SELECT CASE (lower(TRIM(ctype)))
            CASE ("box")
                regions(n)%itype = region_type_box
                CALL rc%get_array("/min", regions(n)%box_min)
                CALL rc%get_array("/max", regions(n)%box_max)
            CASE DEFAULT
                WRITE(*, '(A)') "Invalid region type: "//TRIM(ctype)
                CALL errr(__FILE__, __LINE__)
            END SELECT

            CALL rc%finish()
        END DO

        ! Validate unique id's and names before building the field - a
        ! parse-time mistake should abort here, not surface as a silent
        ! mis-classified cell later
        DO n = 1, nregions
            DO m = n+1, nregions
                IF (regions(n)%id == regions(m)%id) THEN
                    WRITE(*, '(A, I0)') "Duplicate region id: ", regions(n)%id
                    CALL errr(__FILE__, __LINE__)
                END IF
                IF (TRIM(regions(n)%name) == TRIM(regions(m)%name)) THEN
                    WRITE(*, '(A)') "Duplicate region name: "// &
                        TRIM(regions(n)%name)
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO
        END DO

        CALL regionid%init("REGIONID")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(x, "X", igrid)
            CALL get_fieldptr(y, "Y", igrid)
            CALL get_fieldptr(z, "Z", igrid)
            CALL get_fieldptr(bp, "BP", igrid)
            CALL regionid%get_ptr(rid, igrid)

            CALL build_region_grid(kk, jj, ii, x, y, z, bp, rid)
        END DO

        IF (myid == 0) THEN
            WRITE(*, '("REGIONS:")')
            DO n = 1, nregions
                WRITE(*, '(2X, "Name: ", A, ", id: ", I0)') &
                    TRIM(regions(n)%name), regions(n)%id
            END DO
            WRITE(*, '()')
        END IF
    END SUBROUTINE init_region


    SUBROUTINE finish_region()
        IF (has_regions) THEN
            CALL regionid%finish()
        END IF
    END SUBROUTINE finish_region


    ! Resolve a region name to its id. Returns -1 if the name is not
    ! among the defined regions (caller's responsibility to error out -
    ! this module doesn't know whether an unresolved name is fatal for
    ! the caller).
    INTEGER(intk) FUNCTION get_region_id(name) RESULT(id)
        CHARACTER(len=*), INTENT(in) :: name

        INTEGER(intk) :: n

        id = -1
        DO n = 1, nregions
            IF (TRIM(regions(n)%name) == TRIM(name)) THEN
                id = regions(n)%id
                RETURN
            END IF
        END DO
    END FUNCTION get_region_id


    SUBROUTINE build_region_grid(kk, jj, ii, x, y, z, bp, rid)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)
        REAL(realk), INTENT(in) :: bp(kk, jj, ii)
        INTEGER(ifk), INTENT(out) :: rid(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k, n
        INTEGER(intk) :: noverlap

        rid = 0_ifk

        ! Cell-center coordinate box test. This is a purely local,
        ! per-cell test (no inter-level or inter-rank dependency, unlike
        ! IB's seeded flood-fill), so every cell - including this grid's
        ! own ghost/halo layers - is classified directly from its own
        ! (locally known) coordinates. No connect()/halo exchange is
        ! needed: two neighboring ranks independently compute identical
        ! classifications for what is physically the same location.
        DO n = 1, nregions
            IF (regions(n)%itype /= region_type_box) CYCLE

            noverlap = 0
            DO i = 1, ii
                IF (x(i) < regions(n)%box_min(1)) CYCLE
                IF (x(i) > regions(n)%box_max(1)) CYCLE
                DO j = 1, jj
                    IF (y(j) < regions(n)%box_min(2)) CYCLE
                    IF (y(j) > regions(n)%box_max(2)) CYCLE
                    DO k = 1, kk
                        IF (z(k) < regions(n)%box_min(3)) CYCLE
                        IF (z(k) > regions(n)%box_max(3)) CYCLE

                        ! Region-vs-region overlap is a hard error -
                        ! genuinely ambiguous which region's coefficients
                        ! would apply
                        IF (rid(k, j, i) /= 0) THEN
                            WRITE(*, '(A, A)') &
                                "Region overlaps another region: ", &
                                TRIM(regions(n)%name)
                            CALL errr(__FILE__, __LINE__)
                        END IF
                        rid(k, j, i) = INT(regions(n)%id, ifk)

                        ! Region-vs-IB overlap is informational only - a
                        ! Darcy-type region is a homogenized continuum
                        ! (the permeability already lumps whatever solid
                        ! fraction exists in the calibration volume), so
                        ! there's no physical reason to forbid it; this
                        ! only guards against two independently-defined
                        ! features accidentally colliding by mistake.
                        IF (bp(k, j, i) < 0.5) THEN
                            noverlap = noverlap + 1
                        END IF
                    END DO
                END DO
            END DO

            IF (noverlap > 0) THEN
                WRITE(*, '(A, A, A, I0, A)') "Warning: region ", &
                    TRIM(regions(n)%name), " overlaps IB-blocked geometry in ", &
                    noverlap, " cell(s)"
            END IF
        END DO
    END SUBROUTINE build_region_grid
END MODULE region_mod