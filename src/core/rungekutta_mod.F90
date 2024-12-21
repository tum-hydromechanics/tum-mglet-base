MODULE rungekutta_mod
    USE charfunc_mod, ONLY: lower
    USE err_mod, ONLY: errr
    USE precision_mod, ONLY: intk, realk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: maxnrk = 6
    REAL(realk), PARAMETER :: oneeps = (1.0_realk - EPSILON(1.0_realk))

    TYPE, ABSTRACT :: rk_t
        INTEGER(intk) :: nrk = 0
        REAL(realk) :: cflmax = 0.0
        REAL(realk) :: c(maxnrk) = 0.0
    CONTAINS
        PROCEDURE(init_i), DEFERRED :: init
    END TYPE rk_t

    ABSTRACT INTERFACE
        SUBROUTINE init_i(this, ctyp)
            IMPORT :: rk_t
            CLASS(rk_t), INTENT(out) :: this
            CHARACTER(len=*), INTENT(in) :: ctyp
        END SUBROUTINE init_i
    END INTERFACE

    ! 2N storage versions
    ! A more logical name would be 2n_rk_t but it is forbidden to have a
    ! variable name that starts with a digit
    TYPE, EXTENDS(rk_t) :: rk_2n_t
        REAL(realk) :: a(maxnrk) = 0.0
        REAL(realk) :: b(maxnrk) = 0.0
    CONTAINS
        PROCEDURE :: init => init_2n
        GENERIC :: get_coeffs => get_coeffs_a, get_coeffs_b

        PROCEDURE, PRIVATE :: get_coeffs_a, get_coeffs_b
        PROCEDURE, PRIVATE :: init_eeuler
        PROCEDURE, PRIVATE :: init_williamson
        PROCEDURE, PRIVATE :: init_berland
        PROCEDURE, PRIVATE :: init_carpenter
        PROCEDURE, PRIVATE :: init_bernardini
        PROCEDURE, PRIVATE :: comp_c
    END TYPE rk_2n_t

    ! rk type for particle timeintegration
    TYPE, EXTENDS(rk_t) :: bt_rk_t

        ! coefficients represent the common butcher tableau, where
        ! c: intermediate timesteps
        ! a: intermediate weights
        ! b: weights of the final stage

        REAL(realk) :: a(maxnrk, maxnrk) = 0.0
        REAL(realk) :: b(maxnrk) = 0.0

    CONTAINS
        PROCEDURE :: init => init_bt_rk
        PROCEDURE :: get_coeffs => get_bt_coefficients
        PROCEDURE, PRIVATE :: init_euler_bt
        PROCEDURE, PRIVATE :: init_williamson_bt
    END TYPE bt_rk_t

    PUBLIC :: rk_2n_t, bt_rk_t, rkstep, prkstep
CONTAINS

    SUBROUTINE init_2n(this, ctyp)
        ! Function arguments
        CLASS(rk_2n_t), INTENT(out) :: this
        CHARACTER(len=*), INTENT(in) :: ctyp

        SELECT CASE(lower(TRIM(ctyp)))
        CASE("euler")
            CALL this%init_eeuler()
        CASE ("williamson")
            CALL this%init_williamson()
        CASE ("berland")
            CALL this%init_berland()
        CASE ("carpenter")
            CALL this%init_carpenter()
        CASE ("bernardini")
            CALL this%init_bernardini()
        CASE DEFAULT
            WRITE(*, *) "Invalid: ", lower(TRIM(ctyp))
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE init_2n

    SUBROUTINE init_eeuler(rk)
        ! Classic explicit Euler Scheme

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%c = 0.0
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 1
        ! TODO: cflmax
        !rk%cflmax = oneeps*SQRT(3.0)
        rk%c(1:rk%nrk) = [0.0]
        rk%a(1:rk%nrk) = [0.0]
        rk%b(1:rk%nrk) = [1.0]

        ! Compute C from A and B
        !CALL rk%comp_c()
    END SUBROUTINE init_eeuler

    SUBROUTINE init_williamson(rk)
        ! Classic Williamson RK3 scheme
        ! From: Williamson, John H, Low-storage Runge-Kutta schemes. Journal
        ! of computational physics 35.1 (1980)
        ! DOI: https://doi.org/10.1016/0021-9991(80)90033-9

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%c = 0.0
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 3
        rk%cflmax = oneeps*SQRT(3.0)
        rk%c(1:rk%nrk) = [0.0, 1.0/3.0, 3.0/4.0]
        rk%a(1:rk%nrk) = [0.0, -5.0/9.0, -153.0/128.0]
        rk%b(1:rk%nrk) = [1.0/3.0, 15.0/16.0, 8.0/15.0]

        ! Compute C from A and B
        ! CALL rk%comp_c()
    END SUBROUTINE init_williamson


    SUBROUTINE init_berland(rk)
        ! 4th order 6-stage low disipation scheme
        ! From: Berland, Bogey and Bailly, Low-dissipation and low-dispersion
        ! fourth-order Runge–Kutta algorithm, Computers & Fluids 35 (2006),
        ! DOI: 10.1016/j.compfluid.2005.04.003

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 6
        rk%cflmax = 3.80 - 0.01  ! 3.80 from paper, 0.01 round off margin
        rk%a(1:rk%nrk) = [0.0, -0.737101392796, -1.634740794341, &
            -0.744739003780, -1.469897351522, -2.813971388035]
        rk%b(1:rk%nrk) = [0.032918605146, 0.823256998200, 0.381530948900, &
            0.200092213184, 1.718581042715, 0.27]

        ! Compute C from A and B
        CALL rk%comp_c()
    END SUBROUTINE init_berland


    SUBROUTINE init_carpenter(rk)
        ! 4th order 5-stage scheme
        ! Ref: Carpenter & Kennedy, Fourth-Order 2N-Storage Runge-Kutta Schemes,
        ! NASA Technical Memorandum 109112, June 1994.
        !
        ! Online version:
        ! http://www.ece.uvic.ca/~bctill/papers/numacoust/Carpenter_Kennedy_1994.pdf

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 5
        rk%cflmax = 3.34 - 0.01  ! 3.34 given in paper, 0.01 round off margin
        rk%a(1:rk%nrk) = [0.0, -0.4178904745, -1.192151694643, &
            -1.697784692471, -1.514183444257]
        rk%b(1:rk%nrk) = [0.1496590219993, 0.3792103129999, 0.8229550293869, &
            0.6994504559488, 0.1530572479681]

        ! Compute C from A and B
        CALL rk%comp_c()
    END SUBROUTINE init_carpenter


    SUBROUTINE init_bernardini(rk)
        ! 2th order 5-stage scheme
        ! From: Matteo Bernardini and Sergio Pirozzoli, A general strategy for
        ! the optimization of Runge–Kutta schemes for wave propagation
        ! phenomena, Journal of Computational Physics 228 (2009)
        ! DOI: 10.1016/j.jcp.2009.02.032
        !
        ! Coefficients from table 4

        ! Subroutine arguments
        CLASS(rk_2n_t), INTENT(out) :: rk

        ! Set to zero
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        rk%nrk = 5
        rk%cflmax = 3.47 - 0.01  ! 3.47 given in paper, 0.01 round off margin
        rk%a(1:rk%nrk) = [0.0, -1.0, -1.55798, -1.0, -0.45031]
        rk%b(1:rk%nrk) = [0.2, 0.83204, 0.6, 0.35394, 0.2]

        ! Compute C from A and B
        CALL rk%comp_c()
    END SUBROUTINE init_bernardini


    PURE SUBROUTINE get_coeffs_a(this, frhs, fu, dtrk, dtrki, irk)
        CLASS(rk_2n_t), INTENT(in) :: this
        REAL(realk), INTENT(out) :: frhs, fu, dtrk, dtrki
        INTEGER(intk), INTENT(in) :: irk

        IF (irk > this%nrk .OR. irk < 0) ERROR STOP

        CALL this%get_coeffs_b(frhs, fu, irk)

        IF (irk < this%nrk) THEN
            dtrk = this%c(irk+1)
        ELSE
            dtrk = 1.0
        END IF

        dtrki = dtrk - this%c(irk)
    END SUBROUTINE get_coeffs_a


    PURE SUBROUTINE get_coeffs_b(this, frhs, fu, irk)
        CLASS(rk_2n_t), INTENT(in) :: this
        REAL(realk), INTENT(out) :: frhs, fu
        INTEGER(intk), INTENT(in) :: irk

        IF (irk > this%nrk .OR. irk < 0) ERROR STOP

        frhs = this%a(irk)
        fu = this%b(irk)
    END SUBROUTINE get_coeffs_b


    ! Compute coefficient C
    !
    ! Example computation for a six-stage scheme:
    ! C(1) = 0.0
    ! C(2) = B(1)
    ! C(3) = B(1) + B(2)*(A(2) + 1)
    ! C(4) = B(1) + B(2)*(A(2) + 1) + B(3)*(A(3)*(A(2) + 1) + 1)
    ! C(5) = B(1) + B(2)*(A(2) + 1) + B(3)*(A(3)*(A(2) + 1) + 1) + B(4)*(A(4)*(A(3)*(A(2) + 1) + 1) + 1)
    ! C(6) = B(1) + B(2)*(A(2) + 1) + B(3)*(A(3)*(A(2) + 1) + 1) + B(4)*(A(4)*(A(3)*(A(2) + 1) + 1) + 1) + B(5)*(A(5)*(A(4)*(A(3)*(A(2) + 1) + 1) + 1) + 1)
    !
    ! This implementation is generally valid for any number of stages.
    !
    ! Ref: Carpenter & Kennedy, Fourth-Order 2N-Storage Runge-Kutta Schemes,
    ! NASA Technical Memorandum 109112, June 1994.
    !
    ! Online version:
    ! http://www.ece.uvic.ca/~bctill/papers/numacoust/Carpenter_Kennedy_1994.pdf
    SUBROUTINE comp_c(this)
        CLASS(rk_2n_t), INTENT(inout) :: this

        INTEGER(intk) i, j, k
        REAL(realk) :: fak

        ! Make sure C is zero
        this%c = 0.0

        ! Compute C(i)
        DO i = 1, this%nrk
            ! Compute individual terms, i.e. like B(3)*(A(3)*(A(2) + 1) + 1)
            DO j = 1, i-1
                ! Compute factors of individual terms, i.e. like (A(3)*(A(2) + 1) + 1)
                fak = 1.0
                DO k = 2, j
                    fak = this%a(k)*fak + 1.0
                END DO
                this%c(i) = this%c(i) + this%b(j)*fak
            END DO
        END DO
    END SUBROUTINE comp_c


    ! Perform an update of fields in the RK time integration scheme
    ! dU_j = A_j*dU_(j-1) + dt*uo
    ! U_j = U_(j-1) + B_j*dU_j
    PURE SUBROUTINE rkstep(p, dp, rhsp, frhs, dtfu)
        ! Subroutine arguments
        REAL(realk), CONTIGUOUS, INTENT(inout) :: p(:)
        REAL(realk), CONTIGUOUS, INTENT(inout) :: dp(:)
        REAL(realk), CONTIGUOUS, INTENT(in) :: rhsp(:)
        REAL(realk), INTENT(in) :: frhs
        REAL(realk), INTENT(in) :: dtfu

        ! Local variables
        INTEGER(intk) :: i

        ! Perform the update in a manually crafted loop is faster than using an
        ! implicit loop, because of cache effects (dp(i) is already in cache
        ! when p(i) is updated)
        DO i = 1, SIZE(p)
            dp(i) = frhs*dp(i) + rhsp(i)
            p(i) = p(i) + dtfu*dp(i)
        END DO
    END SUBROUTINE rkstep

    ! Compute only particle displacement, not the final/ intermediate position, because
    ! the displacement might have to be influenced by boundary interactions.
    ! The routine that handles boundary interactions cannot be called here (in the core)
    ! dX_{j} = dt * B_{j} * ( A_{j} * dXeff_{j-1} + U(X_{j-1}) )
    PURE SUBROUTINE prkstep(dx_pot, dy_pot, dz_pot, u, v, w, dt, A, B, dx, dy, dz)

        ! Subroutine arguments
        REAL(realk), INTENT(inout) :: dx_pot, dy_pot, dz_pot
        REAL(realk), INTENT(in) :: u, v, w, dt
        REAL(realk), INTENT(in) :: A, B
        REAL(realk), INTENT(out) :: dx, dy, dz

        dx_pot = (A * dx_pot + dt * u)
        dy_pot = (A * dy_pot + dt * v)
        dz_pot = (A * dz_pot + dt * w)

        dx = B * dx_pot
        dy = B * dy_pot
        dz = B * dz_pot

    END SUBROUTINE prkstep

    ! ----- Butcher Table RK from here.... ------

    SUBROUTINE init_bt_rk(this, ctyp)

        ! subroutine arguments
        CLASS(bt_rk_t), INTENT(out) :: this
        CHARACTER(len=*), INTENT(in) :: ctyp

        SELECT CASE(lower(TRIM(ctyp)))
        CASE("euler")
            CALL this%init_euler_bt()
        CASE ("williamson")
            CALL this%init_williamson_bt()
        CASE DEFAULT
            WRITE(*, *) "Invalid: ", lower(TRIM(ctyp))
            CALL errr(__FILE__, __LINE__)
        END SELECT

    END SUBROUTINE init_bt_rk

    SUBROUTINE init_euler_bt(rk)

        ! Subroutine arguments
        CLASS(bt_rk_t), INTENT(out) :: rk

        ! Set to zero
        rk%c = 0.0
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        ! CAUTION: here, a and b refer to coefficients in the butcher table
        ! instead of the low stroage coefficients
        rk%nrk = 1
        rk%cflmax = oneeps*SQRT(3.0)

        rk%b(1) = 1.0

    END SUBROUTINE init_euler_bt

    SUBROUTINE init_williamson_bt(rk)

        ! Subroutine arguments
        CLASS(bt_rk_t), INTENT(out) :: rk

        ! Set to zero
        rk%c = 0.0
        rk%a = 0.0
        rk%b = 0.0

        ! Used coefficients
        ! CAUTION: here, a and b refer to coefficients in the butcher table
        ! instead of the low stroage coefficients
        rk%nrk = 3
        rk%cflmax = oneeps*SQRT(3.0)

        rk%c(1:rk%nrk) = [0.0, 1.0/3.0, 3.0/4.0]

        rk%b(1:rk%nrk) = [1.0/6.0, 3.0/10.0, 8.0/15.0]

        rk%a(2,1) = 1.0/3.0
        rk%a(3,1) = -3.0/16.0
        rk%a(3,2) = 15.0/16.0

    END SUBROUTINE init_williamson_bt

    SUBROUTINE get_bt_coefficients(this, c, b, a, irk)

        ! subroutine arguments
        CLASS(bt_rk_t), INTENT(in) :: this
        REAL(realk), OPTIONAL, ALLOCATABLE, INTENT(out) :: c(:)
        REAL(realk), OPTIONAL, ALLOCATABLE, INTENT(out) :: b(:)
        REAL(realk), OPTIONAL, ALLOCATABLE, INTENT(out) :: a(:,:)
        INTEGER(intk), OPTIONAL, INTENT(in) :: irk

        ! local variables
        INTEGER(intk) :: i, j

        IF (PRESENT(irk)) THEN
            IF (PRESENT(c)) THEN
                c = this%c(irk)
            END IF

            IF (PRESENT(b)) THEN
                b = this%b(irk)
            END IF

            ALLOCATE(a(1, irk-1))
            IF (PRESENT(a)) THEN
                IF (SIZE(a, 2) == 0) THEN
                    CONTINUE
                ELSE
                    DO i = 1, irk -1
                        a(1, i) = this%a(1, i)
                    END DO
                END IF
            END IF
        ELSE
            ALLOCATE(c(this%nrk))
            c = this%c(1:this%nrk)

            ALLOCATE(b(this%nrk))
            b = this%b(1:this%nrk)

            ALLOCATE(a(this%nrk, this%nrk))
            DO i = 1, this%nrk
                DO j = 1, this%nrk
                    a(i, j) = this%a(i, j)
                END DO
            END DO
        END IF

    END SUBROUTINE get_bt_coefficients

END MODULE rungekutta_mod
