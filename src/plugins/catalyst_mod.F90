MODULE catalyst_mod

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_CHAR, C_INT, C_NULL_PTR, C_NULL_CHAR, C_BOOL
    USE core_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    INTERFACE
        SUBROUTINE catalyst_trigger( &
            p_mgdims, p_list_grids_lvl, p_mgbasb, p_get_bbox, p_get_parent, &
            p_get_arrptr, p_get_xyzptr, p_get_dxyzptr, p_get_ddxyzptr, &
            myid, numprocs, istep, &
            nscal, lvlmin, lvlmax ) BIND(C, name="catalyst_trigger")

            USE ISO_C_BINDING, ONLY: c_funptr, c_int

            ! pointers must be passed with VALUE (it is already a pointer)
            TYPE(c_funptr), INTENT(in), VALUE :: &
                p_mgdims, p_list_grids_lvl, p_mgbasb, p_get_bbox, p_get_parent, &
                p_get_arrptr, p_get_xyzptr, p_get_dxyzptr, p_get_ddxyzptr

            ! rest is handled as pointers
            INTEGER(kind=c_int), INTENT(in) :: &
                myid, numprocs, istep, &
                nscal, lvlmin, lvlmax

        END SUBROUTINE catalyst_trigger
    END INTERFACE


    INTERFACE
        SUBROUTINE catalyst_init(file, impl, path, is_repr, myid) BIND(C, name="catalyst_init")
            USE ISO_C_BINDING, ONLY: C_CHAR, C_BOOL, C_INT

            CHARACTER(C_CHAR), INTENT(IN) :: file(*)
            CHARACTER(C_CHAR), INTENT(IN) :: impl(*)
            CHARACTER(C_CHAR), INTENT(IN) :: path(*)
            LOGICAL(C_BOOL), INTENT(IN) :: is_repr
            integer(C_INT), INTENT(IN) :: myid

        END SUBROUTINE catalyst_init
    END INTERFACE

    INTERFACE
        SUBROUTINE catalyst_finish() BIND(C, name="catalyst_finish")
        END SUBROUTINE
    END INTERFACE

    ! Module variables

    LOGICAL :: has_catalyst = .FALSE.
    TYPE(config_t) :: cata_conf

    CHARACTER(len=127) :: path_char
    CHARACTER(len=127) :: script_char
    CHARACTER(len=127) :: implementation_char
    LOGICAL :: repr_logical

    PUBLIC :: init_catalyst, sample_catalyst, finish_catalyst

CONTAINS

    SUBROUTINE init_catalyst(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        IF (.NOT. fort7%exists("/catalyst")) THEN
            RETURN
        END IF

        ! Required values
        CALL fort7%get(cata_conf, "/catalyst")
        ! Catalyst Path
        IF ( cata_conf%is_char("/path") ) THEN
            CALL cata_conf%get_value("/path", path_char)
        ELSE
            WRITE(*,*) "Specifiy path of directory that contains libcatalyst-paraview.so"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Representative dataset
        IF ( cata_conf%is_logical("/repr") ) THEN
            CALL cata_conf%get_value("/repr", repr_logical, .FALSE.)
        ELSE
            WRITE(*,*) "Specify repr to create representative dataset for Catalyst."
            repr_logical = .FALSE.
        END IF

        ! Catalyst Script
        IF ( cata_conf%is_char("/script") ) THEN
            CALL cata_conf%get_value("/script", script_char)
        ELSE
            IF (.NOT. repr_logical) THEN
                WRITE(*,*) "Specifiy script for Catalyst is required when not creating representative dataset."
                CALL errr(__FILE__, __LINE__)
            ELSE
                script_char = ""
            END IF
        END IF
        
        ! Only implemented for Paraview
        implementation_char = "paraview"

        CALL pass_to_init(TRIM(script_char), TRIM(implementation_char), TRIM(path_char), repr_logical)

        ! Auxiliary fields for C-ordered arrays
        CALL set_field("U_C")
        CALL set_field("V_C")
        CALL set_field("W_C")

        ! Declaring initialized
        has_catalyst = .TRUE.

    END SUBROUTINE init_catalyst


    SUBROUTINE pass_to_init(file, impl, path, repr)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: file
        CHARACTER(len=*), INTENT(in) :: impl
        CHARACTER(len=*), INTENT(in) :: path
        LOGICAL, INTENT(in) :: repr

        ! Local variables
        CHARACTER(c_char), DIMENSION(LEN(file)+1) :: c_file
        CHARACTER(c_char), DIMENSION(LEN(impl)+1) :: c_impl
        CHARACTER(c_char), DIMENSION(LEN(path)+1) :: c_path
        LOGICAL(c_bool) :: c_repr
        INTEGER(c_int) :: c_myid

        ! Add trailing C_NULL_CHAR to file
        c_file(1:LEN(file)) = TRANSFER(file, c_file)
        c_file(LEN_TRIM(file)+1) = C_NULL_CHAR

        ! Add trailing C_NULL_CHAR to file
        c_impl(1:LEN(impl)) = TRANSFER(impl, c_impl)
        c_impl(LEN_TRIM(impl)+1) = C_NULL_CHAR        

        ! Add trailing C_NULL_CHAR to path
        c_path(1:LEN(path)) = TRANSFER(path, c_path)
        c_path(LEN_TRIM(path)+1) = C_NULL_CHAR        

        ! Conversion of repr to LOGICAL(C_BOOL)
        c_repr = MERGE(.TRUE., .FALSE., repr)

        ! Conversion of myid to C_INT
        c_myid = INT( myid, kind=C_INT )

        CALL catalyst_init( c_file, c_impl, c_path, c_repr, c_myid )

    END SUBROUTINE pass_to_init


    SUBROUTINE finish_catalyst()
        CALL catalyst_finish()
    END SUBROUTINE finish_catalyst


    SUBROUTINE sample_catalyst(itstep, ittot, timeph, dt)
        USE ISO_C_BINDING, ONLY: C_FUNPTR, C_FUNLOC, C_INT
        IMPLICIT NONE
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        TYPE(field_t), POINTER :: u_c_f, v_c_f, w_c_f

        ! local variables
        TYPE(C_FUNPTR) :: cp_mgdims, cp_iterate_grids_lvl, &
                          cp_mgbasb, cp_get_bbox, cp_get_parent, &
                          cp_get_arrptr, cp_get_xyzptr, &
                          cp_get_dxyzptr, cp_get_ddxyzptr

        INTEGER(kind=C_INT) :: c_myid, c_numprocs, c_istep, &
                               c_nscal, c_lvlmin, c_lvlmax

        IF (.NOT. has_catalyst) RETURN

        ! setting the function pointers
        ! (necessary to transfer MGLET functions to the library)
        cp_mgdims = C_FUNLOC( c_mgdims )
        cp_iterate_grids_lvl = C_FUNLOC( c_iterate_grids_lvl )
        cp_mgbasb = C_FUNLOC( c_mgbasb )
        cp_get_bbox = C_FUNLOC( c_get_bbox )
        cp_get_parent = C_FUNLOC(c_get_parent)

        cp_get_arrptr = C_FUNLOC( c_get_arrptr )
        cp_get_xyzptr = C_FUNLOC( c_get_xyzptr )
        cp_get_dxyzptr = C_FUNLOC( c_get_dxyzptr )
        cp_get_ddxyzptr = C_FUNLOC( c_get_ddxyzptr )

        ! convert the remaining integer to ensure consistency
        c_myid = INT( myid, kind=C_INT )
        c_numprocs = INT( numprocs, kind=C_INT )
        c_istep = INT( itstep, kind=C_INT )
        c_nscal = INT( 1, kind=C_INT )   ! TO DO (!!!)
        c_lvlmin = INT( minlevel, kind=C_INT )
        c_lvlmax = INT( maxlevel, kind=C_INT )
        

        ! preparing the C-ordered fields
        CALL get_field(u_c_f, "U_C")
        CALL get_field(v_c_f, "V_C")
        CALL get_field(w_c_f, "W_C")
        CALL field_to_c(u_c_f, v_c_f, w_c_f)

        ! calling the C function described in the interface
        CALL catalyst_trigger( &
            cp_mgdims, cp_iterate_grids_lvl, cp_mgbasb, cp_get_bbox, cp_get_parent, &
            cp_get_arrptr, cp_get_xyzptr, cp_get_dxyzptr, cp_get_ddxyzptr, &
            c_myid, c_numprocs, c_istep, &
            c_nscal, c_lvlmin, c_lvlmax )

    END SUBROUTINE sample_catalyst


    ! Return grid dimensions with INT(kind=4) arguments
    SUBROUTINE c_mgdims(kk, jj, ii, igrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        USE core_mod, ONLY: intk, get_mgdims
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(out) :: kk, jj, ii
        INTEGER(kind=c_int), INTENT(in) :: igrid
        INTEGER(kind=intk) :: kk_tmp, jj_tmp, ii_tmp, igrid_tmp
        ! Function body
        igrid_tmp = INT(igrid, kind=intk)
        CALL get_mgdims( kk_tmp, jj_tmp, ii_tmp, igrid_tmp )
        kk = INT(kk_tmp, kind=c_int)
        jj = INT(jj_tmp, kind=c_int)
        ii = INT(ii_tmp, kind=c_int)
    END SUBROUTINE c_mgdims


    ! Returns the boundary type at front, back, right, left, bottom, top for igrid
    SUBROUTINE c_mgbasb(cnfro, cnbac, cnrgt, cnlft, cnbot, cntop, cigrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        USE core_mod, ONLY: intk, get_mgbasb
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: cigrid
        INTEGER(kind=c_int), INTENT(out) :: cnfro, cnbac, cnrgt, cnlft, cnbot, cntop
        INTEGER(kind=intk) :: nfro, nbac, nrgt, nlft, nbot, ntop, igrid
        ! Function body
        igrid = INT(cigrid, kind=intk)
        ! calling the MGLET routine
        CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
        ! transforming back for C interoperability
        cnfro = INT(nfro, kind=c_int)
        cnbac = INT(nbac, kind=c_int)
        cnrgt = INT(nrgt, kind=c_int)
        cnlft = INT(nlft, kind=c_int)
        cnbot = INT(nbot, kind=c_int)
        cntop = INT(ntop, kind=c_int)
    END SUBROUTINE c_mgbasb


    ! Facilitates iteration over local grids on certain level with INT(kind=4) arguments
    SUBROUTINE c_iterate_grids_lvl(igrid, num, ilevel) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        USE core_mod, ONLY: intk, nmygridslvl, mygridslvl
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: ilevel
        INTEGER(kind=c_int), INTENT(in) :: num
        INTEGER(kind=c_int), INTENT(out) :: igrid
        INTEGER(kind=intk) :: igrid_tmp, max
        ! Function body
        max = nmygridslvl( INT(ilevel,kind=intk) )
        IF ( INT(num,kind=intk) > max ) THEN
            igrid = INT( -1, kind=c_int ) ! indicates end
        ELSE
            igrid_tmp = mygridslvl(num, ilevel)
            igrid = INT( igrid_tmp, kind=c_int )
        END IF
    END SUBROUTINE c_iterate_grids_lvl


    ! Facilitates iteration over local grids on certain level with INT(kind=4) arguments
    SUBROUTINE c_get_arrptr(c_arr_ptr, name, igrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int, c_ptr, c_char, c_loc
        USE core_mod, ONLY: intk, realk, get_field, field_t
        IMPLICIT NONE
        ! Variables
        TYPE(c_ptr), INTENT(in) :: name
        INTEGER(kind=c_int), INTENT(in) :: igrid
        TYPE(c_ptr), INTENT(out) :: c_arr_ptr
        ! Local variables
        INTEGER(kind=intk) :: igrid_tmp
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_arr_ptr(:,:,:)
        TYPE(field_t), POINTER :: field
        CHARACTER(len=36) :: f_name
        ! Function body
        igrid_tmp = INT(igrid, kind=intk)
        CALL c_string_f_char(name, f_name)
        CALL get_field(field, TRIM(f_name))
        CALL field%get_ptr(f_arr_ptr, igrid_tmp)
        c_arr_ptr = C_LOC( f_arr_ptr )
    END SUBROUTINE c_get_arrptr


    ! Function returns pointers for the 1D fields X, Y, Z
    SUBROUTINE c_get_xyzptr(c_x_ptr, c_y_ptr, c_z_ptr, igrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int, c_ptr, c_loc
        USE core_mod, ONLY: intk, realk, get_field, field_t
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: igrid
        TYPE(c_ptr), INTENT(out) :: c_x_ptr
        TYPE(c_ptr), INTENT(out) :: c_y_ptr
        TYPE(c_ptr), INTENT(out) :: c_z_ptr
        ! Local variables
        INTEGER(kind=intk) :: igrid_tmp
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_x_ptr(:)
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_y_ptr(:)
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_z_ptr(:)
        TYPE(field_t), POINTER :: x_f, y_f, z_f
        ! Function body
        igrid_tmp = INT(igrid, kind=intk)
        CALL get_field(x_f, 'X')
        CALL x_f%get_ptr(f_x_ptr, igrid_tmp)
        c_x_ptr = C_LOC(f_x_ptr)
        CALL get_field(y_f, 'Y')
        CALL y_f%get_ptr(f_y_ptr, igrid_tmp)
        c_y_ptr = C_LOC(f_y_ptr)
        CALL get_field(z_f, 'Z')
        CALL z_f%get_ptr(f_z_ptr, igrid_tmp)
        c_z_ptr = C_LOC(f_z_ptr)
    END SUBROUTINE c_get_xyzptr


    ! Function returns pointers for the 1D fields DX, DY, DZ
    SUBROUTINE c_get_dxyzptr(c_dx_ptr, c_dy_ptr, c_dz_ptr, igrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int, c_ptr, c_loc
        USE core_mod, ONLY: intk, realk, get_field, field_t
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: igrid
        TYPE(c_ptr), INTENT(out) :: c_dx_ptr
        TYPE(c_ptr), INTENT(out) :: c_dy_ptr
        TYPE(c_ptr), INTENT(out) :: c_dz_ptr
        ! Local variables
        INTEGER(kind=intk) :: igrid_tmp
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_dx_ptr(:)
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_dy_ptr(:)
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_dz_ptr(:)
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
        ! Function body
        igrid_tmp = INT(igrid, kind=intk)
        CALL get_field(dx_f, 'DX')
        CALL dx_f%get_ptr(f_dx_ptr, igrid_tmp)
        c_dx_ptr = C_LOC(f_dx_ptr)
        CALL get_field(dy_f, 'DY')
        CALL dy_f%get_ptr(f_dy_ptr, igrid_tmp)
        c_dy_ptr = C_LOC(f_dy_ptr)
        CALL get_field(dz_f, 'DZ')
        CALL dz_f%get_ptr(f_dz_ptr, igrid_tmp)
        c_dz_ptr = C_LOC(f_dz_ptr)
    END SUBROUTINE c_get_dxyzptr


    ! Function returns pointers for the 1D fields DDX, DDY, DDZ
    SUBROUTINE c_get_ddxyzptr(c_ddx_ptr, c_ddy_ptr, c_ddz_ptr, c_igrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int, c_ptr, c_loc
        USE core_mod, ONLY: intk, realk, get_field, field_t
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: c_igrid
        TYPE(c_ptr), INTENT(out) :: c_ddx_ptr
        TYPE(c_ptr), INTENT(out) :: c_ddy_ptr
        TYPE(c_ptr), INTENT(out) :: c_ddz_ptr
        ! Local variables
        INTEGER(kind=intk) :: igrid_tmp
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_ddx_ptr(:)
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_ddy_ptr(:)
        REAL(kind=realk), POINTER, CONTIGUOUS :: f_ddz_ptr(:)
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        ! Function body
        igrid_tmp = INT(c_igrid, kind=intk)
        CALL get_field(ddx_f, 'DDX')
        CALL ddx_f%get_ptr(f_ddx_ptr, igrid_tmp)
        c_ddx_ptr = C_LOC(f_ddx_ptr)
        CALL get_field(ddy_f, 'DDY')
        CALL ddy_f%get_ptr(f_ddy_ptr, igrid_tmp)
        c_ddy_ptr = C_LOC(f_ddy_ptr)
        CALL get_field(ddz_f, 'DDZ')
        CALL ddz_f%get_ptr(f_ddz_ptr, igrid_tmp)
        c_ddz_ptr = C_LOC(f_ddz_ptr)
    END SUBROUTINE c_get_ddxyzptr

    ! Return the parent grid id for igrid
    SUBROUTINE c_get_parent(c_ipargrid, c_igrid) bind(C)
        USE grids_mod, ONLY: iparent
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(IN) :: c_igrid
        INTEGER(kind=c_int), INTENT(out) :: c_ipargrid
        ! Local variables
        INTEGER(kind=intk) :: igrid
        INTEGER(kind=intk) :: ipar
        ! Function body
        igrid = INT(c_igrid, kind=intk)
        ipar = iparent(igrid)
        c_ipargrid = INT(ipar, kind=c_int)
    END SUBROUTINE

    ! Returns the boundary type at front, back, right, left, bottom, top for igrid
    SUBROUTINE c_get_bbox(c_minx, c_maxx, c_miny, c_maxy, c_minz, c_maxz, c_igrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int, c_float
        USE core_mod, ONLY: intk, get_mgbasb, get_bbox
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: c_igrid
        REAL(kind=c_float), INTENT(out) :: c_minx, c_maxx, c_miny, c_maxy, c_minz, c_maxz
        REAL(kind=realk) :: minx, maxx, miny, maxy, minz, maxz
        INTEGER :: igrid
        ! Function body
        igrid = INT(c_igrid, kind=intk)
        ! calling the MGLET routine
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
        ! transforming back for C interoperability
        c_minx = REAL(minx, kind=c_float)
        c_maxx = REAL(maxx, kind=c_float)
        c_miny = REAL(miny, kind=c_float)
        c_maxy = REAL(maxy, kind=c_float)
        c_minz = REAL(minz, kind=c_float)
        c_maxz = REAL(maxz, kind=c_float)
    END SUBROUTINE c_get_bbox


    ! found at https://stackoverflow.com/questions/41247242/c-and-fortran-interoperability-for-strings
    SUBROUTINE c_string_f_char(C_string, F_string)
        USE ISO_C_BINDING
        IMPLICIT NONE
        ! Variables
        type(C_PTR), intent(in) :: C_string
        character(len=*), intent(out) :: F_string
        character(len=1, kind=C_CHAR), dimension(:), pointer :: p_chars
        integer :: i
        if (.not. C_associated(C_string)) then
            F_string = ' '
        else
            call C_F_pointer(C_string, p_chars, [huge(0)])
            do i = 1, len(F_string)
                if (p_chars(i) == C_NULL_CHAR) exit
                    F_string(i:i) = p_chars(i)
                end do
            if (i <= len(F_string)) F_string(i:) = ' '
        end if
    END SUBROUTINE


    SUBROUTINE field_to_c(u_c_f, v_c_f, w_c_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_c_f, v_c_f, w_c_f
        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: u_f, v_f, w_f
        REAL(realk), POINTER, CONTIGUOUS :: arr_c(:, :, :), arr(:, :, :)
        INTEGER, PARAMETER :: nbl = 2
        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            ! treating U
            CALL u_f%get_ptr(arr, igrid)
            CALL u_c_f%get_ptr(arr_c, igrid)
            CALL field_to_c_grid(kk, jj, ii, nbl, arr_c, arr)
            ! treating V
            CALL v_f%get_ptr(arr, igrid)
            CALL v_c_f%get_ptr(arr_c, igrid)
            CALL field_to_c_grid(kk, jj, ii, nbl, arr_c, arr)
            ! treating W
            CALL w_f%get_ptr(arr, igrid)
            CALL w_c_f%get_ptr(arr_c, igrid)
            CALL field_to_c_grid(kk, jj, ii, nbl, arr_c, arr)
        END DO

    END SUBROUTINE field_to_c


    PURE SUBROUTINE field_to_c_grid(kk, jj, ii, nbl, arr_c, arr_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, nbl
        REAL(realk), INTENT(inout) :: arr_c(ii-2*nbl, jj-2*nbl, kk-2*nbl)
        REAL(realk), INTENT(in) :: arr_f(kk, jj, ii)
        ! Local variables
        INTEGER :: k, j, i
        ! Inversion of indices for usage in C
        DO i = 1, ii-2*nbl
            DO j = 1, jj-2*nbl
                DO k = 1, kk-2*nbl
                    arr_c(i, j, k) = arr_f(k+nbl, j+nbl, i+nbl)
                END DO
            END DO
        END DO
    END SUBROUTINE field_to_c_grid


END MODULE catalyst_mod





    ! ! Facilitates iteration over local grids on certain level with INT(kind=4) arguments
    ! SUBROUTINE c_get_arrptr(c_arr_ptr, name, igrid) bind(C)
    !     USE ISO_C_BINDING, ONLY: c_int, c_ptr, c_char, c_loc
    !     USE core_mod, ONLY: intk, nmygridslvl, mygridslvl, get_fieldptr, c_realk
    !     IMPLICIT NONE
    !     ! Variables
    !     CHARACTER(len=1,kind=c_char), INTENT(in) :: name
    !     INTEGER(kind=c_int), INTENT(in) :: igrid
    !     TYPE(c_ptr), VALUE :: c_arr_ptr
    !     ! Local variables
    !     INTEGER(kind=intk) :: igrid_tmp
    !     REAL(realk), POINTER, CONTIGUOUS :: f_arr_ptr(:,:,:)
    !     ! Function body
    !     WRITE(*,*) TRIM(name)
    !     igrid_tmp = INT(igrid, kind=intk)
    !     CALL get_fieldptr(f_arr_ptr, TRIM(name), igrid_tmp)
    !     WRITE(*,*) f_arr_ptr(3,3,3);
    !     c_arr_ptr = C_LOC( f_arr_ptr )
    ! END SUBROUTINE c_get_arrptr




    ! ! Computation of the lambda2 field
    ! SUBROUTINE c_compute_lambda2(ilevel) bind(C)
    !     USE ISO_C_BINDING, ONLY: c_int
    !     USE setmpi_mod, ONLY: nmygridslvl, mygridslvl
    !     USE pointer_mod, ONLY: get_ip1, get_ip3
    !     USE fields_mod, ONLY: u, v, w, dx, dy, dz, ddx, ddy, ddz, hilf
    !     USE connect2_mod, ONLY: connect2 => connect
    !     IMPLICIT NONE
    !     ! Variables
    !     INTEGER(kind=c_int), INTENT(in) :: ilevel
    !     INTEGER :: igrid, i, ii, jj, kk, ip1, ip3
    !     ! Function body
    !     CALL connect2( ilevel, 2, v1=u, v2=v, v3=w, corners=.TRUE. );
    !     hilf = 0.0;
    !     DO i = 1, nmygridslvl(ilevel)
    !         igrid = mygridslvl(i, ilevel)
    !         CALL mgdims(kk, jj, ii, igrid)
    !         CALL get_ip1(ip1, igrid)
    !         CALL get_ip3(ip3, igrid)
    !         CALL lambda2_grid(kk, jj, ii, u(ip3), v(ip3), w(ip3), &
    !         dx(ip1), dy(ip1), dz(ip1), ddx(ip1), ddy(ip1), ddz(ip1), hilf(ip3))
    !         CALL maskbp(kk, jj, ii, hilf(ip3), bp(ip3))
    !     END DO
    !     CALL connect2( ilevel, 2, s1=hilf, corners=.TRUE. )
    ! END SUBROUTINE c_compute_lambda2

    ! SUBROUTINE lambda2_grid( kk, jj, ii, u, v, w, dx, dy, dz, ddx, ddy, ddz, lambda2 )
    !     USE precision_mod, only: realk
    !     IMPLICIT NONE
    !     INTEGER, INTENT(in) :: kk, jj, ii
    !     REAL(kind=realk), INTENT(in) :: u(kk,jj,ii), v(kk,jj,ii), w(kk,jj,ii)
    !     REAL(kind=realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
    !     REAL(kind=realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
    !     REAL(kind=realk), INTENT(out) :: lambda2(kk,jj,ii)
    !     ! local variables
    !     INTEGER :: i, j, k, it, jt
    !     INTEGER, PARAMETER :: d = 1
    !     REAL(kind=realk) :: up, uc, um, vp, vc, vm, wp, wc, wm ! p = plus, c = center, m = minus (in direction of derivative)
    !     REAL(kind=realk) :: e1, e2, e3, eig1, eig2, eig3, p1, p2, p3, q, p, r, phi
    !     REAL(kind=realk), PARAMETER :: PI = ACOS(-1.0)
    !     ! tensors in C indexing (this allows to translate routines from MGTOOLS)
    !     REAL(kind=realk) :: gradient(0:8)
    !     REAL(kind=realk) :: tensorS(0:2,0:2), tensorW(0:2,0:2), A(0:2,0:2), B(0:2,0:2)

    !     lambda2 = 0.0;
    !     DO i = 1+d, ii-d, 1
    !         DO j = 1+d, jj-d, 1
    !             DO k = 1+d, kk-d, 1
    !                 ! compute the velocity gradient tensor
    !                 ! [dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz]
    !                 gradient(0) = ( u(k,j,i) - u(k,j,i-1) ) / ddx(i)
    !                 um = 0.5 * (u(k,j-1,i) + u(k,j-1,i-1))
    !                 uc = 0.5 * (u(k,j,  i) + u(k,j,  i-1))
    !                 up = 0.5 * (u(k,j+1,i) + u(k,j+1,i-1))
    !                 gradient(1) = 0.5 * (uc-um)/dy(j-1) + 0.5 * (up-uc)/dy(j)
    !                 um = 0.5 * (u(k-1,j,i) + u(k-1,j,i-1))
    !                 uc = 0.5 * (u(k  ,j,i) + u(k,  j,i-1))
    !                 up = 0.5 * (u(k+1,j,i) + u(k+1,j,i-1))
    !                 gradient(2) = 0.5 * (uc-um)/dz(k-1) + 0.5 * (up-uc)/dz(k)

    !                 vm = 0.5 * (v(k,j,i-1) + v(k,j-1,i-1))
    !                 vc = 0.5 * (v(k,j,i  ) + v(k,j-1,i  ))
    !                 vp = 0.5 * (v(k,j,i+1) + v(k,j-1,i+1))
    !                 gradient(3) = 0.5 * (vc-vm)/dx(i-1) + 0.5 * (vp-vc)/dx(i)
    !                 gradient(4) = ( v(k,j,i) - v(k,j-1,i) ) / ddy(j)
    !                 vm = 0.5 * (v(k-1,j,i) + v(k-1,j-1,i))
    !                 vc = 0.5 * (v(k  ,j,i) + v(k  ,j-1,i))
    !                 vp = 0.5 * (v(k+1,j,i) + v(k+1,j-1,i))
    !                 gradient(5) = 0.5 * (vc-vm)/dz(k-1) + 0.5 * (vp-vc)/dz(k)

    !                 wm = 0.5 * (w(k,j,i-1) + w(k-1,j,i-1))
    !                 wc = 0.5 * (w(k,j,i  ) + w(k-1,j,i  ))
    !                 wp = 0.5 * (w(k,j,i+1) + w(k-1,j,i+1))
    !                 gradient(6) = 0.5 * (wc-wm)/dx(i-1) + 0.5 * (wp-wc)/dx(i)
    !                 wm = 0.5 * (w(k,j-1,i) + w(k-1,j-1,i))
    !                 wc = 0.5 * (w(k,j,  i) + w(k-1,j,  i))
    !                 wp = 0.5 * (w(k,j+1,i) + w(k-1,j+1,i))
    !                 gradient(7) = 0.5 * (wc-wm)/dy(j-1) + 0.5 * (wp-wc)/dy(j)
    !                 gradient(8) = ( w(k,j,i) - w(k-1,j,i) ) / ddz(k)

    !                 ! Computing strain and rotation rate tensor
    !                 DO jt = 0, 2
    !                     DO it = 0, 2
    !                         tensorS(it,jt) = 0.5 * ( gradient(it*3+jt) + gradient(jt*3+it) )
    !                         tensorW(it,jt) = 0.5 * ( gradient(it*3+jt) - gradient(jt*3+it) )
    !                     END DO
    !                 END DO
    !                 ! Computing matrix A = S2 + Om2
    !                 DO jt = 0, 2
    !                     DO it = 0, 2
    !                         A(it,jt) = &
    !                         tensorS(0,jt)*tensorS(it,0) + tensorS(1,jt)*tensorS(it,1) + tensorS(2,jt)*tensorS(it,2) + &
    !                         tensorW(0,jt)*tensorW(it,0) + tensorW(1,jt)*tensorW(it,1) + tensorW(2,jt)*tensorW(it,2)
    !                     END DO
    !                 END DO
    !                 p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
    !                 IF ( p1 == 0.0 ) THEN
    !                     e1 = A(0,0); e2 = A(1,1); e3 = A(2,2)
    !                     eig2 = e1 + e2 + e3 - MAX(MAX(e1, e2), e3) - MIN(MIN(e1, e2), e3)
    !                 ELSE
    !                     q = ( A(0,0) + A(1,1) + A(2,2) ) / 3.0
    !                     p2 = (A(0,0) - q) * (A(0,0) - q) + &
    !                          (A(1,1) - q) * (A(1,1) - q) + &
    !                          (A(2,2) - q) * (A(2,2) - q) + 2.0 * p1;
    !                     p = SQRT( p2 / 6.0 );
    !                     ! Filling matrix B
    !                     DO jt = 0, 2
    !                         DO it = 0, 2
    !                             IF ( it == jt ) THEN
    !                                 B(it,jt) = ( 1.0 / p ) * ( A(it,jt) - q )
    !                             ELSE
    !                                 B(it,jt) = ( 1.0 / p ) * A(it,jt)
    !                             END IF
    !                         END DO
    !                     END DO
    !                     r =  ( B(0,0) * ( B(1,1)*B(2,2) - B(1,2)*B(2,1) ) &
    !                          - B(0,1) * ( B(1,0)*B(2,2) - B(1,2)*B(2,0) ) &
    !                          + B(0,2) * ( B(1,0)*B(2,1) - B(2,2)*B(2,0) ) ) / 2.0
    !                     ! Correction of possible small precision errors
    !                     IF ( r <= -1.0 ) THEN
    !                         phi = pi / 3.0
    !                     ELSE IF ( r >= 1.0 ) THEN
    !                         phi = 0.0
    !                     ELSE
    !                         phi = ACOS(r) / 3.0
    !                     END IF
    !                     ! Eigenvalues eig3 <= eig2 <= eig1
    !                     eig1 = q + 2.0 * p * COS(phi)
    !                     eig3 = q + 2.0 * p * COS(phi + pi*(2.0/3.0))
    !                     eig2 = 3.0 * q - eig1 - eig3
    !                 END IF
    !                 lambda2(k,j,i) = eig2;
    !             END DO
    !         END DO
    !     END DO
    ! END SUBROUTINE lambda2_grid


    ! ! Computation of the maximum value of finecell inside the grid
    ! SUBROUTINE c_compute_max_finecell(igrid, imax) bind(C)
    !     USE ISO_C_BINDING, ONLY: c_int
    !     USE pointer_mod, ONLY: get_ip3
    !     USE fields_mod, ONLY: finecell
    !     IMPLICIT NONE
    !     ! Variables
    !     INTEGER(kind=c_int), INTENT(in) :: igrid
    !     INTEGER(kind=c_int), INTENT(out) :: imax
    !     INTEGER :: ii, jj, kk, ip3, result
    !     ! Function body
    !     CALL mgdims(kk, jj, ii, INT(igrid))
    !     CALL get_ip3(ip3, INT(igrid))
    !     CALL max_finecell_grid(finecell(ip3), kk, jj, ii, result)
    !     imax = INT( result, kind=c_int )
    ! END SUBROUTINE c_compute_max_finecell

    ! SUBROUTINE max_finecell_grid(finecell, kk, jj, ii, result)
    !     IMPLICIT NONE
    !     ! Variables
    !     INTEGER, INTENT(in) :: kk, jj, ii
    !     REAL(kind=realk), INTENT(in) :: finecell(kk,jj,ii)
    !     INTEGER, INTENT(out) :: result
    !     REAL(kind=realk) :: val
    !     ! Function body
    !     val = MAXVAL( finecell( 3:(kk-2), 3:(jj-2), 3:(ii-2) ) )
    !     result = 0
    !     IF ( val > 0.5 ) THEN
    !         result = 1
    !     END IF
    ! END SUBROUTINE max_finecell_grid
