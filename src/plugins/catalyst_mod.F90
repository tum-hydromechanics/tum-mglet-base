MODULE catalyst_mod
    USE core_mod
    USE catalyst_conduit
    USE catalyst_api
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR
    USE, INTRINSIC :: ISO_FORTRAN_ENV, only: stderr => error_unit

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE(config_t) :: cata_conf
    CHARACTER(len=mglet_filename_max), ALLOCATABLE :: catascripts(:)
    CHARACTER(len=mglet_filename_max) :: catalyst_path
    CHARACTER(len=8) :: CATALYST_IMPL = "paraview"
    INTEGER(intk) :: nscripts
    LOGICAL :: is_repr
    LOGICAL :: has_catalyst = .FALSE.

    PUBLIC :: init_catalyst, sample_catalyst, finish_catalyst
CONTAINS
    SUBROUTINE catalyst_adaptor_initialize()
        ! Local variables
        TYPE(C_PTR) node, script_node
        INTEGER(kind(catalyst_status)) :: err
        INTEGER(intk) :: iscript
        CHARACTER(len=10) :: scriptname

        node = catalyst_conduit_node_create()

        ! Catalyst setup
        CALL catalyst_conduit_node_set_path_char8_str(node, &
            "catalyst_load/implementation", CATALYST_IMPL)
        CALL catalyst_conduit_node_set_path_char8_str(node, &
            "catalyst_load/search_paths/paraview", catalyst_path)

        ! Representative dataset
        IF (is_repr) THEN
            CALL catalyst_conduit_node_set_path_char8_str(node, &
                "catalyst/pipelines/0/type", "io")
            CALL catalyst_conduit_node_set_path_char8_str(node, &
                "catalyst/pipelines/0/filename", "dataout.vthb")
            CALL catalyst_conduit_node_set_path_char8_str(node, &
                "catalyst/pipelines/0/channel", "grid")
        END IF

        ! Catalyst Pipeline Script
        IF (nscripts /= 0) THEN
            script_node = catalyst_conduit_node_fetch(node, "catalyst/scripts")
            DO iscript = 1, nscripts
                CALL iscript_to_string(iscript, scriptname)
                CALL catalyst_conduit_node_set_path_char8_str(script_node, &
                    scriptname // "/filename", catascripts(iscript))
            END DO
        END IF

        CALL catalyst_conduit_node_print(node)

        err = c_catalyst_initialize(node)
        IF (err /= catalyst_status_ok) THEN
            WRITE (stderr, *) "ERROR: Failed to initialize Catalyst: ", err
        END IF

        CALL catalyst_conduit_node_destroy(node)
    END SUBROUTINE catalyst_adaptor_initialize

    SUBROUTINE catalyst_adaptor_execute(itstep, ittot, timeph, dt)
        USE grids_mod, ONLY: minlevel, maxlevel, nmygridslvl
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(C_PTR) exec_node, channel_node, data_node, grid_node, &
            coords_node, topo_node, fields_node
        INTEGER(kind(catalyst_status)) :: err
        INTEGER(intk) :: ilvl, igridlvl, igrid, kk, jj, ii, shiftedlvl
        INTEGER(kind=8) :: nval
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, dx, dy, dz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: u_ptr, v_ptr, w_ptr
        CHARACTER(len=13) :: grid_name
        TYPE(field_t), POINTER :: u_c_f, v_c_f, w_c_f

        CALL get_field(u_c_f, "U_C")
        CALL get_field(v_c_f, "V_C")
        CALL get_field(w_c_f, "W_C")

        ! Catalyst state
        exec_node = catalyst_conduit_node_create()
        CALL catalyst_conduit_node_set_path_int32(exec_node, &
            "catalyst/state/timestep", itstep)
        CALL catalyst_conduit_node_set_path_float32(exec_node, &
            "catalyst/state/time", timeph)

        ! Grid/Channel type
        channel_node = catalyst_conduit_node_fetch(exec_node, &
            "catalyst/channels/grid")
        CALL catalyst_conduit_node_set_path_char8_str(channel_node, &
            "type", "amrmesh")

        data_node = catalyst_conduit_node_fetch(channel_node, "data")

        DO ilvl = minlevel, maxlevel
            DO igridlvl = 1, nmygridslvl(ilvl)
                igrid = mygridslvl(igridlvl, ilvl)
                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

                ! Grids state
                CALL igrid_to_string(igrid, grid_name)
                grid_node = catalyst_conduit_node_fetch(data_node, grid_name)
                CALL catalyst_conduit_node_set_path_int32(grid_node, &
                    "state/domain_id", igrid)
                CALL catalyst_conduit_node_set_path_int32(grid_node, &
                    "state/cycle", itstep)
                CALL catalyst_conduit_node_set_path_float32(grid_node, &
                    "state/time", timeph)
                shiftedlvl = ilvl - minlevel
                CALL catalyst_conduit_node_set_path_int32(grid_node, &
                    "state/level", shiftedlvl)

                ! Grid coordinates
                coords_node = catalyst_conduit_node_fetch(grid_node, &
                    "coordsets/coords")
                CALL catalyst_conduit_node_set_path_char8_str(coords_node, &
                    "type", "uniform")
                CALL catalyst_conduit_node_set_path_int32(coords_node, &
                    "dims/i", ii + 1 - 4)
                CALL catalyst_conduit_node_set_path_int32(coords_node, &
                    "dims/j", jj + 1 - 4)
                CALL catalyst_conduit_node_set_path_int32(coords_node, &
                    "dims/k", kk + 1 - 4)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "origin/x", minx)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "origin/y", miny)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "origin/z", minz)
                    dx = round_to_n_decimals((maxx - minx) / (ii - 4), 8)
                    dy = round_to_n_decimals((maxy - miny) / (jj - 4), 8)
                    dz = round_to_n_decimals((maxz - minz) / (kk - 4), 8)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "spacing/dx", dx)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "spacing/dy", dy)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "spacing/dz", dz)
            
                ! Grid topology
                topo_node = catalyst_conduit_node_fetch(grid_node, &
                    "topologies/mesh")
                CALL catalyst_conduit_node_set_path_char8_str(topo_node, &
                    "type", "uniform")
                CALL catalyst_conduit_node_set_path_char8_str(topo_node, &
                    "coordset", "coords")

                ! Grid fields
                fields_node = catalyst_conduit_node_fetch(grid_node, "fields")
                nval = (ii - 4) * (jj - 4) * (kk - 4)
                CALL u_c_f%get_ptr(u_ptr, igrid)
                CALL v_c_f%get_ptr(v_ptr, igrid)
                CALL w_c_f%get_ptr(w_ptr, igrid)
                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "u/association", "element")
                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "u/topology", "mesh")
                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "u/volume_dependent", "false")
                CALL catalyst_conduit_node_set_path_external_float32_ptr( &
                    fields_node, "u/values", u_ptr, nval)

                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "v/association", "element")
                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "v/topology", "mesh")
                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "v/volume_dependent", "false")
                CALL catalyst_conduit_node_set_path_external_float32_ptr( &
                    fields_node, "v/values", v_ptr, nval)

                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "w/association", "element")
                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "w/topology", "mesh")
                CALL catalyst_conduit_node_set_path_char8_str(fields_node, &
                    "w/volume_dependent", "false")
                CALL catalyst_conduit_node_set_path_external_float32_ptr( &
                    fields_node, "w/values", w_ptr, nval)
            END DO
        END DO

        CALL start_timer(813)
        err = c_catalyst_execute(exec_node)
        CALL stop_timer(813)
        IF (err /= catalyst_status_ok) THEN
            WRITE (stderr, *) "ERROR: Failed to execute Catalyst: ", err
        END IF

        CALL catalyst_conduit_node_destroy(exec_node)
    END SUBROUTINE catalyst_adaptor_execute

    SUBROUTINE catalyst_adaptor_finalize()
        TYPE(C_PTR) node
        INTEGER(kind(catalyst_status)) :: err
    
        node = catalyst_conduit_node_create()
        err = c_catalyst_finalize(node)
        IF (err /= catalyst_status_ok) THEN
          WRITE (stderr, *) "ERROR: Failed to finalize Catalyst: ", err
        END IF
    
        CALL catalyst_conduit_node_destroy(node)
    END SUBROUTINE catalyst_adaptor_finalize

    SUBROUTINE init_catalyst(ittot, mtstep, itint, timeph, dt, tend)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        INTEGER(intk), INTENT(in) :: itint
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: tend

        ! Local variables
        CHARACTER(len=64) :: jsonptr
        INTEGER(intk) :: i

        IF (.NOT. fort7%exists("/catalyst")) THEN
            RETURN
        END IF

        CALL set_timer(810, "CATALYST")
        CALL set_timer(811, "CATALYST_C_FIELDS")
        CALL set_timer(812, "CATALYST_ADAPTOR")
        CALL set_timer(813, "CATALYST_EXECUTE_LIB")

        ! Required values
        CALL fort7%get(cata_conf, "/catalyst")
        ! Catalyst Path
        IF ( cata_conf%is_char("/path") ) THEN
            CALL cata_conf%get_value("/path", catalyst_path)
        ELSE
            WRITE(*,*) "Specifiy directory containing libcatalyst-paraview.so"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Representative dataset
        IF ( cata_conf%is_logical("/repr") ) THEN
            CALL cata_conf%get_value("/repr", is_repr, .FALSE.)
        ELSE
            WRITE(*,*) "Specify repr to create representative dataset"
            is_repr = .FALSE.
        END IF

        ! Catalyst Scripts
        IF (cata_conf%exists("/scripts")) THEN
            CALL cata_conf%get_size("/scripts", nscripts)
            ALLOCATE(catascripts(nscripts))
            DO i = 1, nscripts
                WRITE(jsonptr, '("/scripts/", I0)') i-1
                CALL cata_conf%get_value(jsonptr, catascripts(i))
            END DO
        ELSE
            nscripts = 0
        END IF

        CALL catalyst_adaptor_initialize()

        ! Auxiliary fields for C-ordered arrays
        CALL set_field("U_C")
        CALL set_field("V_C")
        CALL set_field("W_C")

        ! Declaring initialized
        has_catalyst = .TRUE.
    END SUBROUTINE init_catalyst

    SUBROUTINE sample_catalyst(itstep, ittot, timeph, dt)
        IMPLICIT NONE
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        CALL start_timer(810)

        ! Create fields copy for C
        CALL start_timer(811)
        CALL field_to_c()
        CALL stop_timer(811)

        ! Execute adaptor
        CALL start_timer(812)
        CALL catalyst_adaptor_execute(itstep, ittot, timeph, dt)
        CALL stop_timer(812)

        CALL stop_timer(810)
    END SUBROUTINE

    SUBROUTINE finish_catalyst()
        CALL catalyst_adaptor_finalize()
        IF (ALLOCATED(catascripts)) DEALLOCATE(catascripts)
    END SUBROUTINE finish_catalyst

    SUBROUTINE igrid_to_string(igrid, output_string)
        INTEGER, INTENT(in) :: igrid
        CHARACTER(len=13), INTENT(out) :: output_string
        CHARACTER(len=8) :: temp_string
        WRITE(temp_string, '(I8.8)') igrid
        output_string = 'grid_' // temp_string
    END SUBROUTINE igrid_to_string

    SUBROUTINE iscript_to_string(iscript, output_string)
        INTEGER, INTENT(in) :: iscript
        CHARACTER(len=10), INTENT(out) :: output_string
        CHARACTER(len=3) :: temp_string
        WRITE(temp_string, '(I3.3)') iscript
        output_string = 'script_' // temp_string
    END SUBROUTINE iscript_to_string

    FUNCTION round_to_n_decimals(x, decimals) RESULT(rounded)
        REAL(realk), INTENT(in) :: x
        INTEGER(intk), INTENT(in) :: decimals
        REAL(realk) :: rounded
        REAL(realk) :: factor
    
        factor = 10.0 ** decimals    
        rounded = nint(x * factor) / factor
    END FUNCTION round_to_n_decimals

    SUBROUTINE field_to_c()
        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: u_f, v_f, w_f, u_c_f, v_c_f, w_c_f
        REAL(realk), POINTER, CONTIGUOUS :: arr_c(:, :, :), arr(:, :, :)
        INTEGER, PARAMETER :: nbl = 2

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(u_c_f, "U_C")
        CALL get_field(v_c_f, "V_C")
        CALL get_field(w_c_f, "W_C")

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
