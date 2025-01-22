MODULE catalyst_mod
    USE core_mod
    USE catalyst_conduit
    USE catalyst_api
    use, intrinsic :: iso_c_binding, only: C_PTR
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit

    IMPLICIT NONE(type, external)

    PRIVATE
    ! Module variables

    LOGICAL :: has_catalyst = .FALSE.
    TYPE(config_t) :: cata_conf

    CHARACTER(len=127) :: path_char
    CHARACTER(len=127) :: script_char
    CHARACTER(len=127) :: implementation_char
    LOGICAL :: repr_logical

    PUBLIC :: init_catalyst, sample_catalyst, finish_catalyst

CONTAINS
    subroutine catalyst_adaptor_initialize()
        type(C_PTR) node
        integer(kind(catalyst_status)) :: err

        node = catalyst_conduit_node_create()

        ! Catalyst setup
        call catalyst_conduit_node_set_path_char8_str(node, "catalyst_load/implementation", "paraview")
        call catalyst_conduit_node_set_path_char8_str(node, "catalyst_load/search_paths/paraview", "/usr/local/lib/catalyst/")

        ! Representative dataset
        IF (repr_logical) THEN
            call catalyst_conduit_node_set_path_char8_str(node, "catalyst/pipelines/0/type", "io")
            call catalyst_conduit_node_set_path_char8_str(node, "catalyst/pipelines/0/filename", "dataout.vthb")
            call catalyst_conduit_node_set_path_char8_str(node, "catalyst/pipelines/0/channel", "grid")
        END IF

        ! Catalyst Pipeline Script
        IF (.NOT. TRIM(script_char) == "") THEN
            call catalyst_conduit_node_set_path_char8_str(node, "catalyst/scripts/script/filename", script_char)
        END IF

        err = c_catalyst_initialize(node)
        if (err /= catalyst_status_ok) then
        write (stderr, *) "ERROR: Failed to initialize Catalyst: ", err
        end if

        call catalyst_conduit_node_destroy(node)
    end subroutine catalyst_adaptor_initialize

    subroutine catalyst_adaptor_execute(itstep, ittot, timeph, dt)
        USE grids_mod, ONLY: minlevel, maxlevel, nmygridslvl
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        type(C_PTR) exec_node, channel_node, data_node, grid_node, coords_node, topo_node, fields_node
        integer(kind(catalyst_status)) :: err
        INTEGER(intk) :: ilvl, igridlvl, igrid, kk, jj, ii, shiftedlvl
        INTEGER(kind=8) :: nval
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, dx, dy, dz
        character(len=13) :: grid_name
        type(field_t), POINTER :: u_c_f, v_c_f, w_c_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: u_ptr, v_ptr, w_ptr
        CALL get_field(u_c_f, "U_C")
        CALL get_field(v_c_f, "V_C")
        CALL get_field(w_c_f, "W_C")

        ! Function body
        exec_node = catalyst_conduit_node_create()
        call catalyst_conduit_node_set_path_int32(exec_node, "catalyst/state/timestep", itstep)
        call catalyst_conduit_node_set_path_float32(exec_node, "catalyst/state/time", timeph)

        channel_node = catalyst_conduit_node_fetch(exec_node, "catalyst/channels/grid")
        call catalyst_conduit_node_set_path_char8_str(channel_node, "type", "amrmesh")

        data_node = catalyst_conduit_node_fetch(channel_node, "data")

        DO ilvl = minlevel, maxlevel
            DO igridlvl = 1, nmygridslvl(ilvl)
                igrid = mygridslvl(igridlvl, ilvl)
                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

                CALL convert_to_string(igrid, grid_name)
                grid_node = catalyst_conduit_node_fetch(data_node, grid_name)
                call catalyst_conduit_node_set_path_int32(grid_node, "state/domain_id", igrid)
                call catalyst_conduit_node_set_path_int32(grid_node, "state/cycle", itstep)
                call catalyst_conduit_node_set_path_float32(grid_node, "state/time", timeph)
                
                shiftedlvl = ilvl - minlevel
                call catalyst_conduit_node_set_path_int32(grid_node, "state/level", shiftedlvl)
                dx = round_to_n_decimals((maxx - minx) / (ii - 4), 8)
                dy = round_to_n_decimals((maxy - miny) / (jj - 4), 8)
                dz = round_to_n_decimals((maxz - minz) / (kk - 4), 8)
                coords_node = catalyst_conduit_node_fetch(grid_node, "coordsets/coords")
                call catalyst_conduit_node_set_path_char8_str(coords_node, "type", "uniform")
                call catalyst_conduit_node_set_path_int32(coords_node, "dims/i", ii + 1 - 4)
                call catalyst_conduit_node_set_path_int32(coords_node, "dims/j", jj + 1 - 4)
                call catalyst_conduit_node_set_path_int32(coords_node, "dims/k", kk + 1 - 4)
                call catalyst_conduit_node_set_path_float32(coords_node, "origin/x", minx)
                call catalyst_conduit_node_set_path_float32(coords_node, "origin/y", miny)
                call catalyst_conduit_node_set_path_float32(coords_node, "origin/z", minz)
                call catalyst_conduit_node_set_path_float32(coords_node, "spacing/dx", dx)
                call catalyst_conduit_node_set_path_float32(coords_node, "spacing/dy", dy)
                call catalyst_conduit_node_set_path_float32(coords_node, "spacing/dz", dz)
            
                topo_node = catalyst_conduit_node_fetch(grid_node, "topologies/mesh")
                call catalyst_conduit_node_set_path_char8_str(topo_node, "type", "uniform")
                call catalyst_conduit_node_set_path_char8_str(topo_node, "coordset", "coords")

                fields_node = catalyst_conduit_node_fetch(grid_node, "fields")
                nval = (ii - 4) * (jj - 4) * (kk - 4)
                CALL u_c_f%get_ptr(u_ptr, igrid)
                CALL v_c_f%get_ptr(v_ptr, igrid)
                CALL w_c_f%get_ptr(w_ptr, igrid)
                call catalyst_conduit_node_set_path_char8_str(fields_node, "u/association", "element")
                call catalyst_conduit_node_set_path_char8_str(fields_node, "u/topology", "mesh")
                call catalyst_conduit_node_set_path_char8_str(fields_node, "u/volume_dependent", "false")
                call catalyst_conduit_node_set_path_external_float32_ptr(fields_node, "u/values", u_ptr, nval)

                call catalyst_conduit_node_set_path_char8_str(fields_node, "v/association", "element")
                call catalyst_conduit_node_set_path_char8_str(fields_node, "v/topology", "mesh")
                call catalyst_conduit_node_set_path_char8_str(fields_node, "v/volume_dependent", "false")
                call catalyst_conduit_node_set_path_external_float32_ptr(fields_node, "v/values", v_ptr, nval)

                call catalyst_conduit_node_set_path_char8_str(fields_node, "w/association", "element")
                call catalyst_conduit_node_set_path_char8_str(fields_node, "w/topology", "mesh")
                call catalyst_conduit_node_set_path_char8_str(fields_node, "w/volume_dependent", "false")
                call catalyst_conduit_node_set_path_external_float32_ptr(fields_node, "w/values", w_ptr, nval)
            END DO
        END DO

        !call catalyst_conduit_node_print(exec_node)

        err = c_catalyst_execute(exec_node)
        if (err /= catalyst_status_ok) then
            write (stderr, *) "ERROR: Failed to execute Catalyst: ", err
        end if

        call catalyst_conduit_node_destroy(exec_node)
    end subroutine catalyst_adaptor_execute

    subroutine catalyst_adaptor_finalize()
        type(C_PTR) node
        integer(kind(catalyst_status)) :: err
    
        node = catalyst_conduit_node_create()
        err = c_catalyst_finalize(node)
        if (err /= catalyst_status_ok) then
          write (stderr, *) "ERROR: Failed to finalize Catalyst: ", err
        end if
    
        call catalyst_conduit_node_destroy(node)
      end subroutine

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

        CALL set_timer(810, "CATALYST")
        CALL set_timer(811, "CATALYST_C_FIELDS")

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

        CALL catalyst_adaptor_initialize()
        !CALL pass_to_init(TRIM(script_char), TRIM(implementation_char), TRIM(path_char), repr_logical)

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

        CALL field_to_c()
        CALL catalyst_adaptor_execute(itstep, ittot, timeph, dt)
    END SUBROUTINE

    SUBROUTINE finish_catalyst()
        CALL catalyst_adaptor_finalize()
    END SUBROUTINE finish_catalyst

    subroutine convert_to_string(input_num, output_string)
        integer, intent(in) :: input_num
        character(len=13), intent(out) :: output_string
        character(len=8) :: temp_string
        write(temp_string, '(I8.8)') input_num
        output_string = 'grid_' // temp_string
    end subroutine convert_to_string

    function round_to_n_decimals(x, decimals) result(rounded)
        real(realk), intent(in) :: x
        integer(intk), intent(in) :: decimals
        real(realk) :: rounded
        real(realk) :: factor
    
        factor = 10.0 ** decimals    
        rounded = nint(x * factor) / factor
    end function round_to_n_decimals

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
