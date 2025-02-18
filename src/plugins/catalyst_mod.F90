MODULE catalyst_mod
    USE core_mod
    USE catalyst_conduit
    USE catalyst_api
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR
    USE, INTRINSIC :: ISO_FORTRAN_ENV, only: stderr => error_unit

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE :: dummy_field_t
        REAL(realk), ALLOCATABLE :: values(:)
        INTEGER(kind=8) :: nval
    END TYPE dummy_field_t

    TYPE(config_t) :: cata_conf
    TYPE(dummy_field_t), ALLOCATABLE, DIMENSION(:) :: dummy_fields
    CHARACTER(len=mglet_filename_max), ALLOCATABLE :: catascripts(:)
    CHARACTER(len=nchar_name), ALLOCATABLE :: fields(:)
    CHARACTER(len=mglet_filename_max) :: catalyst_path
    CHARACTER(len=8) :: CATALYST_IMPL = "paraview"
    INTEGER(intk) :: nscripts, nfields
    LOGICAL :: is_repr
    LOGICAL :: has_catalyst = .FALSE.

    PUBLIC :: init_catalyst, sample_catalyst, finish_catalyst
CONTAINS
    SUBROUTINE catalyst_adaptor_initialize()
        ! Local variables
        TYPE(C_PTR) node, script_node
        INTEGER(kind(catalyst_status)) :: err
        INTEGER(intk) :: iscript, ndummyfields, ilvl, igriddummy, kk, jj, ii, idummyfield, nval, lvlindex
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
                    scriptname//"/filename", catascripts(iscript))
            END DO
        END IF

        err = c_catalyst_initialize(node)
        IF (err /= catalyst_status_ok) THEN
            WRITE (stderr, *) "ERROR: Failed to initialize Catalyst: ", err
        END IF

        ! Allocate dummy fields
        ! If Catalyst does not see a grid on lvl=0 in an MPI rank, then it will segfault in the conversion to a vtkOverlappingAMR
        ! We can avoid this issue by adding a dummy grid which holds a dummy field to "occupy" the grid in the domain
        ! with some grid dimensions that exist in another MPI rank.
        ALLOCATE(dummy_fields(maxlevel + 1))
        DO ilvl = minlevel, maxlevel
            lvlindex = ilvl - minlevel + 1
            IF (nmygridslvl(ilvl) == 0) THEN
                igriddummy = igrdoflevel(1, ilvl)
                CALL get_mgdims(kk, jj, ii, igriddummy)
                nval = (ii - 4) * (jj - 4) * (kk - 4)
                dummy_fields(lvlindex)%nval = nval
                ALLOCATE(dummy_fields(lvlindex)%values(nval))
                dummy_fields(lvlindex)%values = -1
            END IF
        END DO

        CALL catalyst_conduit_node_destroy(node)
    END SUBROUTINE catalyst_adaptor_initialize

    SUBROUTINE add_dummy_grid(data_node, ilvl, itstep, timeph)
        ! Subroutine arguments
        TYPE(C_PTR), INTENT(INOUT) :: data_node
        INTEGER(intk), INTENT(IN) :: ilvl, itstep
        REAL(realk), INTENT(IN) :: timeph

        ! Local variables
        TYPE(C_PTR) :: grid_node, coords_node, field_node, topo_node
        INTEGER(intk) :: igrid_dummy, kk, jj, ii, lvlindex
        INTEGER(kind=8) :: nval
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, dx, dy, dz
        CHARACTER(len=22) :: grid_name

        ! Get first grid on any rank for ilvl
        igrid_dummy = igrdoflevel(1, ilvl)
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid_dummy)
        CALL get_mgdims(kk, jj, ii, igrid_dummy)

        CALL igrid_to_string(igrid_dummy, myid, grid_name)
        grid_node = catalyst_conduit_node_fetch(data_node, "dummygridhello")
        CALL catalyst_conduit_node_set_path_int32(grid_node, &
            "state/domain_id", igrid_dummy + (1+myid)*100)
        CALL catalyst_conduit_node_set_path_int32(grid_node, &
            "state/cycle", itstep)
        CALL catalyst_conduit_node_set_path_float32(grid_node, &
            "state/time", timeph)
        CALL catalyst_conduit_node_set_path_int32(grid_node, &
            "state/level", ilvl - minlevel)

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
        dx = round_to_n_decimals((maxx - minx) / (ii - 4), 7)
        dy = round_to_n_decimals((maxy - miny) / (jj - 4), 7)
        dz = round_to_n_decimals((maxz - minz) / (kk - 4), 7)
        CALL catalyst_conduit_node_set_path_float32(coords_node, &
            "spacing/dx", dx)
        CALL catalyst_conduit_node_set_path_float32(coords_node, &
            "spacing/dy", dy)
        CALL catalyst_conduit_node_set_path_float32(coords_node, &
            "spacing/dz", dz)

        ! Grid topology
        topo_node = catalyst_conduit_node_fetch(grid_node, &
            "topologies/topo")
        CALL catalyst_conduit_node_set_path_char8_str(topo_node, &
            "type", "uniform")
        CALL catalyst_conduit_node_set_path_char8_str(topo_node, &
            "coordset", "coords")
        CALL catalyst_conduit_node_set_path_int32(topo_node, &
            "elements/origin/i0", 0)
        CALL catalyst_conduit_node_set_path_int32(topo_node, &
            "elements/origin/j0", 0)
        CALL catalyst_conduit_node_set_path_int32(topo_node, &
            "elements/origin/k0", 0)

        ! Grid fields
        field_node = catalyst_conduit_node_fetch(grid_node, "fields/DUMMY")
        CALL catalyst_conduit_node_set_path_char8_str(field_node, &
            "association", "element")
        CALL catalyst_conduit_node_set_path_char8_str(field_node, &
            "topology", "topo")
        CALL catalyst_conduit_node_set_path_char8_str(field_node, &
            "volume_dependent", "false")
        lvlindex = ilvl - minlevel + 1
        CALL catalyst_conduit_node_set_path_external_float32_ptr( &
            field_node, "values", dummy_fields(lvlindex)%values, dummy_fields(lvlindex)%nval)

        !CALL catalyst_conduit_node_print(data_node)
    END SUBROUTINE add_dummy_grid

    SUBROUTINE catalyst_adaptor_execute(itstep, ittot, timeph, dt)
        USE grids_mod, ONLY: minlevel, maxlevel, nmygridslvl, get_children
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(C_PTR) exec_node, channel_node, data_node, grid_node, &
            coords_node, topo_node, fields_node, field_node, parent_node, &
            nest_node, child_node
        INTEGER(kind(catalyst_status)) :: err
        INTEGER(intk) :: ilvl, igridlvl, igrid, kk, jj, ii, shiftedlvl, &
            ifield, parentid
        INTEGER(intk) :: parent_kk, parent_jj, parent_ii, child_kk, child_jj, child_ii
        INTEGER(kind=8) :: nval
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, dx, dy, dz
        REAL(realk) :: parent_minx, parent_maxx, parent_miny, parent_maxy, &
            parent_minz, parent_maxz, parent_dx, parent_dy, parent_dz
        REAL(realk) :: child_minx, child_maxx, child_miny, child_maxy, &
            child_minz, child_maxz, child_dx, child_dy, child_dz
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: field_ptr
        CHARACTER(len=22) :: grid_name
        CHARACTER(len=15) :: window_name
        TYPE(field_t), POINTER :: field

        INTEGER, POINTER :: childarray(:)
        INTEGER :: nchild, igrdchild, ichild

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
            IF (nmygridslvl(ilvl) == 0) THEN
                !CALL add_dummy_grid(data_node, ilvl, itstep, timeph)
            END IF

            DO igridlvl = 1, nmygridslvl(ilvl)
                igrid = mygridslvl(igridlvl, ilvl)
                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

                ! Grids state
                CALL igrid_to_string(igrid, myid, grid_name)
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
                dx = round_to_n_decimals((maxx - minx) / (ii - 4), 7)
                dy = round_to_n_decimals((maxy - miny) / (jj - 4), 7)
                dz = round_to_n_decimals((maxz - minz) / (kk - 4), 7)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "spacing/dx", dx)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "spacing/dy", dy)
                CALL catalyst_conduit_node_set_path_float32(coords_node, &
                    "spacing/dz", dz)
            
                ! Grid topology
                topo_node = catalyst_conduit_node_fetch(grid_node, &
                    "topologies/topo")
                CALL catalyst_conduit_node_set_path_char8_str(topo_node, &
                    "type", "uniform")
                CALL catalyst_conduit_node_set_path_char8_str(topo_node, &
                    "coordset", "coords")
                !CALL catalyst_conduit_node_set_path_int32(topo_node, &
                !    "elements/origin/i0", 0)
                !CALL catalyst_conduit_node_set_path_int32(topo_node, &
                !    "elements/origin/j0", 0)
                !CALL catalyst_conduit_node_set_path_int32(topo_node, &
                !    "elements/origin/k0", 0)

                ! Comment out nestsets, this implementation provides the skeleton for child-parent relationships
                ! However, as per https://gitlab.kitware.com/vtk/vtk/-/blob/master/IO/CatalystConduit/vtkConduitToDataObject.cxx
                ! nestset relationships don't actually do anything in Catalyst for an amrmesh
                ! Nestsets
                !nest_node = catalyst_conduit_node_fetch(grid_node, &
                !    "nestsets/nest")
                !CALL catalyst_conduit_node_set_path_char8_str(nest_node, &
                !    "association", "element")
                !CALL catalyst_conduit_node_set_path_char8_str(nest_node, &
                !    "topology", "topo")
                !! Grids with a parent
                !IF (ilvl > minlevel) THEN
                !    parentid = iparent(igrid)
                !    CALL get_window_string(parentid, window_name)
                !    parent_node = catalyst_conduit_node_fetch(nest_node, &
                !        "windows/"//window_name)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "domain_id", parentid)
                !    CALL catalyst_conduit_node_set_path_char8_str(parent_node, &
                !        "domain_type", "parent")
                !    
                !    CALL get_mgdims(parent_kk, parent_jj, parent_ii, parentid)
                !    CALL get_bbox(parent_minx, parent_maxx, parent_miny, parent_maxy, parent_minz, parent_maxz, parentid)
                !    parent_dx = round_to_n_decimals((parent_maxx - parent_minx) / (parent_ii - 4), 7)
                !    parent_dy = round_to_n_decimals((parent_maxy - parent_miny) / (parent_jj - 4), 7)
                !    parent_dz = round_to_n_decimals((parent_maxz - parent_minz) / (parent_kk - 4), 7)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "origin/i", 0)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "origin/j", 0)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "origin/k", 0)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "dims/i", 0)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "dims/j", 0)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "dims/k", 0)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "ratio/i", 2)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "ratio/j", 2)
                !    CALL catalyst_conduit_node_set_path_int32(parent_node, &
                !        "ratio/k", 2)
                !END IF
                !! Grids with at least one child
                !IF (ilvl < maxlevel) THEN
                !    CALL get_children(igrid, childarray, nchild)
                !    DO igrdchild = 1, nchild
                !        ichild = childarray(igrdchild)
                !        CALL get_window_string(ichild, window_name)
                !        child_node = catalyst_conduit_node_fetch(nest_node, &
                !            "windows/"//window_name)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !            "domain_id", ichild)
                !        CALL catalyst_conduit_node_set_path_char8_str(child_node, &
                !            "domain_type", "child")
                !        
                !        CALL get_mgdims(child_kk, child_jj, child_ii, ichild)
                !        CALL get_bbox(child_minx, child_maxx, child_miny, child_maxy, child_minz, child_maxz, ichild)
                !        child_dx = round_to_n_decimals((child_maxx - child_minx) / (child_ii - 4), 7)
                !        child_dy = round_to_n_decimals((child_maxy - child_miny) / (child_jj - 4), 7)
                !        child_dz = round_to_n_decimals((child_maxz - child_minz) / (child_kk - 4), 7)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "origin/i", 0)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "origin/j", 0)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "origin/k", 0)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "dims/i", 0)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "dims/j", 0)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "dims/k", 0)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "ratio/i", 2)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "ratio/j", 2)
                !        CALL catalyst_conduit_node_set_path_int32(child_node, &
                !                "ratio/k", 2)
                !    END DO
                !END IF
                !                 Grid fields
                fields_node = catalyst_conduit_node_fetch(grid_node, "fields")
                nval = (ii - 4) * (jj - 4) * (kk - 4)
                DO ifield = 1, nfields
                    field_node = catalyst_conduit_node_fetch(fields_node, &
                        TRIM(fields(ifield)))
                    CALL get_field(field, TRIM(fields(ifield))//"_C")
                    CALL field%get_ptr(field_ptr, igrid)
                    CALL catalyst_conduit_node_set_path_char8_str(field_node, &
                        "association", "element")
                    CALL catalyst_conduit_node_set_path_char8_str(field_node, &
                        "topology", "topo")
                    CALL catalyst_conduit_node_set_path_char8_str(field_node, &
                        "volume_dependent", "false")
                    CALL catalyst_conduit_node_set_path_external_float32_ptr( &
                        field_node, "values", field_ptr, nval)
                END DO
            END DO
        END DO

        IF (myid == 0) THEN
            !CALL catalyst_conduit_node_print(exec_node)
        END IF

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

        IF (ALLOCATED(dummy_fields)) DEALLOCATE(dummy_fields)
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
            WRITE(*,*) "Specifiy directory containing libcatalyst-paraview.so!"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Representative dataset
        IF ( cata_conf%is_logical("/repr") ) THEN
            CALL cata_conf%get_value("/repr", is_repr, .FALSE.)
        ELSE
            WRITE(*,*) "Specify whether to create a representative dataset!"
            is_repr = .FALSE.
        END IF

        ! Fields
        IF (cata_conf%exists("/fields")) THEN
            CALL cata_conf%get_size("/fields", nfields)
            ALLOCATE(fields(nfields))
            DO i = 1, nfields
                WRITE(jsonptr, '("/fields/", I0)') i-1
                CALL cata_conf%get_value(jsonptr, fields(i))
            END DO
        ELSE
            WRITE(*,*) "Fields for Catalyst not specified!"
            CALL errr(__FILE__, __LINE__)
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
        DO i = 1, nfields
            CALL set_field(TRIM(fields(i))//"_C")
        END DO

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

    SUBROUTINE igrid_to_string(igrid, irank, output_string)
        INTEGER, INTENT(in)              :: igrid, irank
        CHARACTER(len=22), INTENT(out)   :: output_string
        CHARACTER(len=8)                 :: temp_string_igrid, temp_string_myid
        WRITE(temp_string_igrid, '(I8.8)') igrid
        WRITE(temp_string_myid,  '(I8.8)') irank
        output_string = 'grid_' // temp_string_igrid // '_' // temp_string_myid
    END SUBROUTINE igrid_to_string

    SUBROUTINE get_window_string(parentid, output_string)
        INTEGER, INTENT(in) :: parentid
        CHARACTER(len=15), INTENT(out) :: output_string
        CHARACTER(len=8) :: temp_string
        WRITE(temp_string, '(I8.8)') parentid
        output_string = 'window_'//temp_string
    END SUBROUTINE get_window_string

    SUBROUTINE iscript_to_string(iscript, output_string)
        INTEGER, INTENT(in) :: iscript
        CHARACTER(len=10), INTENT(out) :: output_string
        CHARACTER(len=3) :: temp_string
        WRITE(temp_string, '(I3.3)') iscript
        output_string = 'script_'//temp_string
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
        INTEGER(intk) :: i, igrid, ifield
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: f_t, f_t_c
        REAL(realk), POINTER, CONTIGUOUS :: arr_c(:, :, :), arr(:, :, :)
        INTEGER, PARAMETER :: nbl = 2

        DO ifield = 1, nfields
            CALL get_field(f_t, TRIM(fields(ifield)))
            CALL get_field(f_t_c, TRIM(fields(ifield))//"_C")
            DO i = 1, nmygrids
                igrid = mygrids(i)
                CALL get_mgdims(kk, jj, ii, igrid)
                CALL f_t%get_ptr(arr, igrid)
                CALL f_t_c%get_ptr(arr_c, igrid)
                CALL field_to_c_grid(kk, jj, ii, nbl, arr_c, arr)
            END DO
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
