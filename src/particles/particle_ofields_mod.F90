MODULE particle_ofields_mod
    USE precision_mod, ONLY: intk, realk
    USE pointers_mod, ONLY: ip3d, ip1d
    USE grids_mod, ONLY: nmygrids, get_mgdims, get_mgbasb, get_bbox, get_bc_ctyp
    USE err_mod, ONLY: errr
    USE fields_mod
    USE realfield_mod
    USE flowcore_mod
    USE ib_mod

    USE particle_config_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Module constants
    INTEGER(intk), PARAMETER :: N_DIMS = 3
    INTEGER(intk), PARAMETER :: N_BASB = 6
    INTEGER(intk), PARAMETER :: N_FACES = 6
    INTEGER(intk), PARAMETER :: N_BC_RANGE = 2

    ! ┌────────────────────────────────────────────────────────────────────────────┐
    ! | Keeps a pointer to the data that is required on the target device          |
    ! | WHY POINTERS?                                                              |
    ! |     - Keep omp directives central for the sake of code readability         |
    ! |     - Prevent any unwanted intereference with the core flow implementation |
    ! |     - Allows to directly map field data without omp directives in fields   |
    ! └────────────────────────────────────────────────────────────────────────────┘
    
    ! ----- Pointers to fields -----
    ! Grid parameters
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:) :: ip3d_offload, ip1d_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x_offload, y_offload, z_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx_offload, dy_offload, dz_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx_offload, ddy_offload, ddz_offload
    
    ! Flow/Scalar fields
    REAL(realk), ALLOCATABLE, TARGET :: unull(:), vnull(:), wnull(:)
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: u_offload, v_offload, w_offload
    
    ! ----- Newly encoded or global arrays -----
    REAL(realk), ALLOCATABLE :: bbox_offload(:)
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:) :: mgdims_offload

#ifdef __GFORTRAN__
    !$omp declare target(ip3d_offload, ip1d_offload, mgdims_offload, bbox_offload)
    !$omp declare target(x_offload, y_offload, z_offload)
    !$omp declare target(dx_offload, dy_offload, dz_offload, ddx_offload, ddy_offload, ddz_offload)
    !$omp declare target(u_offload, v_offload, w_offload)
#endif

    ! Public subroutines for host
    PUBLIC :: offload_fields, finish_offload_fields

    ! Public variables for host

    ! Public subroutines for device
    PUBLIC :: ptr_to_grid_x, ptr_to_grid_y, ptr_to_grid_z, ptr_to_grid3, get_bbox_target, get_mgdims_target

    ! Public variables for device
    PUBLIC :: x_offload, y_offload, z_offload, &
        dx_offload, dy_offload, dz_offload, ddx_offload, ddy_offload, ddz_offload, &
        u_offload, v_offload, w_offload, bbox_offload, unull, vnull, wnull

CONTAINS

    !> @brief Sets up field pointers for target device
    !!
    !! All field pointers required by the scalar computation are made available for use on the target device.
    !! Includes both setting up a new pointer to fields and maps the data to the target device.
    SUBROUTINE offload_fields()        
        CALL map_grid_data()
        CALL map_constant_grid_fields()
        
        ! particle boundary conditions are handled by thee particle_boundaries_mod
        ! (for particles, flow/sclalar boundary information is only needed at initialization)

        CALL map_flow()
    END SUBROUTINE offload_fields


    !> @brief Sets up grid data for target device
    !!
    !! Sets up a custom data representation for mgdims suitable for mapping to the target device.
    !! Sets up pointers to ip3d and ip1d.
    !! Maps fields to the target device.
    SUBROUTINE map_grid_data()

        ! Local variables
        INTEGER(intk) :: igrid, i, mgdims_arr_size, bbox_arr_size, kk, jj, ii
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        ! Create grids_mod copy to offload
        mgdims_arr_size = 3 * ngrid
        ALLOCATE(mgdims_offload(mgdims_arr_size))
        DO igrid = 1, ngrid
            i = (igrid - 1) * 3 + 1
            CALL get_mgdims(kk, jj, ii, igrid)
            mgdims_offload(i) = ii
            mgdims_offload(i+1) = jj
            mgdims_offload(i+2) = kk
        END DO

        ! Create grids_mod copy to offload
        bbox_arr_size = 6 * ngrid
        ALLOCATE(bbox_offload(bbox_arr_size))
        DO igrid = 1, ngrid
            i = (igrid - 1) * 6 + 1
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
            bbox_offload(i) = minx
            bbox_offload(i+1) = maxx
            bbox_offload(i+2) = miny
            bbox_offload(i+3) = maxy
            bbox_offload(i+4) = minz
            bbox_offload(i+5) = maxz
        END DO

        ! Create pointers to pointers_mod fields just to have all omp directives to map data in this file
        ip3d_offload => ip3d
        ip1d_offload => ip1d
        
        !$omp target enter data map(to: ip3d_offload, ip1d_offload, mgdims_offload, bbox_offload)
    END SUBROUTINE


    !> @brief Sets up constant grid fields for target device
    !!
    !! Sets up pointers to constant grid fields (e.g. spacings or derivatives in spatial directions).
    !! Maps fields to the target device.
    SUBROUTINE map_constant_grid_fields()
        ! Local variables
        TYPE(field_t), POINTER :: x_f, y_f, z_f, dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f, bt_f

        ! Create copy for grid constants
        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        x_offload => x_f%arr
        y_offload => y_f%arr
        z_offload => z_f%arr
        dx_offload => dx_f%arr
        dy_offload => dy_f%arr
        dz_offload => dz_f%arr
        ddx_offload => ddx_f%arr
        ddy_offload => ddy_f%arr
        ddz_offload => ddz_f%arr

        !$omp target enter data map(to: x_offload, y_offload, z_offload)
        !$omp target enter data map(to: dx_offload, dy_offload, dz_offload, ddx_offload, ddy_offload, ddz_offload)
    END SUBROUTINE


    !> @brief Sets up flow fields for target device
    !!
    !! Sets up pointers to flow fields U, V, W
    !! Maps fields to the target device.
    SUBROUTINE map_flow()
        ! Local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f, sca_f, g_f

        IF (dadvection) THEN
            IF (duse_avg_flow) THEN
                ! use the point values deduced from the average flow field
                CALL get_field(u_f, "PWU_AVG")
                CALL get_field(v_f, "PWV_AVG")
                CALL get_field(w_f, "PWW_AVG")
            ELSE
                IF (ib%type == "GHOSTCELL") THEN
                    CALL get_field(u_f, "PWU")
                    CALL get_field(v_f, "PWV")
                    CALL get_field(w_f, "PWW")
                ELSE
                    CALL get_field(u_f, "U")
                    CALL get_field(v_f, "V")
                    CALL get_field(w_f, "W")
                END IF
            END IF
            u_offload => u_f%arr
            v_offload => v_f%arr
            w_offload => w_f%arr
        ELSE 
            ! TODO: remove this dirty workaround
            ALLOCATE(unull(1))
            ALLOCATE(vnull(1))
            ALLOCATE(wnull(1))
            u_offload => unull
            v_offload => vnull
            w_offload => wnull
        END IF
        
        !$omp target enter data map(to: u_offload, v_offload, w_offload)
    END SUBROUTINE


    !> @brief Finish offloaded fields
    !!
    !! All newly allocated field pointers required by the scalar computation are deallocated.
    !! Exits all fields on the target device
    SUBROUTINE finish_offload_fields()
        !$omp target exit data map(delete: bbox_offload, mgdims_offload, ip3d_offload, ip1d_offload)
        !$omp target exit data map(delete: x_offload, y_offload, z_offload)
        !$omp target exit data map(delete: dx_offload, dy_offload, dz_offload, ddx_offload, ddy_offload, ddz_offload)
        !$omp target exit data map(delete: u_offload, v_offload, w_offload)

        DEALLOCATE(mgdims_offload)
        DEALLOCATE(bbox_offload)
    END SUBROUTINE finish_offload_fields


    !> @brief Get mgdims - grid dimensions
    !!
    !! Mirrors core get_mgdims method but uses fields available on the target device
    !! @param[out] kk    Z-direction dimensions
    !! @param[out] jj    Y-direction dimensions
    !! @param[out] ii    X-direction dimensions
    !! @param[in]  igrid Grid index
    SUBROUTINE get_mgdims_target(kk, jj, ii, igrid)
        
        !$omp declare target
        
        INTEGER(intk), INTENT(OUT) :: kk, jj, ii
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        ii = mgdims_offload(i)
        jj = mgdims_offload(i+1)
        kk = mgdims_offload(i+2)
    END SUBROUTINE get_mgdims_target


    SUBROUTINE get_bbox_target(minx, maxx, miny, maxy, minz, maxz, igrid)

        !$omp declare target

        REAL(realk), INTENT(OUT) :: minx, maxx, miny, maxy, minz, maxz
        INTEGER(intk), INTENT(IN) :: igrid

        ! local variables
        INTEGER(intk) :: i

        i = (igrid - 1) * 6 + 1
        minx = bbox_offload(i)
        maxx = bbox_offload(i+1)
        miny = bbox_offload(i+2)
        maxy = bbox_offload(i+3)
        minz = bbox_offload(i+4)
        maxz = bbox_offload(i+5)
    END SUBROUTINE get_bbox_target


    !> @brief Get field length 1d in x-direction
    !!
    !! Mirrors core get_len method but uses fields available on the target device.
    !! @param[out] len   Field length
    !! @param[in]  igrid Grid index
    SUBROUTINE get_len_ii_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i)
    END SUBROUTINE get_len_ii_target


    !> @brief Get field length 1d in y-direction
    !!
    !! Mirrors core get_len method but uses fields available on the target device.
    !! @param[out] len   Field length
    !! @param[in]  igrid Grid index
    SUBROUTINE get_len_jj_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i+1)
    END SUBROUTINE get_len_jj_target


    !> @brief Get field length 1d in z-direction
    !!
    !! Mirrors core get_len method but uses fields available on the target device.
    !! @param[out] len   Field length
    !! @param[in]  igrid Grid index
    SUBROUTINE get_len_kk_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i+2)
    END SUBROUTINE get_len_kk_target


    !> @brief Get ip1 (1 dimensional fields)
    !!
    !! Mirrors core get_ip1 method but uses fields available on the target device.
    !! @param[out] ip1   Index for pointer to nth grid
    !! @param[in]  igrid Grid index
    SUBROUTINE get_ip1_target(ip1, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip1

        ip1 = ip1d_offload(igrid)
    END SUBROUTINE get_ip1_target


    !> @brief Get ip3 (3 dimensional fields)
    !!
    !! Mirrors core get_ip3 method but uses fields available on the target device.
    !! @param[out] ip3   Index for pointer to nth grid
    !! @param[in]  igrid Grid index
    SUBROUTINE get_ip3_target(ip3, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip3

        ip3 = ip3d_offload(igrid)
    END SUBROUTINE get_ip3_target


    !> @brief Get pointer to 1d grid-specific field (only x direction)
    !!
    !! Mirrors core get_grid1 method but uses fields available on the target device.
    !! Apply only for x-direction.
    !! @param[in]     arr_ptr  Pointer to field data
    !! @param[in]     n_grid   Grid index
    !! @param[in,out] grid_ptr 1d-pointer to grid-specific field data
    SUBROUTINE ptr_to_grid_x(arr_ptr, n_grid, grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(inout) :: grid_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid

        ! Local variables
        INTEGER(intk) :: ip, len
        
        CALL get_ip1_target(ip, n_grid)
        CALL get_len_ii_target(len, n_grid)
        grid_ptr(1:len) => arr_ptr(ip:ip+len-1)
    END SUBROUTINE ptr_to_grid_x


    !> @brief Get pointer to 1d grid-specific field (only y direction)
    !!
    !! Mirrors core get_grid1 method but uses fields available on the target device.
    !! Apply only for y-direction.
    !! @param[in]     arr_ptr  Pointer to field data
    !! @param[in]     n_grid   Grid index
    !! @param[in,out] grid_ptr 1d-pointer to grid-specific field data
    SUBROUTINE ptr_to_grid_y(arr_ptr, n_grid, grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(inout) :: grid_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid

        ! Local variables
        INTEGER(intk) :: ip, len
        
        CALL get_ip1_target(ip, n_grid)
        CALL get_len_jj_target(len, n_grid)
        grid_ptr(1:len) => arr_ptr(ip:ip+len-1)
    END SUBROUTINE ptr_to_grid_y


    !> @brief Get pointer to 1d grid-specific field (only z direction)
    !!
    !! Mirrors core get_grid1 method but uses fields available on the target device.
    !! Apply only for z-direction.
    !! @param[in]     arr_ptr  Pointer to field data
    !! @param[in]     n_grid   Grid index
    !! @param[in,out] grid_ptr 1d-pointer to grid-specific field data
    SUBROUTINE ptr_to_grid_z(arr_ptr, n_grid, grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(inout) :: grid_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid
        ! Result variables
        ! Local variables
        INTEGER(intk) :: ip, len
        
        CALL get_ip1_target(ip, n_grid)
        CALL get_len_kk_target(len, n_grid)
        grid_ptr(1:len) => arr_ptr(ip:ip+len-1)
    END SUBROUTINE ptr_to_grid_z


    !> @brief Get pointer to 3d grid-specific field (only x direction)
    !!
    !! Mirrors core get_grid3 method but uses fields available on the target device.
    !! @param[in]     arr_ptr  Pointer to field data
    !! @param[in]     n_grid   Grid index
    !! @param[in,out] grid_ptr 3d-pointer to grid-specific field data
    SUBROUTINE ptr_to_grid3(arr_ptr, n_grid, grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(inout) :: grid_ptr(:, :, :)
        INTEGER(intk), INTENT(in) :: n_grid
        ! Result variables
        ! Local variables
        INTEGER(intk) :: ip, ii, jj, kk
        
        CALL get_mgdims_target(kk, jj, ii, n_grid)
        CALL get_ip3_target(ip, n_grid)
        grid_ptr(1:kk, 1:jj, 1:ii) => arr_ptr(ip:ip+kk*jj*ii-1)
    END SUBROUTINE ptr_to_grid3

END MODULE particle_ofields_mod
