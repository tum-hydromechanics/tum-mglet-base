MODULE offload_helper_mod
    USE precision_mod, ONLY: intk, realk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: N_CELLS_PER_DIM = 36
    INTEGER(intk), PARAMETER :: N_DIMS = 3

    INTEGER(intk), POINTER, CONTIGUOUS :: ip3d_offload(:), ip1d_offload(:)
    INTEGER(intk), POINTER, CONTIGUOUS :: mgdims_offload(:)
    REAL(realk), POINTER, CONTIGUOUS :: rddx_offload(:), rddy_offload(:), rddz_offload(:)

    !$omp declare target(ip3d_offload, ip1d_offload, mgdims_offload, rddx_offload, rddy_offload, rddz_offload)

    PUBLIC :: N_CELLS_PER_DIM, ptr_to_grid_x, ptr_to_grid_y, ptr_to_grid_z, ptr_to_grid3, offload_constants, finish_offload_constants, rddx_offload, rddy_offload, rddz_offload, get_mgdims_target
CONTAINS
    SUBROUTINE offload_constants()
        USE pointers_mod, ONLY: ip3d, ip1d
        USE grids_mod, ONLY: nmygrids, get_mgdims
        USE fields_mod
        USE realfield_mod

        ! Local variablees
        INTEGER(intk) :: igrid, i, mgdims_arr_size, kk, jj, ii
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f

        ! Create grids_mod copy to offload
        mgdims_arr_size = N_DIMS * nmygrids
        ALLOCATE(mgdims_offload(mgdims_arr_size))
        DO igrid = 1, nmygrids
            i = (igrid - 1) * N_DIMS + 1
            CALL get_mgdims(kk, jj, ii, igrid)
            mgdims_offload(i) = ii
            mgdims_offload(i+1) = jj
            mgdims_offload(i+2) = kk
        END DO
        print *, mgdims_offload

        ! Create pointers_mod copy to offload
        ALLOCATE(ip3d_offload(SIZE(ip3d)))
        ALLOCATE(ip1d_offload(SIZE(ip1d)))
        ip3d_offload = ip3d
        ip1d_offload = ip1d

        ! Create copy for grid constants
        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")
        ALLOCATE(rddx_offload(SIZE(rddx_f%arr)))
        ALLOCATE(rddy_offload(SIZE(rddy_f%arr)))
        ALLOCATE(rddz_offload(SIZE(rddz_f%arr)))

        rddx_offload = rddx_f%arr
        rddy_offload = rddy_f%arr
        rddz_offload = rddz_f%arr

        !$omp target enter data map(to: ip3d_offload, ip1d_offload, mgdims_offload, rddx_offload, rddy_offload, rddz_offload)
    END SUBROUTINE offload_constants

    SUBROUTINE finish_offload_constants()
        !$omp target exit data map(delete: ip3d_offload, ip1d_offload, mgdims_offload, rddx_offload, rddy_offload, rddz_offload)

        DEALLOCATE(ip3d_offload)
        DEALLOCATE(ip1d_offload)
        DEALLOCATE(mgdims_offload)
        DEALLOCATE(rddx_offload)
        DEALLOCATE(rddy_offload)
        DEALLOCATE(rddz_offload)
    END SUBROUTINE finish_offload_constants

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

    SUBROUTINE get_len_ii_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i)
    END SUBROUTINE get_len_ii_target

    SUBROUTINE get_len_jj_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i+1)
    END SUBROUTINE get_len_jj_target

    SUBROUTINE get_len_kk_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i+2)
    END SUBROUTINE get_len_kk_target

    SUBROUTINE get_ip1_target(ip1, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip1

        ip1 = ip1d_offload(igrid)
    END SUBROUTINE get_ip1_target

    SUBROUTINE get_ip3_target(ip3, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip3

        ip3 = ip3d_offload(igrid)
    END SUBROUTINE get_ip3_target

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
END MODULE offload_helper_mod
