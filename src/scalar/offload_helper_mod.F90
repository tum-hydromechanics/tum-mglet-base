MODULE offload_helper_mod
    USE precision_mod, ONLY: intk, realk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: N_CELLS_PER_DIM = 36

    INTEGER(intk), ALLOCATABLE :: ip3d_offload(:), ip1d_offload(:)
    INTEGER(intk), ALLOCATABLE :: mgdims_offload(:)

    !$omp declare target(ip3d_offload, ip1d_offload, mgdims_offload)

    PUBLIC :: N_CELLS_PER_DIM, ptr_to_grid1, ptr_to_grid3, offload_constants, release_constants
CONTAINS
    SUBROUTINE offload_constants()
        USE pointers_mod, ONLY: ip3d, ip1d
        USE grids_mod, ONLY: nmygrids, get_mgdims

        INTEGER(intk) :: i, kk, jj, ii

        ! Make necessary grids_mod copy
        ALLOCATE(mgdims_offload(3*nmygrids))
        DO i = 1, 3*nmygrids, 3
            CALL get_mgdims(kk, jj, ii, i)
            mgdims_offload(i) = ii
            mgdims_offload(i+1) = jj
            mgdims_offload(i+2) = kk
        END DO

        ! Make necessary pointers_mod copies
        ip3d_offload = ip3d
        ip1d_offload = ip1d

        !$omp target enter data map(to: ip3d_offload, ip1d_offload, mgdims_offload)

    END SUBROUTINE offload_constants

    SUBROUTINE release_constants()
        !$omp target exit data map(delete: ip3d_offload, ip1d_offload, mgdims_offload)
        DEALLOCATE(ip3d_offload)
        DEALLOCATE(ip1d_offload)
        DEALLOCATE(mgdims_offload)
    END SUBROUTINE release_constants

    SUBROUTINE get_mgdims_target(kk, jj, ii, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: kk, jj, ii
        INTEGER(intk), INTENT(IN) :: igrid

        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        ii = mgdims_offload(i)
        jj = mgdims_offload(i+1)
        kk = mgdims_offload(i+2)
    END SUBROUTINE get_mgdims_target

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

    FUNCTION ptr_to_grid1(arr_ptr, n_grid) RESULT(grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid
        ! Result variables
        REAL(realk), POINTER, CONTIGUOUS :: grid_ptr(:)
        ! Local variables
        INTEGER(intk) :: ip
        
        CALL get_ip1_target(ip, n_grid)
        grid_ptr(1:N_CELLS_PER_DIM) => arr_ptr(ip:ip+N_CELLS_PER_DIM-1)
    END FUNCTION ptr_to_grid1
    
    FUNCTION ptr_to_grid3(arr_ptr, n_grid) RESULT(grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid
        ! Result variables
        REAL(realk), POINTER, CONTIGUOUS :: grid_ptr(:, :, :)
        ! Local variables
        INTEGER(intk) :: ip, ii, jj, kk
        
        CALL get_mgdims_target(kk, jj, ii, n_grid)
        CALL get_ip3_target(ip, n_grid)
        grid_ptr(1:kk, 1:jj, 1:ii) => arr_ptr(ip:ip+kk*jj*ii-1)
    END FUNCTION ptr_to_grid3
END MODULE offload_helper_mod
