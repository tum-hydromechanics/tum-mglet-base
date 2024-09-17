MODULE offload_helper_mod
    USE precision_mod, ONLY: intk, realk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: N_CELLS_PER_DIM = 36

    PUBLIC :: N_CELLS_PER_DIM, ptr_to_grid1, ptr_to_grid3
CONTAINS
    FUNCTION ptr_to_grid1(arr_ptr, n_grid) RESULT(grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid
        ! Result variables
        REAL(realk), POINTER, CONTIGUOUS :: grid_ptr(:)
        ! Local variables
        INTEGER(intk) :: ip
        
        ip = (n_grid - 1) * N_CELLS_PER_DIM + 1
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
        INTEGER(intk) :: ip
        
        ip = (n_grid - 1) * N_CELLS_PER_DIM**3 + 1
        grid_ptr(1:N_CELLS_PER_DIM, 1:N_CELLS_PER_DIM, 1:N_CELLS_PER_DIM) => arr_ptr(ip:ip+N_CELLS_PER_DIM**3-1)
    END FUNCTION ptr_to_grid3
END MODULE offload_helper_mod
