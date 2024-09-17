MODULE offload_helper_mod
    USE precision_mod, ONLY: intk, realk, mglet_hdf5_real, mglet_mpi_real
    USE field_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: N_CELLS_PER_DIM = 36

    PUBLIC :: N_CELLS_PER_DIM, ptr_to_grid1, ptr_to_grid3
CONTAINS
    FUNCTION ptr_to_grid1(ptr_a, n_grid) RESULT(ptr_grid)
        !$omp declare target
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: ptr_a(:)
        INTEGER(intk), INTENT(in) :: n_grid

        REAL(realk), POINTER, CONTIGUOUS :: ptr_grid(:)
        INTEGER(intk) :: ip
        ip = (n_grid - 1) * N_CELLS_PER_DIM + 1

        ptr_grid(1:N_CELLS_PER_DIM) => ptr_a(ip:ip+N_CELLS_PER_DIM-1)
    END FUNCTION ptr_to_grid1

    FUNCTION ptr_to_grid3(ptr_a, n_grid) RESULT(ptr_grid)
        !$omp declare target
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: ptr_a(:)
        INTEGER(intk), INTENT(in) :: n_grid

        REAL(realk), POINTER, CONTIGUOUS :: ptr_grid(:, :, :)
        INTEGER(intk) :: ip
        ip = (n_grid - 1) * N_CELLS_PER_DIM**3 + 1

        ptr_grid(1:N_CELLS_PER_DIM, 1:N_CELLS_PER_DIM, 1:N_CELLS_PER_DIM) => ptr_a(ip:ip+N_CELLS_PER_DIM**3-1)
    END FUNCTION ptr_to_grid3

END MODULE offload_helper_mod
