MODULE offload_realfield_wrapper_mod
    USE precision_mod, ONLY: intk, realk, mglet_hdf5_real, mglet_mpi_real
    USE field_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: N_CELLS_PER_DIM = 36

    TYPE :: offload_realfield
        REAL(realk), POINTER, CONTIGUOUS :: data(:)
    CONTAINS
        GENERIC, PUBLIC :: set_data_ptr => set_data_ptr_from_arr, set_data_ptr_from_field_t
        PROCEDURE, PRIVATE :: set_data_ptr_from_arr, set_data_ptr_from_field_t

        PROCEDURE :: get_grid_data_1d
    END TYPE offload_realfield

    PUBLIC :: offload_realfield

CONTAINS
    SUBROUTINE set_data_ptr_from_field_t(this, field)
        ! Subroutine arguments
        CLASS(offload_realfield), TARGET, INTENT(inout) :: this
        TYPE(field_t), POINTER, INTENT(in) :: field

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: field_arr_ptr(:)

        CALL field%get_arr_ptr(field_arr_ptr)

        this%data => field_arr_ptr
    END SUBROUTINE set_data_ptr_from_field_t

    SUBROUTINE set_data_ptr_from_arr(this, data_ptr)
        ! Subroutine arguments
        CLASS(offload_realfield), TARGET, INTENT(inout) :: this
        REAL(realk), POINTER, INTENT(in) :: data_ptr(:)

        this%data => data_ptr
    END SUBROUTINE set_data_ptr_from_arr

    SUBROUTINE get_grid_data_1d(this, n_grid, grid_data)
        !$omp declare target

        ! Subroutine arguments
        CLASS(offload_realfield), TARGET, INTENT(inout) :: this
        REAL(realk), POINTER, INTENT(inout) :: grid_data(:)
        INTEGER(intk), INTENT(IN) :: n_grid
        
        grid_data(1:N_CELLS_PER_DIM) => this%data((n_grid - 1) * N_CELLS_PER_DIM + 1 : N_CELLS_PER_DIM * n_grid)
    END SUBROUTINE get_grid_data_1d

END MODULE offload_realfield_wrapper_mod
