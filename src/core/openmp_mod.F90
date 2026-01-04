MODULE openmp_mod
    USE precision_mod
    USE omp_lib

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: init_openmp

CONTAINS

    SUBROUTINE init_openmp()
        USE comms_mod, ONLY: myid

        INTEGER(intk) :: num_devices, host_device_nbr, device_nbr

        num_devices = omp_get_num_devices()

        host_device_nbr = omp_get_initial_device()
        device_nbr = omp_get_device_num()
        IF (myid == 0) THEN
            WRITE(*, '("OPENMP INFORMATION:")')
            WRITE(*, '("    Target devices:     ", I0)') num_devices
            WRITE(*, '("    Host device:        ", I0)') host_device_nbr
            WRITE(*, '()')
            WRITE(*, '("    Current device:        ", I0)') device_nbr
            WRITE(*, '()')
        END IF

    END SUBROUTINE init_openmp

END MODULE openmp_mod
