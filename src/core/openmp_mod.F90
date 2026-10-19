MODULE openmp_mod
    USE precision_mod

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: init_openmp

CONTAINS

    SUBROUTINE init_openmp()
#ifdef _MGLET_OPENMP_
        USE omp_lib
        USE comms_mod, ONLY: myid

        INTEGER(intk) :: num_devices, host_device_nbr, device_nbr

        num_devices = omp_get_num_devices()

        host_device_nbr = omp_get_initial_device()
        device_nbr = omp_get_device_num()
        IF (myid == 0) THEN
            WRITE(*, '("OPENMP INFORMATION:")')
            WRITE(*, '("    Host (ID):                      ", I0)') host_device_nbr
            WRITE(*, '("    Number of Target devices:       ", I0)') num_devices
            WRITE(*, '()')
        END IF
#else
        ! Keep this module link-safe when OpenMP support is disabled.
        CONTINUE
#endif

    END SUBROUTINE init_openmp

END MODULE openmp_mod
