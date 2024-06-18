MODULE openmp_mod
    use omp_lib

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_openmp

CONTAINS
    SUBROUTINE init_openmp()
        USE comms_mod, ONLY: myid

        INTEGER :: num_devices
        num_devices = omp_get_num_devices()

        IF (myid == 0) THEN
            WRITE(*, '("OPENMP INFORMATION:")')
            WRITE(*, '("    Target devices:     ", I0)') num_devices
            WRITE(*, '()')
        END IF

    END SUBROUTINE init_openmp

END MODULE openmp_mod
