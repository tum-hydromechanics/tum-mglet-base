message(STATUS "Checking if OpenMP threading works")

if ( THREADS )

    find_package(OpenMP QUIET)
    set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_FLAGS})

    check_fortran_source_runs("
    PROGRAM main
        USE omp_lib
        IMPLICIT NONE
        INTEGER :: n

        CALL omp_set_num_threads(2)

        !$omp parallel private(n)
        n = omp_get_num_threads()
        IF ( n /= 2 ) THEN
            STOP 9
        END IF
        !$omp end parallel

    END PROGRAM" ASSUMED_OPENMP_OK SRC_EXT "F90")
    if(NOT ASSUMED_OPENMP_OK)
        message(WARNING "Your Fortran compiler does not support OpenMP")
    endif()

endif()
