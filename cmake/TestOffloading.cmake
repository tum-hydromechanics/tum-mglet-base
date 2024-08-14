message(STATUS "Checking if Fortran compiler supports GPU offloading")
check_fortran_source_runs("
PROGRAM main
    USE omp_lib

    IMPLICIT NONE(type, external)
    LOGICAL :: is_initial_device
    is_initial_device = .FALSE.

    is_initial_device = omp_is_initial_device()
    IF (is_initial_device) THEN
        WRITE(*, *) 'Program starts with initial device'
    ELSE
        WRITE(*, *) 'Device is not initial device'
        STOP 9
    END IF

    !$omp target map(tofrom:is_initial_device)
    is_initial_device = omp_is_initial_device()
    !$omp end target

    IF (is_initial_device) THEN
        WRITE(*, *) 'No offloaded device used'
        STOP 9
    ELSE
        WRITE(*, *) 'Offloaded device is not initial host device'
    END IF
END PROGRAM" OFFLOADING_OK SRC_EXT "F90")
if(NOT OFFLOADING_OK)
    message(WARNING "GPU offloading is not supported")
    file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" CMAKE_ERROR_LOG)
    message(STATUS "CMake Error Log:\n${CMAKE_ERROR_LOG}")
endif()
