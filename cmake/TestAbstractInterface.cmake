message(STATUS "Checking if Fortran compiler supports abstract interface")
check_fortran_source_runs("
MODULE abstract_mod

    ABSTRACT INTERFACE
        SUBROUTINE add(a,b,c)
            INTEGER, INTENT(in) :: a, b
            INTEGER, INTENT(out) :: c
        END SUBROUTINE
    END INTERFACE

END MODULE abstract_mod" ABSTRACT_TEST_OK SRC_EXT "F90")
if(NOT ABSTRACT_TEST_OK)
    message(FATAL_ERROR "Compiler fails at 'abstract interface'")
endif()
