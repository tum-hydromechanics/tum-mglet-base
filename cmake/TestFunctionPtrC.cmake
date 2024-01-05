message(STATUS "Checking if Fortran compiler supports C function pointers")
check_fortran_source_runs("
MODULE ptr_mod

    abstract interface
        subroutine add_f_type(a, b, c)
            implicit none
            integer, intent(inout) :: a, b, c
        end subroutine add_f_type
    end interface

contains

    subroutine ptr_conversion()
        use, intrinsic :: iso_c_binding, only : c_funloc, c_funptr
        procedure(add_f_type) :: add_f
        type(c_funptr) :: add_c_fptr
        add_c_fptr = c_funloc(add_f)
    end subroutine ptr_conversion

END MODULE ptr_mod" CPTR_TEST_OK SRC_EXT "F90")
if(NOT CPTR_TEST_OK)
    message(FATAL_ERROR "Compiler fails at C function pointers")
endif()
