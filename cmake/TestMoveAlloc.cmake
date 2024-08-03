message(STATUS "Checking if your Fortran compiler supports 'MOVE_ALLOC'")
check_fortran_source_runs("
PROGRAM test_move_alloc
    integer, allocatable :: a(:), b(:)
    allocate(a(3))
    allocate(b(5))
    a = [ 1, 2, 3 ]
    b = 0
    b(1:3)= a
    call move_alloc(b,a)
END PROGRAM" MOVE_ALLOC_OK SRC_EXT "F90")
if(NOT MOVE_ALLOC_OK)
    message(WARNING "Your Fortran compiler does not support 'MOVE_ALLOC'")
endif()
