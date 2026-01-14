MODULE particle_rng_mod
    
    ! see
    ! https://vectrx.substack.com/p/lcg-xs-fast-gpu-rng#footnote-7-159943017
    
    USE, INTRINSIC :: ISO_FORTRAN_ENV
    USE, INTRINSIC :: ISO_C_BINDING

    USE omp_lib

    USE precision_mod

    IMPLICIT NONE

    INTEGER(c_int) :: particle_base_seed = 9891477

    INTEGER(c_int64_t), ALLOCATABLE :: lcg_parameters(:) 
    
    !$omp declare target(lcg_parameters)

CONTAINS

    SUBROUTINE init_parallel_lcg()

        ALLOCATE(lcg_parameters(2))

        !lcg_multiplier
        lcg_parameters(1) = 747796405_c_int64_t
        !lcg_increment
        lcg_parameters(2) = 2891336453_c_int64_t

        !$omp target enter data map(to: lcg_parameters)

    END SUBROUTINE init_parallel_lcg


    SUBROUTINE finish_parallel_lcg()

        !$omp target exit data map(delete: lcg_parameters)

    END SUBROUTINE finish_parallel_lcg


    SUBROUTINE lcg(seed, rn)
        
        !$omp declare target

        ! subroutine arguments
        INTEGER(c_int), INTENT(inout) :: seed
        REAL(c_float), INTENT(out) :: rn

        ! local variables
        INTEGER(c_int64_t) :: tmp
        INTEGER(c_int64_t) :: mask32 = int(Z'FFFFFFFF', c_int64_t)
        REAL(c_float) :: max_int = 16777216.0_c_float

        ! Compute RNG step in 64-bit and mask to 32-bit unsigned
        tmp = iand(int(seed, c_int64_t) * lcg_parameters(1) + &
                    lcg_parameters(2), mask32)

        ! Write back lower 32 bits
        seed = int(tmp, c_int)

        ! Logical right shift by 8 bits, mask to 24 bits, scale
        rn = real( iand(ishft(tmp, -8), int(Z'00FFFFFF', c_int64_t)), c_float ) &
                / max_int

    END SUBROUTINE lcg


END MODULE particle_rng_mod