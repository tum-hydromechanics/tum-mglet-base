MODULE particle_rng_mod
    
    ! https://vectrx.substack.com/p/lcg-xs-fast-gpu-rng#footnote-7-159943017
    
    USE, INTRINSIC :: ISO_FORTRAN_ENV
    USE, INTRINSIC :: ISO_C_BINDING

    USE omp_lib

    USE precision_mod

    IMPLICIT NONE

    INTEGER(c_int) :: particle_base_seed = 9891477

    INTEGER(c_int64_t) :: lcg_parameters(2) 

#if defined __INTEL_COMPILER
    !$omp declare target link(lcg_parameters)
#else
    !$omp declare target(lcg_parameters)
#endif

CONTAINS

    SUBROUTINE init_parallel_lcg()

        !lcg_multiplier
        lcg_parameters(1) = 747796405_c_int64_t
        !lcg_increment
        lcg_parameters(2) = 2891336453_c_int64_t

        !$omp target enter data map(always, to: lcg_parameters)

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
        INTEGER(c_int64_t) :: tmp = 0_c_int64_t

        REAL(c_float) :: max_int = 16777216.0_c_float
        
        ! convert to 64-bit signed (sign extension) and nullification of bits 32(33)-63(64)
        tmp = iand(int(seed, c_int64_t), Z'FFFFFFFF')
        ! compute RNG step in 64-bit
        tmp = tmp * lcg_parameters(1) + lcg_parameters(2)

        ! update seed
        seed = int(tmp, c_int)

        ! logical right shift by 8 bits, mask to 24 bits, scale
        rn = real(iand(ishft(tmp, -8), int(Z'00FFFFFF', c_int64_t)), c_float) / max_int

    END SUBROUTINE lcg


END MODULE particle_rng_mod