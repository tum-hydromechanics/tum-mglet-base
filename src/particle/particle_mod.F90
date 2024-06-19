MODULE particle_mod

    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)
    USE particle_timeintegration_mod

    IMPLICIT NONE(type, external)

CONTAINS
    SUBROUTINE init_particle()
        USE particle_timeintegration_mod
        USE particle_io_mod

        WRITE(*,*) "HELLO WORLD PARTICLE(s)---"
        CALL init_particles()

    END SUBROUTINE init_particle


    SUBROUTINE finish_particle

        WRITE(*,*) "TO BE DONE---"

    END SUBROUTINE finish_particle

END MODULE particle_mod
