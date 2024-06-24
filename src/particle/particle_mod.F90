MODULE particle_mod

    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)
    USE particle_timeintegration_mod
    USE particle_snapshot_mod

    IMPLICIT NONE(type, external)

CONTAINS

    SUBROUTINE init_particles()
        USE particle_timeintegration_mod
        USE particle_snapshot_mod

        WRITE(*,*) "HELLO WORLD PARTICLE(s)---"
        CALL init_particle_list()

    END SUBROUTINE init_particles


    SUBROUTINE finish_particles

        WRITE(*,*) "TO BE DONE---"

    END SUBROUTINE finish_particles

END MODULE particle_mod
