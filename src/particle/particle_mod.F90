MODULE particle_mod

    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)
    USE particle_timeintegration_mod
    USE particle_snapshot_mod
    USE particle_connect_mod

    IMPLICIT NONE(type, external)

    !LOGICAL :: dsim_particles = .TRUE.

CONTAINS

    SUBROUTINE init_particles()

        USE particle_timeintegration_mod
        USE particle_snapshot_mod

        CALL init_particle_config()

        IF (dsim_particles) THEN

            WRITE(*,*) "PARTICLE SIMULATION IS INITIALIZED ..."

            CALL init_particle_list()

            ! initializing intra-level communication
            CALL init_particle_connect()

        ELSE

            WRITE(*,*) "NO PARTICLE SIMULATION"

        END IF

    END SUBROUTINE init_particles


    SUBROUTINE finish_particles

        WRITE(*,*) "TO BE DONE---"

    END SUBROUTINE finish_particles

END MODULE particle_mod
