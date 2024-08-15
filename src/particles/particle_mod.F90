MODULE particle_mod

    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)
    USE particle_timeintegration_mod
    USE particle_snapshot_mod
    USE particle_exchange_mod

    IMPLICIT NONE(type, external)

    !LOGICAL :: dsim_particles = .TRUE.

CONTAINS

    SUBROUTINE init_particles()

        USE particle_timeintegration_mod
        USE particle_snapshot_mod

        CALL init_particle_config()

        IF (dsim_particles) THEN

            IF (myid == 0) THEN
                WRITE(*,*) "PARTICLE SIMULATION STARTED"
            END IF

            CALL set_timer(900, 'PARTICLES')
            CALL set_timer(910, 'PARTICLE_TIMEINTEGRATION')
            CALL set_timer(920, 'PARTICLE_SNAPSHOTS')

            CALL start_timer(900)

            CALL init_particle_list()

            ! initializing intra-level communication
            CALL init_particle_exchange()

            CALL stop_timer(900)

        ELSE

            IF (myid == 0) THEN
                WRITE(*,*) "NO PARTICLE SIMULATION"
            END IF

        END IF

    END SUBROUTINE init_particles


    SUBROUTINE finish_particles()

        ! TODO: DEALLOCATIONS ETC. ...

        WRITE(*,*) "PARTICLE SIMULATION FINISHED"

    END SUBROUTINE finish_particles

END MODULE particle_mod
