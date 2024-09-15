MODULE particle_mod

    ! This module is responsible for:
    ! Initialization of the particle simulation.
    ! Finishing of the particle simulation.


    USE particle_timeintegration_mod
    USE particle_snapshot_mod

    IMPLICIT NONE(type, external)

CONTAINS

    SUBROUTINE init_particles()

        CALL init_particle_config()

        IF (dsim_particles) THEN

            IF (myid == 0) THEN
                WRITE(*,*) "PARTICLE SIMULATION STARTED."
            END IF

            CALL set_timer(900, 'PARTICLES')
            CALL set_timer(910, 'PARTICLE_TIMEINTEGRATION')
            CALL set_timer(920, 'PARTICLE_SNAPSHOTS')

            CALL start_timer(900)

            CALL init_particle_list()

            CALL init_particle_boundaries()

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

        IF (myid == 0) THEN
            WRITE(*,*) "PARTICLE SIMULATION COMPLETED."
        END IF

    END SUBROUTINE finish_particles

END MODULE particle_mod
