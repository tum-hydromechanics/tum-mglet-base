MODULE particle_mod

    ! This module is responsible for:
    ! Initialization of the particle simulation.
    ! Finishing of the particle simulation.

    USE fields_mod

    USE particle_runtimestat_mod
    USE particle_timeintegration_mod
    USE particle_statistics_mod
    USE particle_snapshot_mod
    USE particle_fields_mod

    IMPLICIT NONE(type, external)

CONTAINS

    SUBROUTINE init_particles()

        CALL init_particle_config()

        IF (dsim_particles) THEN

            IF (myid == 0) THEN
                WRITE(*, '("PARTICLE SIMULATION STARTED.")')
                WRITE(*, '()')
            END IF

            ! --- TIMERS ---
            ! TODO: restructure
            ! 900: particles
            !  910: initialization
            !   911: init_particle_boundaries
            !   ...
            !  920: timeintegration
            !   ...
            !  930: timeloop extras
            !   931: fields
            !   932: statistics
            !   933: snapshots
            !  940: finishing
            !   ...

            ! PARTICLES : includes all other timers
            CALL set_timer(900, 'PARTICLE_SIMULATION')

            ! PARTICLE SIMULATION CORE INITIALIZATION:
            ! init_particle_boundaries()
            ! init_particle_list()
            ! init_particle_diffusion()
            ! init_particle_timeintegration()
            ! init_particle_exchange()
            ! init_particle_field()
            CALL set_timer(910, 'PSIM_CORE_INIT')

            ! PARTICLE TIMEINTEGRATION:
            CALL set_timer(920, 'PSIM_TIMEINTEGRATION')
                CALL set_timer(921, 'ADV_VELOCITY')
                CALL set_timer(922, 'ADV_MOTION')
                CALL set_timer(923, 'DIF_COEFF')
                CALL set_timer(924, 'DIF_RW_GENERATION')
                CALL set_timer(925, 'DIF_MOTION')

            ! PARTICLE BOUNDARY INTERACTION:
            ! ... see TIMEINTEGRATION

            ! PARTICLE EXCHANGE
            CALL set_timer(940, 'PSIM_EXCHANGE')

            ! PARTICLE STATISTICS
            ! everything relatet to statistics (incl. init/finish_particle_statistics)
            CALL set_timer(950, 'PSIM_STATISTICS')

            ! PARTICLE SNAPSHOTS
            ! everything related to snapshots (incl. init/finish_particle_snapshots)
            CALL set_timer(960, 'PSIM_SNAPSHOTS')

            ! PARTICLE SIMULATION CORE FINISHING:
            ! finish_particle_boundaries()
            ! finish_particle_list()
            ! finish_particle_exchange()
            ! finish_particle_timeintegration()
            ! finish_particle_config()
            CALL set_timer(990, 'PSIM_CORE_FINISH')

            ! determine particle boundaries and their normal vectors
            CALL init_particle_boundaries()

            ! read or generate particles and init particle list
            CALL init_particle_list()

            ! generate diffusion field
            !CALL init_particle_diffusion() EDIT: MUST BE DONE IN INIT TIMELOOP AFTER INIT_STATISTICS
            ! init backup fields for particle timeintegration
            !CALL init_particle_timeintegration() EDIT: MUST BE DONE IN INIT TIMELOOP AFTER INIT_STATISTICS

            ! determine particle exchange connections and init particle exchange
            CALL init_particle_exchange()

            ! set particle concentration field
            CALL init_particle_field()

            !  DIFFUSION, TIMEINTEGRATION, STATISTICS AND SNAPSHOT INITIALIZATION IN TIMELOOP

        ELSE

            IF (myid == 0) THEN
                WRITE(*,*) "NO PARTICLE SIMULATION"
            END IF

        END IF

    END SUBROUTINE init_particles

    SUBROUTINE finish_particles()

        IF (dsim_particles) THEN

            CALL finish_particle_snapshots()

            CALL finish_particle_statistics()

            CALL finish_particle_timeintegration()

            CALL finish_particle_exchange()

            CALL finish_particle_list()

            CALL finish_particle_boundaries()

            IF (myid == 0) THEN
                WRITE(*,*) "PARTICLE SIMULATION FINISHED SUCCESSFULLY."
            END IF

        END IF

        CALL finish_particle_config()

    END SUBROUTINE finish_particles

END MODULE particle_mod
