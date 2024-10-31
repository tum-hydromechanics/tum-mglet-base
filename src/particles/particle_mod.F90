MODULE particle_mod

    ! This module is responsible for:
    ! Initialization of the particle simulation.
    ! Finishing of the particle simulation.

    USE fields_mod
    USE particle_timeintegration_mod
    USE particle_snapshot_mod
    USE particle_statistics_mod

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
            ! PARTICLES : includes all other timers
            CALL set_timer(900, 'PARTICLES')

            ! PARTICLE_INITIALIZATION:
            ! init_particle_list (particle reading/particle distribution, particle list initialization)
            ! init_particle_boundaries (initialization of particle grid boundaries, reading of obstacles)
            ! init_particle_timeintegration (initialization of velocity backup field)
            CALL set_timer(910, 'PARTICLE_CORE')

            ! PARTICLE_TIMEINTEGRATION:
            ! update_backup_fields (copying of velocity fields at beginning of each timestep)
            ! timeintegrate_particles (particle timeintegration)
            CALL set_timer(920, 'PARTICLE_TIMEINTEGRATION')

                ! VELOCITY BACKUP MANIPULATION:
                ! update_backup_fields (copying of velocity fields at beginning of each timestep)
                ! time_interpolate_field (time interpolation of the the velocity field for particle velocity deduction)
                CALL set_timer(921, 'VELOCITY_BACKUP_MANIPULATION')

                ! PARTICLE DIFFSUION GENERATION:
                ! get_particle_diffusion (random walk computation)
                CALL set_timer(922, 'PARTICLE_DIFFUSION_GENERATION')

                ! PARTICLE ADVECTION GENERATION:
                ! get_particle_uvw (particle velocity deduction)
                ! nearest_particle_uvw (particle velocity deduction)
                ! gobert_particle_uvw (particle velocity deduction with interpolation)
                CALL set_timer(923, 'PARTICLE_ADVECTION_GENERATION')

                ! PARTICLE BOUNDARY INTERACTION:
                ! move_particle (particle motion and interaction with boundaries)
                CALL set_timer(924, 'PARTICLE_BOUNDARY_INTERACTION')

            ! PARTICLE EXCHANGE
            ! exchange particles (particle target grid and proc determination, particle exchange between procs)
            CALL set_timer(930, 'PARTICLE_EXCHANGE')

            ! PARTICLE STATISTICS
            ! init_particle_gridstat (initialization of particle grid statistics)
            ! advance_np_counter (copy the number of particles of each grid from last timestep)
            ! init_particle_gridstat (initialization of particle grid statistics)
            ! register_particle/ deregister_particle (register/ deregister particle on grid)
            CALL set_timer(940, 'PARTICLE_STATISTICS')

            ! PARTICLE SNAPSHOTS
            ! init_psnapshots (initialization of particle snaphots, directory creation)
            ! write_psnapshot (particle snapshot writing)
            ! write_psnapshot_timeinfo (meta info writing)
            CALL set_timer(950, 'PARTICLE_SNAPSHOTS')

            CALL init_particle_list()

            ! determine particle boundaries and their normal vectors
            CALL init_particle_boundaries()

            ! init backup fields for particle timeintegration
            CALL init_particle_timeintegration()

            ! determine particle exchange connections and init particle exchange
            CALL init_particle_exchange()

        ELSE

            IF (myid == 0) THEN
                WRITE(*,*) "NO PARTICLE SIMULATION"
            END IF

        END IF

    END SUBROUTINE init_particles


    SUBROUTINE finish_particles()

        CALL finish_particle_snapshots()

        CALL finish_particle_statistics()

        CALL finish_particle_timeintegration()

        CALL finish_particle_exchange()

        CALL finish_particle_boundaries()

        CALL finish_particle_list()

        IF (myid == 0) THEN
            WRITE(*,*) "PARTICLE SIMULATION COMPLETED."
        END IF

    END SUBROUTINE finish_particles

END MODULE particle_mod
