MODULE particle_mod

    ! This module is responsible for:
    ! Initialization of the particle simulation.
    ! Finishing of the particle simulation.

    USE fields_mod

    USE particle_list_mod
    USE particle_runtimestat_mod
    USE particle_timeintegration_mod
    USE particle_statistics_mod
    USE particle_snapshot_mod
    USE particle_io_mod
    USE particle_loadbalance_mod
    
    IMPLICIT NONE

CONTAINS

    SUBROUTINE init_particles()

        CALL init_particle_config()

        IF (dsim_particles) THEN

            IF (myid == 0) THEN
                WRITE(*, '("PARTICLE SIMULATION STARTED.")')
                WRITE(*, '()')
            END IF

            ! --- TIMERS ---

            ! PARTICLES: includes all other timers
            CALL set_timer(900, 'PARTICLE_SIMULATION')

            ! PARTICLE SIMULATION CORE INITIALIZATION:
            ! init_particle_boundaries()
            ! init_particle_list(); read_particles_h5()
            ! init_particle_diffusion()
            ! init_particle_timeintegration()
            ! init_particle_exchange()
            CALL set_timer(910, 'PSIM_CORE_INIT')

            ! PARTICLE TIMEINTEGRATION:
            CALL set_timer(920, 'PSIM_TIMEINTEGRATION')
                CALL set_timer(921, 'ADV_VELOCITY')
                CALL set_timer(922, 'ADV_MOTION')
                CALL set_timer(924, 'DIF_RN_GENERATION')
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
            ! finish_particle_list(); write_particles_h5()
            ! finish_particle_exchange()
            ! finish_particle_timeintegration()
            ! finish_particle_config()
            CALL set_timer(990, 'PSIM_CORE_FINISH')

#ifdef _MGLET_OPENMP_
            CALL offload_fields()
#endif

            CALL init_particle_utils()

            ! determine particle boundaries and their normal vectors
            CALL init_particle_boundaries()

            ! read or generate particles and init particle list
            CALL init_particle_list()

            IF (dread_particles_h5) THEN
                CALL read_particles_h5("particles.h5")

                CALL MPI_Allreduce(my_particle_list%active_np, global_np, 1, mglet_mpi_int, MPI_SUM, MPI_COMM_WORLD)

                IF (myid == 0) THEN
                    IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                        WRITE(*, '("INITIALIZATION OF ", I0, " PARTICLE(S) SUCCESSFULLY COMPLETED.")') global_np
                        WRITE(*, '()')
                    END IF
                END IF

            END IF
            
            CALL init_particle_loadbalance()

            ! determine particle exchange connections and init particle exchange
            CALL init_particle_exchange()

            ! DIFFUSION, TIMEINTEGRATION, STATISTICS AND SNAPSHOT INITIALIZATION IN TIMELOOP

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

            CALL finish_particle_boundaries()

            CALL finish_particle_utils()

#ifdef _MGLET_OPENMP_
            CALL finish_offload_fields()
#endif

            ! stupid test case for read / write

            IF (dwrite_particles_h5) THEN
                IF (myid == 0) THEN
                    WRITE(*,*) "Writing the particles"
                END IF
                CALL write_particles_h5("particles.h5")
            END IF

            CALL finish_particle_list()

            IF (myid == 0) THEN
                WRITE(*,*) "PARTICLE SIMULATION FINISHED SUCCESSFULLY."
            END IF

        END IF

        CALL finish_particle_config()

    END SUBROUTINE finish_particles

END MODULE particle_mod
