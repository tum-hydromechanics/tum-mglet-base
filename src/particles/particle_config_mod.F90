MODULE particle_config_mod

    ! This module is responsible for:
    ! Reading and storage of particle simulation parameters from parameters.json file.

    USE MPI_f08
    USE comms_mod
    USE err_mod
    USE precision_mod
    USE timer_mod

    IMPLICIT NONE

    ! ON/OFF SWITCH FOR PARTICLE SIMULATION
    LOGICAL :: dsim_particles = .FALSE.

    ! INPUT
    LOGICAL :: dread_particles_dict
    LOGICAL :: dread_obstacles

    ! OUTPUT
    CHARACTER(len = 7) :: particle_terminal
    LOGICAL :: dwrite_npcfield
    LOGICAL :: dwrite_psnapshots
    INTEGER(intk) :: psnapshot_step

    ! RANDOM NUMBER GENERATION
    LOGICAL :: dput_seed
    INTEGER(int32), ALLOCATABLE :: particle_seed(:)

    ! PARTICLE LIST
    INTEGER(intk) :: init_npart ! only relevant if particles are not read
    INTEGER(intk) :: plist_len
    LOGICAL :: list_limit = .FALSE.

    ! ADVECTION
    CHARACTER(len=16) :: prkmethod
    LOGICAL :: dinterp_padvection = .TRUE.

    ! DIFFUSION / RANDOM WALK
    LOGICAL :: ddiffusion
    LOGICAL :: dturb_diff
    LOGICAL :: dinterp_pdiffsuion = .TRUE.
    CHARACTER(len = 16) :: random_walk_mode
    REAL(realk) :: truncation_limit
    REAL(realk) :: D(3) = 0.0_realk

    ! STATISTICS (GRID AND SLICE SAMPLES)
    INTEGER(intk) :: nsamples
    CHARACTER(len = 1) :: slice_dir = "N" ! X, Y, Z or N for None ! MAX ONE DIRECTION !
    INTEGER(intk) :: nslice_levels
    INTEGER(intk), ALLOCATABLE :: nslices(:)
    REAL(realk), ALLOCATABLE :: slice_levels(:)

CONTAINS

    SUBROUTINE init_particle_config()

        USE config_mod
        USE fort7_mod

        TYPE(config_t) :: pconf
        TYPE(config_t) :: timeconf_temp
        INTEGER(intk) :: i, j, mtstep_temp, seed_n, pseed_n, dummy
        INTEGER(int32), ALLOCATABLE ::  seed(:)

        != = = = = = = = = = PARTICLE SIM ON/OFF = = = = = = = = = =

        IF (.NOT. fort7%exists("/particles")) THEN

            IF (myid == 0) THEN
                WRITE(*, '("NO PARTICLE CONFIGURATION DATA AVAILABLE. NO PARTICLE SIMULATION WILL BE CONDUCTED.")')
                WRITE(*, '()')
            END IF

            RETURN

        END IF

        dsim_particles = .TRUE.

        CALL fort7%get(pconf, "/particles")

        != = = = = = = = = = INPUT = = = = = = = = = =

        CALL pconf%get_value("/read", dread_particles_dict, .FALSE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/read_obst", dread_obstacles, .FALSE.)

        != = = = = = = = = = OUTPUT = = = = = = = = = =

        CALL pconf%get_value("/terminal", particle_terminal, "normal")

        IF (TRIM(particle_terminal) == "none" .OR. TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
            CONTINUE
        ELSEIF (myid == 0) THEN
            WRITE(*, *) "WARNING: Unknown Terminal Output Specifier for Particle Modules. Using normal instead"
            WRITE(*, '()')
        END IF

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/dwrite_npc", dwrite_npcfield, .FALSE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/snapshot_step", psnapshot_step, 0)

        dwrite_psnapshots = .TRUE.
        IF (psnapshot_step < 1_intk) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: Snapshot step should be an integer greater or equal to 1. Not writing Particle Snapshots."
                    WRITE(*, '()')
                END IF
            END IF

            psnapshot_step = 0
            dwrite_psnapshots = .FALSE.

        END IF

        != = = = = = = = = = RANDOM SEED = = = = = = = = = =

        ! initialize random number generation
        CALL RANDOM_SEED(size = seed_n)
        ALLOCATE(seed(seed_n))

        dput_seed = .TRUE.
        IF (fort7%exists("/particles/particle_seed")) THEN
            CALL pconf%get_size("/particle_seed", pseed_n)
            IF (pseed_n == (seed_n * numprocs)) THEN
                ALLOCATE(particle_seed(pseed_n))
                CALL pconf%get_array("/particle_seed", particle_seed)
            ELSE
                dput_seed = .FALSE.
            END IF
        ELSE
            dput_seed = .FALSE.
        END IF

        IF (dput_seed) THEN
            DO i = 1, seed_n
                seed(i) = particle_seed(seed_n * myid + i)
            END DO
            CALL RANDOM_SEED(put = seed)
        ELSE
            ! TODO: use RANDOM_INIT instead?
            CALL RANDOM_SEED(get = seed)
            !ALLOCATE(particle_seed(seed_n * numprocs))
            !DO i = 1, seed_n
            !    particle_seed(seed_n * myid + i) = seed(i)
            !END DO
        END IF

        != = = = = = = = = = PARTICLE LIST = = = = = = = = = =

        CALL pconf%get_value("/list_len", plist_len, -1)

        list_limit = .TRUE.

        IF (plist_len <= 0_intk) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: Maximum Particle List Length must be a positve Integer. Using automatic List Length instead."
                    WRITE(*, '()')
                END IF
            END IF

            plist_len = -1
            list_limit = .FALSE.

        END IF

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/init_np", init_npart, 100_intk)

        IF (init_npart <= 0_intk) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: Initial Number of Particles must be positve. Using default Number ", init_npart, "instead."
                    WRITE(*, '()')
                END IF
            END IF

            init_npart = 100_intk

        END IF

        != = = = = = = = = = ADVECTION = = = = = = = = = =

        CALL pconf%get_value("/interp", dinterp_padvection, .TRUE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/prkmethod", prkmethod, "euler")

        != = = = = = = = = = DIFFUSION = = = = = = = = = =

        CALL pconf%get_value("/dturb_diff", dturb_diff, .FALSE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/interp", dinterp_padvection, .TRUE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/random_walk_mode", random_walk_mode, "uniform")

        CALL pconf%get_value("/truncation_limit", truncation_limit, 2.0)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_array("/D", D)

        IF (D(1) < 0.0_realk) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: Diffusion Constant Dx must be positve. Using Dx = 0 instead."
                    WRITE(*, '()')
                END IF
            END IF

            D(1) = 0.0_realk

        END IF

        IF (D(2) < 0.0_realk) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: Diffusion Constant Dy must be positve. Using Dy = 0 instead."
                    WRITE(*, '()')
                END IF
            END IF

            D(2) = 0.0_realk

        END IF

        IF (D(3) < 0.0_realk) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: Diffusion Constant Dz must be positve. Using Dz = 0 instead."
                    WRITE(*, '()')
                END IF
            END IF

            D(3) = 0.0_realk

        END IF

        !- - - - - - - - - - - - - - - - - -

        IF (D(1) == 0.0 .AND. D(2) == 0.0 .AND. D(3) == 0.0 .AND. .NOT. dturb_diff) THEN
            ddiffusion = .FALSE.
            dput_seed = .FALSE.
        END IF

        IF (.NOT. dturb_diff) THEN
            dinterp_pdiffsuion = .FALSE.
        END IF

        != = = = = = = = = = STATISTICS = = = = = = = = = =

        CALL fort7%get(timeconf_temp, "/time")
        CALL timeconf_temp%get_value("/mtstep", mtstep_temp)

        CALL pconf%get_value("/nsamples", nsamples, mtstep_temp)

        IF (nsamples < mtstep_temp) THEN

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: The number of samples must be an integer greater or equal to mtstep. Using nsamples = mtstep instead."
                    WRITE(*, '()')
                END IF
            END IF

            nsamples = mtstep_temp

        END IF

        !- - - - - - - - - - - - - - - - - -

        IF (fort7%exists("/particles/slice_dir")) THEN

            CALL pconf%get_value("/slice_dir", slice_dir, "N")

            IF (slice_dir /= "X" .AND. slice_dir /= "Y" .AND. slice_dir /= "Z" .AND. slice_dir /= "N") THEN

                IF (myid == 0) THEN
                    IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                        WRITE(*, *) "WARNING: Slice direction is not valid. Using slice_dir = N instead (no slice statistics)."
                        WRITE(*, '()')
                    END IF
                END IF

                slice_dir = "N"

            END IF

            !- - - - - - - - - - - - - - - - - -

            IF (slice_dir /= "N") THEN

                CALL pconf%get_value("/nslice_levels", nslice_levels, 1)

                IF (nslice_levels < 1) THEN

                    IF (myid == 0) THEN
                        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                            WRITE(*, *) "WARNING: The number of slice levels must be an integer greater than 0. Using nslice_levels = 1 instead."
                            WRITE(*, '()')
                        END IF
                    END IF

                    nslice_levels = 1

                END IF

                !- - - - - - - - - - - - - - - - - -

                ALLOCATE(nslices(nslice_levels))

                CALL pconf%get_array("/nslices", nslices)

                DO i = 1, nslice_levels

                    IF (nslices(i) < 1) THEN

                        IF (myid == 0) THEN
                            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                                WRITE(*, *) "WARNING: The number of slices must be an integer greater than 0. Using nslices(i) = 1 for all i instead."
                                WRITE(*, '()')
                            END IF
                        END IF

                        DO j = 1, nslice_levels
                            nslices(i) = 1
                        END DO

                        EXIT

                    END IF

                END DO

                !- - - - - - - - - - - - - - - - - -

                ALLOCATE(slice_levels(nslice_levels))

                CALL pconf%get_array("/slice_levels", slice_levels)

                DO i = 1, nslice_levels

                    IF (i == 1) THEN

                        IF (slice_levels(i) < 0.0 .OR. slice_levels(i) > 1.0) THEN

                            IF (myid == 0) THEN
                                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                                    WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                    WRITE(*, '()')
                                END IF
                            END IF

                            slice_levels(1) = 1/SUM(nslices)
                            DO j = 2, nslice_levels
                                slice_levels(j) = slice_levels(j-1) + 1/SUM(nslices)
                            END DO

                            EXIT

                        END IF

                    ELSEIF (i > 1) THEN

                        IF (slice_levels(i) < 0.0 .OR. slice_levels(i) > 1.0) THEN

                            IF (myid == 0) THEN
                                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                                    WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                    WRITE(*, '()')
                                END IF
                            END IF

                            slice_levels(1) = 1/SUM(nslices)
                            DO j = 2, nslice_levels
                                slice_levels(j) = slice_levels(j-1) + 1/SUM(nslices)
                            END DO

                            EXIT

                        END IF

                    END IF

                END DO


                IF (SUM(slice_levels) /= 1.0) THEN

                    IF (myid == 0) THEN
                        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                                WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                WRITE(*, '()')
                        END IF
                    END IF

                    slice_levels(1) = 1/SUM(nslices)
                    DO j = 2, nslice_levels
                        slice_levels(j) = slice_levels(j-1) + 1/SUM(nslices)
                    END DO

                END IF

            END IF

        END IF

        != = = = = = = = = = PRINTOUT = = = = = = = = = =

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                IF (dsim_particles) THEN
                    WRITE(*, '("PARTICLE CONFIGURATION:")')
                    WRITE(*, '("        Particle Simulation:              ", L12)') dsim_particles
                    ! READING
                    WRITE(*, '("    Input:")')
                    WRITE(*, '("        Reading ParticlesDict:            ", L12)') dread_particles_dict
                    WRITE(*, '("        Reading ObstaclesDict:            ", L12)') dread_obstacles
                    ! OUTPUT
                    WRITE(*, '("    Output:")')
                    WRITE(*, '("        Terminal Output:                  ", A12)') TRIM(particle_terminal)
                    WRITE(*, '("        Writing NPC Field:                ", L12)') dwrite_npcfield
                    WRITE(*, '("        Writing Particle Snapshots:       ", L12)') dwrite_psnapshots
                    IF (dwrite_psnapshots) THEN
                    WRITE(*, '("        Particle Snapshot Step:           ", I12)') psnapshot_step
                    END IF
                    ! PARTICLE LIST
                    WRITE(*, '("    Particle List:")')
                    WRITE(*, '("        Limiting initial List Length:     ", L12)') list_limit
                    IF (list_limit) THEN
                    WRITE(*, '("        Max Particle List Length:         ", I12)') plist_len
                    END IF
                    IF (.NOT. dread_particles_dict) THEN
                    WRITE(*, '("        Initial number of Particles:      ", I12)') init_npart
                    END IF
                    ! ADVECTION
                    WRITE(*, '("    Advection:")')
                    WRITE(*, '("        Interpolating Advection:          ", L12)') dinterp_padvection
                    WRITE(*, '("        Runge Kutta Method:               ", A12)') TRIM(prkmethod)
                    ! DIFFUSION
                    WRITE(*, '("    Diffusion:")')
                    WRITE(*, '("        Simulating Diffusion:             ", L12)') ddiffusion
                    IF (ddiffusion) THEN
                    WRITE(*, '("        Simulating turbulent Diffusion:   ", L12)') dturb_diff
                    IF (dturb_diff) THEN
                    WRITE(*, '("        Interpolating Diffusion:          ", L12)') dinterp_pdiffsuion
                    END IF
                    WRITE(*, '("        Random Walk Mode:                 ", A12)') TRIM(random_walk_mode)
                    WRITE(*, '("        Symmetric Truncation Limit:       ", F12.3)') truncation_limit
                    WRITE(*, '("        Global Diffusion Const. Dx:       ", E12.3)') D(1)
                    WRITE(*, '("        Global Diffusion Const. Dy:       ", E12.3)') D(2)
                    WRITE(*, '("        Global Diffusion Const. Dz:       ", E12.3)') D(3)
                    END IF
                    ! STATISTICS
                    WRITE(*, '("    Statistics:")')
                    WRITE(*, '("        Number of Stat. Samples:          ", I12)') nsamples
                    WRITE(*, '("        Slice Direction:                  ", A12)') slice_dir
                    IF (slice_dir == "X" .OR. slice_dir == "Y" .OR. slice_dir == "Z") THEN
                    WRITE(*, '("        Number of Slice Levels:           ", I12)') nslice_levels
                    DO i = 1, SIZE(nslices)
                    WRITE(*, '("            Number of Slices on Level ", I2 ,": ", I12)') i, nslices(i)
                    WRITE(*, '("            Width of Level ", I2 ,":          ", E12.3)') i, slice_levels(i)
                    END DO
                    END IF
                    WRITE(*, '()')
                ELSE
                    WRITE(*, '("PARTICLE CONFIGURATION:")')
                    WRITE(*, '("Particle Simulation:          ", L12)') dsim_particles
                    WRITE(*, '()')
                END IF
            END IF
        END IF

        ! RANDOM NUMBER GENERATION
        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
            IF (myid /= 0) THEN
                CALL MPI_Recv(dummy, 1, mglet_mpi_int, myid - 1, 900, &
                MPI_COMM_WORLD, MPI_STATUS_IGNORE)
            END IF

            IF (myid == 0) THEN
                WRITE(*, '("    Random Number Generation:")')
            END IF
            WRITE(*, '("        Particle Seed on Process ", I0, ":")') myid
            WRITE(*, *) seed
            WRITE(*, '()')

            IF (myid /= numprocs - 1) THEN
                CALL MPI_Send(dummy, 1, mglet_mpi_int, myid + 1, 900, &
                MPI_COMM_WORLD)
            END IF
        END IF

        DEALLOCATE(seed)

        CALL MPI_Barrier(MPI_COMM_WORLD)

    END SUBROUTINE init_particle_config

    SUBROUTINE finish_particle_config()

        CALL start_timer(900)
        CALL start_timer(910)

        IF (ALLOCATED(nslices)) DEALLOCATE(nslices)
        IF (ALLOCATED(slice_levels)) DEALLOCATE(slice_levels)
        IF (ALLOCATED(particle_seed)) DEALLOCATE(particle_seed)

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE finish_particle_config

END MODULE particle_config_mod