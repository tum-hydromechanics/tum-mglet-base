MODULE particle_config_mod

! This module is responsible for:
! Reading of particle parameters from parameters.json file.

USE comms_mod, ONLY: myid, numprocs
USE err_mod
USE precision_mod
USE timer_mod

IMPLICIT NONE

LOGICAL :: dsim_particles = .FALSE.

CHARACTER(len = 7) :: particle_terminal

LOGICAL :: dread_particles
LOGICAL :: dread_obstacles
LOGICAL :: dinterp_particles
LOGICAL :: dwrite_particles

INTEGER(intk) :: init_npart
INTEGER(intk) :: plist_len
INTEGER(intk) :: psnapshot_step

REAL(realk) :: D(3) = 0.0_realk

CHARACTER(len=16) :: prkmethod

CONTAINS    !===================================

    SUBROUTINE init_particle_config()

        USE config_mod
        USE fort7_mod

        !- - - - - - - - - - - - - - - - - -

        TYPE(config_t) :: pconf

        IF (.NOT. fort7%exists("/particles")) THEN

            IF (myid == 0) THEN
                WRITE(*, '("NO PARTICLE CONFIGURATION DATA AVAILABLE. NO PARTICLE SIMULATION WILL BE CONDUCTED.")')
                WRITE(*, '()')
            END IF

            RETURN

        END IF

        dsim_particles = .TRUE.

        CALL fort7%get(pconf, "/particles")

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/terminal", particle_terminal, "normal")

        IF (TRIM(particle_terminal) == "none" .OR. TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN

            CONTINUE

        ELSEIF (myid == 0) THEN

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: Unknown Terminal Output Specifier for Particle Modules. Using normal instead"
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) "WARNING: Unknown Terminal Output Specifier for Particle Modules. Using normal instead"
                    WRITE(*, '()')
            END SELECT

        END IF

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/read", dread_particles, .FALSE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/read_obst", dread_obstacles, .FALSE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/interp", dinterp_particles, .TRUE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/prkmethod", prkmethod, "williamson")

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/write", dwrite_particles, .TRUE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/list_len", plist_len, 1000_intk)


        IF (plist_len <= 0_intk) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: Particle List Length must be positve. Using default Number ", plist_len, "instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: Particle List Length must be positve. Using default Number ", plist_len, "instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            plist_len = 1000_intk

        END IF

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/init_np", init_npart, 100_intk)

        IF (init_npart <= 0_intk) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: Initial Number of Particles must be positve. Using default Number ", init_npart, "instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: Initial Number of Particles must be positve. Using default Number ", init_npart, "instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            init_npart = 100_intk

        ELSEIF (init_npart > plist_len) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: Initial Number of Particles must be smaller than the maximum List Length. Using default Number ", &
                         init_npart, "instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: Initial Number of Particles must be smaller than the maximum List Length. Using default Number ", &
                         init_npart, "instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            init_npart = 100_intk

        END IF

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_array("/D", D)

        IF (D(1) < 0.0_realk) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: Diffusion Constant Dx must be positve. Using Dx = 0 instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: Diffusion Constant Dx must be positve. Using Dx = 0 instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            D(1) = 0.0_realk

        END IF

        IF (D(2) < 0.0_realk) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: Diffusion Constant Dy must be positve. Using Dy = 0 instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: Diffusion Constant Dy must be positve. Using Dy = 0 instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            D(2) = 0.0_realk

        END IF

        IF (D(3) < 0.0_realk) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: Diffusion Constant Dz must be positve. Using Dz = 0 instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: Diffusion Constant Dz must be positve. Using Dz = 0 instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            D(3) = 0.0_realk

        END IF

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/snapshot_step", psnapshot_step, 10_intk)

        IF (psnapshot_step < 1_intk) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: Snapshot step should be an integer greater or equal to 1. Using snapshot_step = 10 instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: Snapshot step should be an integer greater or equal to 1. Using snapshot_step = 10 instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            psnapshot_step = 10_intk

        END IF

        !- - - - - - - - - - - - - - - - - -

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, '("PARTICLE CONFIGURATION:")')
                        WRITE(*, '("Terminal Output:              ", A12)') particle_terminal
                        WRITE(*, '("Reading ParticlesDict:        ", L12)') dread_particles
                        WRITE(*, '("Reading ObstaclesDict:        ", L12)') dread_obstacles
                        WRITE(*, '("Field Interpolation:          ", L12)') dinterp_particles
                        WRITE(*, '("Writing Particle Snapshots:   ", L12)') dwrite_particles
                        WRITE(*, '("Particle Snapshot Step:       ", I12)') psnapshot_step
                        WRITE(*, '("Max Particle List Length:     ", I12)') plist_len
                        WRITE(*, '("Initial number of Particles:  ", I12)') init_npart
                        WRITE(*, '("Diffusion Constant Dx:        ", E12.3)') D(1)
                        WRITE(*, '("Diffusion Constant Dy:        ", E12.3)') D(2)
                        WRITE(*, '("Diffusion Constant Dz:        ", E12.3)') D(3)
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, '("PARTICLE CONFIGURATION:")')
                        WRITE(*, '("Terminal Output:              ", A12)') particle_terminal
                        WRITE(*, '("Reading ParticlesDict:        ", L12)') dread_particles
                        WRITE(*, '("Reading ObstaclesDict:        ", L12)') dread_obstacles
                        WRITE(*, '("Field Interpolation:          ", L12)') dinterp_particles
                        WRITE(*, '("Writing Particle Snapshots:   ", L12)') dwrite_particles
                        WRITE(*, '("Particle Snapshot Step:       ", I12)') psnapshot_step
                        WRITE(*, '("Max Particle List Length:     ", I12)') plist_len
                        WRITE(*, '("Initial number of Particles:  ", I12)') init_npart
                        WRITE(*, '("Diffusion Constant Dx:        ", E12.3)') D(1)
                        WRITE(*, '("Diffusion Constant Dy:        ", E12.3)') D(2)
                        WRITE(*, '("Diffusion Constant Dz:        ", E12.3)') D(3)
                        WRITE(*, '()')
            END SELECT
        END IF

    END SUBROUTINE init_particle_config

END MODULE particle_config_mod