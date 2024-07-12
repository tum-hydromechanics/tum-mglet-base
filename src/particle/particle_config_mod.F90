MODULE particle_config_mod

USE core_mod

IMPLICIT NONE

LOGICAL :: dsim_particles

CHARACTER(len = 7) :: particle_terminal

LOGICAL :: dread_particles
LOGICAL :: dinterp_particles
LOGICAL :: dwrite_particles

INTEGER(intk) :: init_npart
INTEGER(intk) :: plist_len
INTEGER(intk) :: psnapshot_step

REAL(realk) :: D

CONTAINS

    SUBROUTINE init_particle_config()

        TYPE(config_t) :: pconf

        IF (.NOT. fort7%exists("/particles")) THEN
            IF (myid == 0) THEN
                WRITE(*, '("NO PARTICLE CONFIGURATION DATA. NO PARTICLE SIMULATION.")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF

        dsim_particles = .TRUE.

        CALL fort7%get(pconf, "/particles")


        CALL pconf%get_value("/terminal", particle_terminal, "normal")

        IF (TRIM(particle_terminal) == "none" .OR. TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN

            CONTINUE

        ELSE

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: Unknown Terminal Output Specifier for Particle Modules. Using normal instead"
                CASE ("verbose")
                    WRITE(*, *) "WARNING: Unknown Terminal Output Specifier for Particle Modules. Using normal instead"
            END SELECT

        END IF


        CALL pconf%get_value("/read", dread_particles, .FALSE.)

        CALL pconf%get_value("/interp", dinterp_particles, .FALSE.)

        CALL pconf%get_value("/write", dwrite_particles, .TRUE.)


        CALL pconf%get_value("/list_len", plist_len, 1000_intk)

        IF (plist_len <= 0_intk) THEN

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: Particle List Length must be positve. Using default Number ", plist_len, "instead."
                CASE ("verbose")
                    WRITE(*, *) "WARNING: Particle List Length must be positve. Using default Number ", plist_len, "instead."
            END SELECT

            plist_len = 1000_intk

        END IF


        CALL pconf%get_value("/init_np", init_npart, 100_intk)

        IF (init_npart <= 0_intk) THEN

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: Initial Number of Particles must be positve. Using default Number ", init_npart, "instead."
                CASE ("verbose")
                    WRITE(*, *) "WARNING: Initial Number of Particles must be positve. Using default Number ", init_npart, "instead."
            END SELECT

            init_npart = 100_intk

        ELSEIF (init_npart > plist_len) THEN

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: Initial Number of Particles must be smaller than the maximum List Length. Using default Number ", init_npart, "instead."
                CASE ("verbose")
                    WRITE(*, *) "WARNING: Initial Number of Particles must be smaller than the maximum List Length. Using default Number ", init_npart, "instead."
            END SELECT

            init_npart = 100_intk

        END IF


        CALL pconf%get_value("/D", D, 0.0_realk)

        IF (D < 0.0_realk) THEN

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: Diffusion Constant must be positve. Using D = 0 instead (pure Advection)."
                CASE ("verbose")
                    WRITE(*, *) "WARNING: Diffusion Constant must be positve. Using D = 0 instead (pure Advection)."
            END SELECT

            D = 0.0_realk

        END IF

        SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, *) 'PARTICLE CONFIGURATION: '
                    WRITE(*, *) 'Particle related Terminal Output: ', particle_terminal
                    WRITE(*, *) 'Reading of ParticlesDict.txt: ', dread_particles
                    WRITE(*, *) 'Interpolation of flow field: ', dinterp_particles
                    WRITE(*, *) 'Writing of Particle Snapshots: ', dwrite_particles
                    WRITE(*, *) 'Particcle Snapshot Step: ', psnapshot_step
                    WRITE(*, *) 'Maximum Particle List Length: ', plist_len
                    WRITE(*, *) 'Initial number of particles for automated Generation: ', init_npart
                    WRITE(*, *) 'Diffusion Constant: ', D
        END SELECT


        CALL pconf%get_value("/snapshot_step", psnapshot_step, 10_intk)

        IF (psnapshot_step < 1_intk) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "Invalid snapshot step. Should be integer greater or equal to 1. Using snapshot_step = 10 instead."
                CASE ("verbose")
                    WRITE(*, *) "Invalid snapshot step. Should be integer greater or equal to 1. Using snapshot_step = 10 instead."
            END SELECT

            psnapshot_step = 10_intk

        END IF

    END SUBROUTINE init_particle_config

END MODULE particle_config_mod