MODULE particle_config_mod

USE core_mod

IMPLICIT NONE

LOGICAL :: dsim_particles
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
                WRITE(*, '("NO PARTICLE CONFIGURATION DATA")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF

        dsim_particles = .TRUE.

        CALL fort7%get(pconf, "/particles")

        CALL pconf%get_value("/read", dread_particles, .FALSE.)

        CALL pconf%get_value("/interp", dinterp_particles, .FALSE.)

        CALL pconf%get_value("/write", dwrite_particles, .TRUE.)

        CALL pconf%get_value("/init_np", init_npart, 100_intk)

        IF (init_npart <= 0_intk) THEN
            WRITE(*, *) "Initial Number of Particles must be positve. Using default Number instead."
            init_npart = 100_intk
        END IF

        CALL pconf%get_value("/list_len", plist_len, 1000_intk)

        IF (plist_len <= 0_intk) THEN
            WRITE(*, *) "Particle List Length must be positve. Using default Number instead."
            plist_len = 1000_intk
        END IF

        CALL pconf%get_value("/D", D, 0.0_realk)

        IF (D <= 0.0_realk) THEN
            WRITE(*, *) "Diffusion Constant must be positve. Using D = 0 instead (pure Advection)."
            D = 0.0_realk
        ELSE
            WRITE(*, *) 'Diffusion Constant: ', D
        END IF

        CALL pconf%get_value("/snapshot_step", psnapshot_step, 10_intk)

        IF (psnapshot_step < 1_intk) THEN
            WRITE(*, *) "Invalid snapshot step. Should be integer greater or equal to 1. Using snapshot_step = 10 instead."
            psnapshot_step = 10_intk
        END IF

    END SUBROUTINE init_particle_config

END MODULE particle_config_mod