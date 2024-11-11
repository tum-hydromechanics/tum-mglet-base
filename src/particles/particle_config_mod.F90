MODULE particle_config_mod

! This module is responsible for:
! Reading of particle parameters from parameters.json file.

USE comms_mod, ONLY: myid, numprocs
USE err_mod
USE precision_mod
USE timer_mod

IMPLICIT NONE

! ON/OFF SWITCH FOR PARTICLE SIMULATION
LOGICAL :: dsim_particles = .FALSE.

! INPUT
LOGICAL :: dread_particles
LOGICAL :: dread_obstacles

! LIST SPECIFICATION
INTEGER(intk) :: init_npart
INTEGER(intk) :: plist_len
LOGICAL :: list_limit = .FALSE.

! SIMULATION SPECIFICATION
LOGICAL :: dinterp_particles
CHARACTER(len=16) :: prkmethod
REAL(realk) :: D(3) = 0.0_realk

! OUTPUT
CHARACTER(len = 7) :: particle_terminal

INTEGER(intk) :: nsamples
CHARACTER(len = 1) :: slice_dir = "N" ! X, Y, Z or N for None
INTEGER(intk) :: nslice_levels
INTEGER(intk), ALLOCATABLE :: nslices(:)
REAL(realk), ALLOCATABLE :: slice_levels(:)

LOGICAL :: dwrite_particles
INTEGER(intk) :: psnapshot_step


CONTAINS    !===================================

    SUBROUTINE init_particle_config()

        USE config_mod
        USE fort7_mod

        INTEGER(intk) :: i, j, mtstep_temp
        TYPE(config_t) :: pconf
        TYPE(config_t) :: timeconf_temp

        !- - - - - - - - - - - - - - - - - -

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

        CALL pconf%get_value("/prkmethod", prkmethod, "euler")

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/write", dwrite_particles, .TRUE.)

        !- - - - - - - - - - - - - - - - - -

        CALL pconf%get_value("/list_len", plist_len, -1)

        list_limit = .TRUE.

        IF (plist_len <= 0_intk) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: Maximum Particle List Length must be a positve Integer. Using automatic List Length instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: Maximum Particle List Length must be a positve Integer. Using automatic List Length instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            plist_len = -1
            list_limit = .FALSE.

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

        CALL fort7%get(timeconf_temp, "/time")
        CALL timeconf_temp%get_value("/mtstep", mtstep_temp)

        CALL pconf%get_value("/nsamples", nsamples, mtstep_temp)

        IF (nsamples < mtstep_temp) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: The number of samples must be an integer greater or equal to mtstep. Using nsamples = mtstep instead."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: The number of samples must be an integer greater or equal to mtstep. Using nsamples = mtstep instead."
                        WRITE(*, '()')
                END SELECT
            END IF

            nsamples = mtstep_temp

        END IF

        !- - - - - - - - - - - - - - - - - -

        IF (fort7%exists("/particles/slice_dir")) THEN

            CALL pconf%get_value("/slice_dir", slice_dir, "N")

            IF (slice_dir /= "X" .AND. slice_dir /= "Y" .AND. slice_dir /= "Z" .AND. slice_dir /= "N") THEN

                IF (myid == 0) THEN
                    SELECT CASE (TRIM(particle_terminal))
                        CASE ("none")
                            CONTINUE
                        CASE ("normal")
                            WRITE(*, *) "WARNING: Slice direction is not valid. Using slice_dir = N instead (no slice statistics)."
                            WRITE(*, '()')
                        CASE ("verbose")
                            WRITE(*, *) "WARNING: Slice direction is not valid. Using slice_dir = N instead (no slice statistics)."
                            WRITE(*, '()')
                    END SELECT
                END IF

                slice_dir = "N"

            END IF

            !- - - - - - - - - - - - - - - - - -

            IF (slice_dir /= "N") THEN

                CALL pconf%get_value("/nslice_levels", nslice_levels, 1)

                IF (nslice_levels < 1) THEN

                    IF (myid == 0) THEN
                        SELECT CASE (TRIM(particle_terminal))
                            CASE ("none")
                                CONTINUE
                            CASE ("normal")
                                WRITE(*, *) "WARNING: The number of slice levels must be an integer greater than 0. Using nslice_levels = 1 instead."
                                WRITE(*, '()')
                            CASE ("verbose")
                                WRITE(*, *) "WARNING: The number of slice levels must be an integer greater than 0. Using nslice_levels = 1 instead."
                                WRITE(*, '()')
                        END SELECT
                    END IF

                    nslice_levels = 1

                END IF

                !- - - - - - - - - - - - - - - - - -

                ALLOCATE(nslices(nslice_levels))

                CALL pconf%get_array("/nslices", nslices)

                DO i = 1, nslice_levels

                    IF (nslices(i) < 1) THEN

                        IF (myid == 0) THEN
                            SELECT CASE (TRIM(particle_terminal))
                                CASE ("none")
                                    CONTINUE
                                CASE ("normal")
                                    WRITE(*, *) "WARNING: The number of slices must be an integer greater than 0. Using nslices(i) = 1 for all i instead."
                                    WRITE(*, '()')
                                CASE ("verbose")
                                    WRITE(*, *) "WARNING: The number of slices must be an integer greater than 0. Using nslices(i) = 1 for all i instead."
                                    WRITE(*, '()')
                            END SELECT
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
                                SELECT CASE (TRIM(particle_terminal))
                                    CASE ("none")
                                        CONTINUE
                                    CASE ("normal")
                                        WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                        WRITE(*, '()')
                                    CASE ("verbose")
                                        WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                        WRITE(*, '()')
                                END SELECT
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
                                SELECT CASE (TRIM(particle_terminal))
                                    CASE ("none")
                                        CONTINUE
                                    CASE ("normal")
                                        WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                        WRITE(*, '()')
                                    CASE ("verbose")
                                        WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                        WRITE(*, '()')
                                END SELECT
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
                        SELECT CASE (TRIM(particle_terminal))
                            CASE ("none")
                                CONTINUE
                            CASE ("normal")
                                WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                WRITE(*, '()')
                            CASE ("verbose")
                                WRITE(*, *) "WARNING: Slice levels are invalid. Using automatically generated levels instead."
                                WRITE(*, '()')
                        END SELECT
                    END IF

                    slice_levels(1) = 1/SUM(nslices)
                    DO j = 2, nslice_levels
                        slice_levels(j) = slice_levels(j-1) + 1/SUM(nslices)
                    END DO

                END IF

            END IF

        END IF

        !- - - - - - - - - - - - - - - - - -


        !- - - - - - - - - - - - - - - - - -

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, '("PARTICLE CONFIGURATION:")')
                        WRITE(*, '("Reading ParticlesDict:        ", L12)') dread_particles
                        WRITE(*, '("Reading ObstaclesDict:        ", L12)') dread_obstacles
                        WRITE(*, '("Max Particle List Length:     ", I12)') plist_len
                        WRITE(*, '("Initial number of Particles:  ", I12)') init_npart
                        WRITE(*, '("Field Interpolation:          ", L12)') dinterp_particles
                        WRITE(*, '("Runge Kutta Method:           ", A12)') TRIM(prkmethod)
                        WRITE(*, '("Diffusion Constant Dx:        ", E12.3)') D(1)
                        WRITE(*, '("Diffusion Constant Dy:        ", E12.3)') D(2)
                        WRITE(*, '("Diffusion Constant Dz:        ", E12.3)') D(3)
                        WRITE(*, '("Terminal Output:              ", A12)') TRIM(particle_terminal)
                        WRITE(*, '("Number of Stat. Samples:      ", I12)') nsamples
                        WRITE(*, '("Slice Direction:              ", A12)') slice_dir
                        IF (slice_dir == "X" .OR. slice_dir == "Y" .OR. slice_dir == "Z") THEN
                            WRITE(*, '("Number of Slice Levels:       ", I12)') nslice_levels
                            DO i = 1, SIZE(nslices)
                                WRITE(*, '("Number of Slices on Level ", I2 ,": ", I12)') i, nslices(i)
                                WRITE(*, '("Width of Level ", I2 ,":            ", E12.3)') i, slice_levels(i)
                            END DO
                        END IF
                        WRITE(*, '("Writing Particle Snapshots:   ", L12)') dwrite_particles
                        WRITE(*, '("Particle Snapshot Step:       ", I12)') psnapshot_step
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, '("PARTICLE CONFIGURATION:")')
                        WRITE(*, '("Reading ParticlesDict:        ", L12)') dread_particles
                        WRITE(*, '("Reading ObstaclesDict:        ", L12)') dread_obstacles
                        WRITE(*, '("Max Particle List Length:     ", I12)') plist_len
                        WRITE(*, '("Initial number of Particles:  ", I12)') init_npart
                        WRITE(*, '("Field Interpolation:          ", L12)') dinterp_particles
                        WRITE(*, '("Runge Kutta Method:           ", A12)') TRIM(prkmethod)
                        WRITE(*, '("Diffusion Constant Dx:        ", E12.3)') D(1)
                        WRITE(*, '("Diffusion Constant Dy:        ", E12.3)') D(2)
                        WRITE(*, '("Diffusion Constant Dz:        ", E12.3)') D(3)
                        WRITE(*, '("Terminal Output:              ", A12)') TRIM(particle_terminal)
                        WRITE(*, '("Number of Stat. Samples:      ", I12)') nsamples
                        WRITE(*, '("Slice Direction:              ", A12)') slice_dir
                        IF (slice_dir == "X" .OR. slice_dir == "Y" .OR. slice_dir == "Z") THEN
                            WRITE(*, '("Number of Slice Levels:       ", I12)') nslice_levels
                            DO i = 1, SIZE(nslices)
                                WRITE(*, '("Number of Slices on Level ", I2 ,": ", I12)') i, nslices(i)
                                WRITE(*, '("Width of Level ", I2 ,":            ", E12.3)') i, slice_levels(i)
                            END DO
                        END IF
                        WRITE(*, '("Writing Particle Snapshots:   ", L12)') dwrite_particles
                        WRITE(*, '("Particle Snapshot Step:       ", I12)') psnapshot_step
                        WRITE(*, '()')
            END SELECT
        END IF

    END SUBROUTINE init_particle_config

END MODULE particle_config_mod