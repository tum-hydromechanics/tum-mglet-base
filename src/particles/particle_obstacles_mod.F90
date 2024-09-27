MODULE particle_obstacles_mod

    USE particle_config_mod

    IMPLICIT NONE

    TYPE :: obstacle_t

        REAL(realk) :: x
        REAL(realk) :: y
        REAL(realk) :: z
        REAL(realk) :: radius

    END TYPE obstacle_t

    TYPE(obstacle_t), ALLOCATABLE :: obstacles(:)

CONTAINS    !===================================

    ! CAUTION: each process reads all obstacles
    SUBROUTINE read_obstacles()

        ! local variables
        INTEGER(intk) :: unit = 1, dict_len, iobst

        IF (.NOT. dread_obstacles) THEN
            ALLOCATE(obstacles(0))
            RETURN
        END IF

        INQUIRE(file = 'ObstaclesDict.txt', exist = dread_obstacles)

        IF (.NOT. dread_obstacles) THEN

            ALLOCATE(obstacles(0))

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: No file for reading obstacles detected!"
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) ' '
                        WRITE(*, *) "WARNING: No file for reading obstacles detected!"
                        WRITE(*, '()')
                END SELECT
            END IF

            RETURN

        END IF

        OPEN(unit, file = 'ObstaclesDict.txt', status = 'OLD', action = 'READ')

        READ(unit, fmt = *) dict_len

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, '("READING ", I0, " OBSTACLE(S) ON ", I0, " PROCESSES EACH.")') dict_len, numprocs
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, '("READING ", I0, " OBSTACLE(S) ON ", I0, " PROCESSES EACH.")') dict_len, numprocs
                    WRITE(*, '()')
            END SELECT
        END IF

        ALLOCATE(obstacles(dict_len))

        DO iobst = 1, dict_len

            READ(unit, fmt = *) obstacles(iobst)%x, obstacles(iobst)%y, &
             obstacles(iobst)%z, obstacles(iobst)%radius

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*,'("Obstacle read: x/y/z/r = ", 4F12.6)') obstacles(iobst)%x, obstacles(iobst)%y, &
                         obstacles(iobst)%z, obstacles(iobst)%radius
                END SELECT
            END IF

        END DO

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, '()')
                    WRITE(*, '("READING OF OBSTACLES SUCCESSFULLY COMPLETED.")')
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, '()')
                    WRITE(*, '("READING OF OBSTACLES SUCCESSFULLY COMPLETED.")')
                    WRITE(*, '()')
            END SELECT
        END IF

    END SUBROUTINE read_obstacles

END MODULE particle_obstacles_mod