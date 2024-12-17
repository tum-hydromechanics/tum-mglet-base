MODULE particle_obstacles_mod

    ! This module is responsible for:
    ! Reading and storage of a sphere pack (obstacles that make up a porous domain)

    USE MPI_f08
    USE comms_mod
    USE grids_mod
    USE field_mod
    USE fields_mod

    USE particle_utils_mod

    IMPLICIT NONE

    ! TODO: make an abstract parent class
    TYPE :: obstacle_t

        INTEGER(intk) :: iobst

        REAL(realk) :: x
        REAL(realk) :: y
        REAL(realk) :: z

        REAL(realk) :: radius

        CONTAINS
            PROCEDURE :: in_zone

    END TYPE obstacle_t

    ! TODO: rename
    TYPE :: obstacle_pointer_t

        ! list that holds pointers to all obstacles that are relevant on this process for the associated grid
        ! (not neccessarily all obstacles on the grid, if the grid is not part of mygrids)
        INTEGER(intk), ALLOCATABLE :: grid_obstacles(:)

    END TYPE obstacle_pointer_t

    ! list of all obstacles that are on any grid or neighbouring grid of this process
    TYPE(obstacle_t), ALLOCATABLE :: my_obstacles(:)

    ! list that holds one obstacle_pointer_t per grid (all grids, not limited to those on this process)
    TYPE(obstacle_pointer_t), ALLOCATABLE :: my_obstacle_pointers(:)

    ! TODO: change this value?
    ! factor to compute the minimum distance between obstacles for which no intermediate obstacle is generated
    ! using EPSILON(realk)
    REAL(realk), PARAMETER :: min_space_factor = 100

    ! ratio of intermediate (filling) obstacles over readius of regular obstacles
    REAL(realk), PARAMETER :: radius_ratio = 0.348

CONTAINS    !===================================

    ! CAUTION: each process reads all obstacles
    SUBROUTINE read_obstacles()

        ! local variables
        TYPE(obstacle_t), ALLOCATABLE :: obstacles_src(:), obstacles_itm(:) ! temporary obtscle lists
        !(src -> input/regular/source obstacles (read), itm -> intermediate/filling obstacles (generated))

        INTEGER(intk) :: unit, dict_len, iobst, igrid, i, j, k, counter, dummy
        INTEGER(intk) :: neighbours(26)
        INTEGER(intk), ALLOCATABLE :: counter_array(:) ! number of obstacles that is relevant on this process per grid

        LOGICAL :: dcycle, dexit
        LOGICAL, ALLOCATABLE :: grid_processed(:) ! array to indicate if a grid has been considered for a certain obstacle already
        LOGICAL, ALLOCATABLE :: is_relevant_src(:), is_relevant_itm(:) ! indicates if an obstacle is relevant on this process

        REAL(realk) :: dist
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        CHARACTER(12) :: dummy_char

        ! TODO: optimize this routine

        IF (.NOT. dread_obstacles) THEN
            ALLOCATE(my_obstacles(0))
            RETURN
        END IF

        INQUIRE(file = 'ObstaclesDict.txt', exist = dread_obstacles)

        IF (.NOT. dread_obstacles) THEN
            ALLOCATE(my_obstacles(0))
            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "WARNING: No file for reading obstacles detected!"
                    WRITE(*, '()')
                END IF
            END IF

            RETURN

        END IF

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("READING OBSTACLES ...")')
                WRITE(*, '()')
            END IF
        END IF

        ALLOCATE(counter_array(ngrid))
        counter_array = 0
        ALLOCATE(my_obstacle_pointers(ngrid))
        ALLOCATE(grid_processed(ngrid))

        OPEN(newunit = unit, file = 'ObstaclesDict.txt', status = 'OLD', action = 'READ')

        READ(unit, fmt = *) dummy_char, dummy_char, dummy_char, dict_len

        IF (dict_len < 1) THEN
            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, '()')
                    WRITE(*, '("ERROR in particle_obstacles_mod: Number of Obstacles given in ObstaclesDict.txt must be an integer larger than 0.")')
                    WRITE(*, '("If no Obstacles should be registered, set read_obst to FALSE.")')
                    WRITE(*, '()')
                END IF
            END IF
            CALL errr(__FILE__, __LINE__)
        END IF

        ALLOCATE(obstacles_src(dict_len))
        ALLOCATE(is_relevant_src(dict_len))
        is_relevant_src = .FALSE.

        DO iobst = 1, dict_len

            READ(unit, fmt = *) dummy_char

            READ(unit, fmt = *) obstacles_src(iobst)%radius, obstacles_src(iobst)%x, obstacles_src(iobst)%y, &
             obstacles_src(iobst)%z

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("Obstacle read: x/y/z/r = ", 4F12.6)') obstacles_src(iobst)%x, obstacles_src(iobst)%y, &
                     obstacles_src(iobst)%z, obstacles_src(iobst)%radius
                END IF
            END IF

            obstacles_src(iobst)%iobst = iobst

        END DO

        CLOSE(unit)

        ! generate intermediate obstacles to fill space between touching obstacles

        ! https://mathworld.wolfram.com/KissingNumber.html

        counter = 0
        DO i = 1, dict_len - 1
            DO j = i + 1, dict_len

                dist = SQRT((obstacles_src(i)%x - obstacles_src(j)%x)**2 + &
                 (obstacles_src(i)%y - obstacles_src(j)%y)**2 + &
                 (obstacles_src(i)%z - obstacles_src(j)%z)**2)

                ! TODO: change minimum space between obstacles?
                IF (dist < (obstacles_src(i)%radius + obstacles_src(j)%radius + min_space_factor * EPSILON(obstacles_src(i)%radius))) THEN
                    counter = counter + 1
                END IF

            END DO
        END DO

        ALLOCATE(obstacles_itm(counter))
        ALLOCATE(is_relevant_itm(counter))
        is_relevant_itm = .FALSE.

        counter = 1
        DO i = 1, dict_len - 1
            DO j = i + 1, dict_len

                dist = SQRT((obstacles_src(i)%x - obstacles_src(j)%x)**2 + &
                 (obstacles_src(i)%y - obstacles_src(j)%y)**2 + &
                 (obstacles_src(i)%z - obstacles_src(j)%z)**2)

                IF (dist < (obstacles_src(i)%radius + obstacles_src(j)%radius + min_space_factor * EPSILON(obstacles_src(i)%radius))) THEN

                    obstacles_itm(counter)%iobst = counter + dict_len
                    obstacles_itm(counter)%x = obstacles_src(i)%x + (obstacles_src(j)%x - obstacles_src(i)%x) * obstacles_src(i)%radius / dist
                    obstacles_itm(counter)%y = obstacles_src(i)%y + (obstacles_src(j)%y - obstacles_src(i)%y) * obstacles_src(i)%radius / dist
                    obstacles_itm(counter)%z = obstacles_src(i)%z + (obstacles_src(j)%z - obstacles_src(i)%z) * obstacles_src(i)%radius / dist
                    obstacles_itm(counter)%radius = radius_ratio * obstacles_src(i)%radius

                    counter = counter + 1

                END IF

            END DO
        END DO

        ! count regular/source obstacles that are relevant for this process and determine for which grids they are relevant
        counter = 0
        DO iobst = 1, dict_len

            grid_processed = .FALSE.

            DO i = 1, nmygrids

                igrid = mygrids(i)

                IF (obstacles_src(iobst)%in_zone(igrid) .AND. .NOT. grid_processed(igrid)) THEN
                    grid_processed(igrid) = .TRUE.
                    counter_array(igrid) = counter_array(igrid) + 1
                    is_relevant_src(iobst) = .TRUE.
                END IF

                CALL get_neighbours(neighbours, igrid)

                DO j = 1, 26
                    IF (obstacles_src(iobst)%in_zone(neighbours(j)) .AND. .NOT. grid_processed(neighbours(j))) THEN
                        grid_processed(neighbours(j)) = .TRUE.
                        counter_array(neighbours(j)) = counter_array(neighbours(j)) + 1
                        is_relevant_src(iobst) = .TRUE.
                    END IF
                END DO

                IF (is_relevant_src(iobst)) THEN
                    counter = counter + 1
                END IF

            END DO

        END DO

        ! count intermediate/filling obstacles that are relevant for this process and determine for which grids they are relevant
        DO iobst = 1, SIZE(obstacles_itm)

            grid_processed = .FALSE.

            DO i = 1, nmygrids

                igrid = mygrids(i)

                IF (obstacles_itm(iobst)%in_zone(igrid) .AND. .NOT. grid_processed(igrid)) THEN
                    grid_processed(igrid) = .TRUE.
                    counter_array(igrid) = counter_array(igrid) + 1
                    is_relevant_itm(iobst) = .TRUE.
                END IF

                CALL get_neighbours(neighbours, igrid)

                DO j = 1, 26
                    IF (obstacles_itm(iobst)%in_zone(neighbours(j)) .AND. .NOT. grid_processed(neighbours(j))) THEN
                        grid_processed(neighbours(j)) = .TRUE.
                        counter_array(neighbours(j)) = counter_array(neighbours(j)) + 1
                        is_relevant_itm(iobst) = .TRUE.
                    END IF
                END DO

                IF (is_relevant_src(iobst)) THEN
                    counter = counter + 1
                END IF

            END DO

        END DO

        DO igrid = 1, ngrid
            ALLOCATE(my_obstacle_pointers(igrid)%grid_obstacles(counter_array(igrid)))
        END DO

        ALLOCATE(my_obstacles(counter))

        counter_array = 1
        counter = 1

        ! copy relevant obstacles for this process into heap memory
        DO iobst = 1, dict_len

            IF (is_relevant_src(iobst)) THEN

                my_obstacles(counter) = obstacles_src(iobst)

                DO igrid = 1, ngrid
                    IF (obstacles_src(iobst)%in_zone(igrid)) THEN
                        my_obstacle_pointers(igrid)%grid_obstacles(counter_array(igrid)) = counter
                        counter_array(igrid) = counter_array(igrid) + 1
                    END IF
                END DO

                counter = counter + 1

            END IF

        END DO

        ! copy relevant intermediate/filling obstacles for this process into heap memory
        DO iobst = 1, SIZE(obstacles_itm)

            IF (is_relevant_itm(iobst)) THEN

                my_obstacles(counter) = obstacles_itm(iobst)

                DO igrid = 1, ngrid
                    IF (obstacles_itm(iobst)%in_zone(igrid)) THEN
                        my_obstacle_pointers(igrid)%grid_obstacles(counter_array(igrid)) = counter
                        counter_array(igrid) = counter_array(igrid) + 1
                    END IF
                END DO

                counter = counter + 1

            END IF

        END DO

        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
            IF (myid /= 0) THEN
                CALL MPI_Recv(dummy, 1, mglet_mpi_int, myid - 1, 900, &
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE)
            END IF

            WRITE(*, '("Obstacels on Process ", I0, ":")') myid
            WRITE(*, '()')

            DO igrid = 1, ngrid

                WRITE(*,'(I0, " Obstacles registered for Grid ", I0, ".")') SIZE(my_obstacle_pointers(igrid)%grid_obstacles), igrid
                WRITE(*, '()')

                DO i = 1, SIZE(my_obstacle_pointers(igrid)%grid_obstacles)

                    WRITE(*,'("Obstacle ", I0, ":")') my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%iobst
                    WRITE(*,'("x/y/z/r = ", 4F12.6)') my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%x, &
                    my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%y, &
                    my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%z, &
                    my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%radius
                    WRITE(*, '()')

                END DO

            END DO

            IF (myid /= numprocs - 1) THEN
                CALL MPI_Send(dummy, 1, mglet_mpi_int, myid + 1, 900, &
                 MPI_COMM_WORLD)
            END IF
        END IF


        DEALLOCATE(counter_array)
        DEALLOCATE(grid_processed)
        DEALLOCATE(obstacles_src)
        DEALLOCATE(is_relevant_src)
        DEALLOCATE(obstacles_itm)
        DEALLOCATE(is_relevant_itm)

        ! TODO: remove barrier?
        CALL MPI_Barrier(MPI_COMM_WORLD)

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("READING OF OBSTACLES SUCCESSFULLY COMPLETED.")')
                WRITE(*, '()')
            END IF
        END IF

    END SUBROUTINE read_obstacles

    SUBROUTINE finish_obstacles()

        ! local variables
        INTEGER(intk) :: i

        IF (ALLOCATED(my_obstacles)) DEALLOCATE(my_obstacles)
        IF (ALLOCATED(my_obstacle_pointers)) THEN
            DO i = 1, SIZE(my_obstacle_pointers)
                IF (ALLOCATED(my_obstacle_pointers(i)%grid_obstacles)) THEN
                    DEALLOCATE(my_obstacle_pointers(i)%grid_obstacles)
                END IF
            END DO
            DEALLOCATE(my_obstacle_pointers)
        END IF


    END SUBROUTINE finish_obstacles

    LOGICAL FUNCTION in_zone(this, igrid, overlap_f) result(res)

        ! subroutine arguments
        CLASS(obstacle_t), INTENT(in) :: this
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), OPTIONAL, INTENT(in) :: overlap_f

        ! local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: of = 0.0

        IF (PRESENT(overlap_f)) THEN
            of = overlap_f
        END IF

        res = .FALSE.

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        IF(maxx + (maxx - minx) * of + this%radius + EPSILON(maxx) < this%x) THEN
            RETURN
        END IF

        IF(minx - (maxx - minx) * of - this%radius - EPSILON(minx) > this%x) THEN
            RETURN
        END IF

        IF(maxy + (maxy - miny) * of + this%radius + EPSILON(maxy) < this%y) THEN
            RETURN
        END IF

        IF(miny - (maxy - miny) * of - this%radius - EPSILON(miny) > this%y) THEN
            RETURN
        END IF

        IF(maxz + (maxz - minz) * of + this%radius + EPSILON(maxz) < this%z) THEN
            RETURN
        END IF

        IF(minz - (maxz - minz) * of - this%radius - EPSILON(minz) > this%z) THEN
            RETURN
        END IF

        res = .TRUE.

    END FUNCTION in_zone

END MODULE particle_obstacles_mod