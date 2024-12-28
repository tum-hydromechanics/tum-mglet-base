MODULE particle_obstacles_mod

    ! This module is responsible for:
    ! Reading and storage of a sphere pack (obstacles that make up a porous domain)

    ! TODO: optimize this module

    USE MPI_f08
    USE comms_mod
    USE grids_mod
    USE field_mod
    USE fields_mod
    USE utils_mod

    USE particle_utils_mod

    IMPLICIT NONE

    ! TODO: make an abstract parent class
    TYPE :: obstacle_t

        INTEGER(intk) :: iobst = -1

        REAL(realk) :: x = 0.0
        REAL(realk) :: y = 0.0
        REAL(realk) :: z = 0.0

        REAL(realk) :: radius = 0.0

        CONTAINS
            PROCEDURE :: in_grid_zone

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
    REAL(realk), PARAMETER :: min_space = 0.05

    ! ratio of intermediate (filling) obstacles over readius of regular obstacles
    REAL(realk), PARAMETER :: radius_ratio = 0.348

CONTAINS    !===================================

    ! CAUTION: each process reads all obstacles
    SUBROUTINE read_obstacles()

        ! local variables
        TYPE(obstacle_t), ALLOCATABLE :: obstacles_src(:), obstacles_itm(:) ! temporary obtscle lists
        !(src -> input/regular/source obstacles (read), itm -> intermediate/filling obstacles (generated))

        INTEGER(intk) :: unit, dict_len, iobst, igrid, h, i, j, k, counter, dummy
        INTEGER(intk) :: neighbours(26)
        INTEGER(intk), ALLOCATABLE :: proc_neigbhours(:) ! array to store all grids ONCE that are neighbours to any grid of this proc
        INTEGER(intk), ALLOCATABLE :: counter_array(:) ! number of obstacles that is relevant on this process per grid

        LOGICAL :: dcycle, dexit
        LOGICAL, ALLOCATABLE :: grid_processed(:) ! array to indicate if a grid has been considered for a certain obstacle already
        LOGICAL, ALLOCATABLE :: is_relevant_src(:), is_relevant_itm(:) ! indicates if an obstacle is relevant on this process

        REAL(realk) :: dist
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        CHARACTER(12) :: dummy_char

        IF (.NOT. dread_obstacles) THEN
            ALLOCATE(my_obstacles(0))
            RETURN
        END IF

        INQUIRE(file = 'ObstaclesDict.txt', exist = dread_obstacles)

        IF (.NOT. dread_obstacles) THEN
            WRITE(*, *) "ERROR: No file for reading obstacles detected!"
            CALL errr(__FILE__,__LINE__)
        END IF

        ALLOCATE(grid_processed(ngrid))

        grid_processed = .FALSE.
        counter = 0
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_neighbours(neighbours, igrid)
            DO j = 1, 26
                IF (neighbours(j) < 1 .OR. neighbours(j) > ngrid) THEN
                    CYCLE
                ELSEIF (idprocofgrd(neighbours(j)) /= myid .AND. .NOT. grid_processed(neighbours(j))) THEN
                    grid_processed(neighbours(j)) = .TRUE.
                    counter = counter + 1
                END IF
            END DO
        END DO

        ALLOCATE(proc_neigbhours(counter))

        counter = 1
        DO i = 1, ngrid
            IF (grid_processed(i)) THEN
                IF (counter > SIZE(proc_neigbhours)) THEN
                    CALL errr(__FILE__,__LINE__)
                END IF
                proc_neigbhours(counter) = i
                counter = counter + 1
            END IF
        END DO

        ALLOCATE(counter_array(ngrid))
        counter_array = 0
        ALLOCATE(my_obstacle_pointers(ngrid))

        OPEN(newunit = unit, file = 'ObstaclesDict.txt', status = 'OLD', action = 'READ')

        READ(unit, fmt = *) dummy_char, dummy_char, dummy_char, dict_len

        IF (dict_len < 1) THEN
            WRITE(*, '("ERROR in particle_obstacles_mod: Number of Obstacles given in ObstaclesDict.txt must be an integer larger than 0.")')
            WRITE(*, '("If no Obstacles should be registered, set read_obst to FALSE.")')
            WRITE(*, '()')
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("READING ", I0, " OBSTACLES:")') dict_len
                WRITE(*, '()')
            END IF
        END IF

        ALLOCATE(obstacles_src(dict_len))
        ALLOCATE(is_relevant_src(dict_len))
        is_relevant_src = .FALSE.

        DO h = 1, dict_len

            READ(unit, fmt = *) dummy_char

            READ(unit, fmt = *) obstacles_src(h)%radius, obstacles_src(h)%x, obstacles_src(h)%y, &
             obstacles_src(h)%z

            IF (myid == 0) THEN
                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'("Obstacle read: x/y/z/r = ", 4F12.6)') obstacles_src(h)%x, obstacles_src(h)%y, &
                     obstacles_src(h)%z, obstacles_src(h)%radius
                END IF
            END IF

            obstacles_src(h)%iobst = h

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
                IF (dist < (obstacles_src(i)%radius + obstacles_src(j)%radius + min_space + EPSILON(obstacles_src(i)%radius))) THEN
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

                IF (dist < (obstacles_src(i)%radius + obstacles_src(j)%radius + min_space + EPSILON(obstacles_src(i)%radius))) THEN

                    obstacles_itm(counter)%iobst = counter + dict_len
                    obstacles_itm(counter)%x = obstacles_src(i)%x + (obstacles_src(j)%x - obstacles_src(i)%x) * obstacles_src(i)%radius / dist
                    obstacles_itm(counter)%y = obstacles_src(i)%y + (obstacles_src(j)%y - obstacles_src(i)%y) * obstacles_src(i)%radius / dist
                    obstacles_itm(counter)%z = obstacles_src(i)%z + (obstacles_src(j)%z - obstacles_src(i)%z) * obstacles_src(i)%radius / dist
                    obstacles_itm(counter)%radius = radius_ratio * MAX(obstacles_src(i)%radius, obstacles_src(j)%radius)

                    counter = counter + 1

                END IF

            END DO
        END DO

        ! count regular/source obstacles that are relevant for this process and determine for which grids they are relevant
        counter = 0
        DO h = 1, dict_len

                DO i = 1, nmygrids
                    igrid = mygrids(i)
                    IF (obstacles_src(h)%in_grid_zone(igrid)) THEN
                        counter_array(igrid) = counter_array(igrid) + 1
                        is_relevant_src(h) = .TRUE.
                    END IF
                END DO

                DO i = 1, SIZE(proc_neigbhours)
                    igrid = proc_neigbhours(i)
                    IF (obstacles_src(h)%in_grid_zone(igrid)) THEN
                        counter_array(igrid) = counter_array(igrid) + 1
                        is_relevant_src(h) = .TRUE.
                    END IF
                END DO

                IF (is_relevant_src(h)) THEN
                    counter = counter + 1
                END IF

        END DO

        ! count intermediate/filling obstacles that are relevant for this process and determine for which grids they are relevant
        DO h = 1, SIZE(obstacles_itm)

            DO i = 1, nmygrids
                igrid = mygrids(i)
                IF (obstacles_itm(h)%in_grid_zone(igrid)) THEN
                    counter_array(igrid) = counter_array(igrid) + 1
                    is_relevant_itm(h) = .TRUE.
                END IF
            END DO

            DO i = 1, SIZE(proc_neigbhours)
                igrid = proc_neigbhours(i)
                IF (obstacles_itm(h)%in_grid_zone(igrid)) THEN
                    counter_array(igrid) = counter_array(igrid) + 1
                    is_relevant_itm(h) = .TRUE.
                END IF
            END DO

            IF (is_relevant_itm(h)) THEN
                counter = counter + 1
            END IF

        END DO

        DO igrid = 1, ngrid
            ALLOCATE(my_obstacle_pointers(igrid)%grid_obstacles(counter_array(igrid)))
        END DO

        ALLOCATE(my_obstacles(counter))

        counter_array = 1
        counter = 1

        ! copy relevant input/source obstacles for this process into heap memory
        DO h = 1, dict_len

            IF (is_relevant_src(h)) THEN

                my_obstacles(counter) = obstacles_src(h)

                DO i = 1, nmygrids
                    igrid = mygrids(i)
                    IF (obstacles_src(h)%in_grid_zone(igrid)) THEN
                        IF (counter_array(igrid) > SIZE(my_obstacle_pointers(igrid)%grid_obstacles)) THEN
                            CALL errr(__FILE__, __LINE__)
                        END IF
                        my_obstacle_pointers(igrid)%grid_obstacles(counter_array(igrid)) = counter
                        counter_array(igrid) = counter_array(igrid) + 1
                    END IF
                END DO

                DO i = 1, SIZE(proc_neigbhours)
                    igrid = proc_neigbhours(i)
                    IF (obstacles_src(h)%in_grid_zone(igrid)) THEN
                        IF (counter_array(igrid) > SIZE(my_obstacle_pointers(igrid)%grid_obstacles)) THEN
                            CALL errr(__FILE__, __LINE__)
                        END IF
                        my_obstacle_pointers(igrid)%grid_obstacles(counter_array(igrid)) = counter
                        counter_array(igrid) = counter_array(igrid) + 1
                    END IF
                END DO

                counter = counter + 1

            END IF

        END DO

        ! copy relevant intermediate/filling obstacles for this process into heap memory
        DO h = 1, SIZE(obstacles_itm)

            IF (is_relevant_itm(h)) THEN

                my_obstacles(counter) = obstacles_itm(h)

                DO i = 1, nmygrids
                    igrid = mygrids(i)
                    IF (obstacles_itm(h)%in_grid_zone(igrid)) THEN
                        IF (counter_array(igrid) > SIZE(my_obstacle_pointers(igrid)%grid_obstacles)) THEN
                            CALL errr(__FILE__, __LINE__)
                        END IF
                        my_obstacle_pointers(igrid)%grid_obstacles(counter_array(igrid)) = counter
                        counter_array(igrid) = counter_array(igrid) + 1
                    END IF
                END DO

                DO i = 1, SIZE(proc_neigbhours)
                    igrid = proc_neigbhours(i)
                    IF (obstacles_itm(h)%in_grid_zone(igrid)) THEN
                        IF (counter_array(igrid) > SIZE(my_obstacle_pointers(igrid)%grid_obstacles)) THEN
                            CALL errr(__FILE__, __LINE__)
                        END IF
                        my_obstacle_pointers(igrid)%grid_obstacles(counter_array(igrid)) = counter
                        counter_array(igrid) = counter_array(igrid) + 1
                    END IF
                END DO

                counter = counter + 1

            END IF

        END DO

        IF (TRIM(particle_terminal) == "verbose") THEN
            IF (myid /= 0) THEN
                CALL MPI_Recv(dummy, 1, mglet_mpi_int, myid - 1, 900, &
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE)
            END IF

            WRITE(*, '("Obstacels on Process ", I0, ":")') myid
            WRITE(*, '()')

            DO igrid = 1, ngrid

                IF ((SIZE(my_obstacle_pointers(igrid)%grid_obstacles) > 0 .AND. TRIM(particle_terminal) == "normal") &
                 .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*,'(I0, " Obstacles registered for Grid ", I0, ".")') SIZE(my_obstacle_pointers(igrid)%grid_obstacles), igrid
                    WRITE(*, '()')
                END IF

                IF (TRIM(particle_terminal) == "verbose") THEN
                    DO i = 1, SIZE(my_obstacle_pointers(igrid)%grid_obstacles)
                        WRITE(*,'("Obstacle ", I0, ":")') my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%iobst
                        WRITE(*,'("x/y/z/r = ", 4F12.6)') my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%x, &
                        my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%y, &
                        my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%z, &
                        my_obstacles(my_obstacle_pointers(igrid)%grid_obstacles(i))%radius
                        WRITE(*, '()')
                    END DO
                END IF

            END DO

            IF (myid /= numprocs - 1) THEN
                CALL MPI_Send(dummy, 1, mglet_mpi_int, myid + 1, 900, &
                 MPI_COMM_WORLD)
            END IF
        END IF

        ! the following vtk output is optional and can be removed
        CALL write_obstacles()
        CALL write_grids(mygrids, nmygrids, "   ")
        CALL write_grids(proc_neigbhours, SIZE(proc_neigbhours), "nbr")

        IF (ALLOCATED(proc_neigbhours)) DEALLOCATE(proc_neigbhours)
        IF (ALLOCATED(counter_array)) DEALLOCATE(counter_array)
        IF (ALLOCATED(grid_processed)) DEALLOCATE(grid_processed)
        IF (ALLOCATED(obstacles_src)) DEALLOCATE(obstacles_src)
        IF (ALLOCATED(is_relevant_src)) DEALLOCATE(is_relevant_src)
        IF (ALLOCATED(obstacles_itm)) DEALLOCATE(obstacles_itm)
        IF (ALLOCATED(is_relevant_itm)) DEALLOCATE(is_relevant_itm)

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

    ! TODO: merge this with particle snapshots writing or use library
    ! write obstacles vtk to validate if they have been registered properly
    SUBROUTINE write_obstacles()

        ! local variables
        INTEGER(intk) :: h, i, j, k, unit, igrid
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        CHARACTER(len = mglet_filename_max) :: filename, my_nobst_char, igrid_char

        IF (myid == 0) THEN
            CALL create_directory("Particle_Obstacles") ! ! ! realtive to working directory ! ! !
        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

        ! write obstacles (one file per proc)
        WRITE(filename, '("Particle_Obstacles/obstacles_proc", I0, ".vtp")') myid

        WRITE(my_nobst_char, '(I0)') SIZE(my_obstacles)

        OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

        WRITE(unit, '(A)') '<?xml version="1.0"?>'
        WRITE(unit, '(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
        WRITE(unit, '(A)') '  <PolyData>'
        WRITE(unit, '(A)') '    <Piece NumberOfPoints="' // TRIM(my_nobst_char) // &
                              '" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
        WRITE(unit, '(A)') '      <PointData Name="obstacles">'
        WRITE(unit, '(A)') '        <DataArray type="Int32" format="ascii" NumberOfComponents="1" Name="iobst">'

        DO i = 1, SIZE(my_obstacles)
            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, '(I0)') my_obstacles(i)%iobst
        END DO

        WRITE(unit, '(A)') '        </DataArray>'
        WRITE(unit, '(A)') '        <DataArray type="Float32" format="ascii" NumberOfComponents="1" Name="radius">'

        DO i = 1, SIZE(my_obstacles)
            WRITE(unit, '("          ")', advance="no")
            WRITE(unit,  vtk_float_format) my_obstacles(i)%radius
        END DO

        WRITE(unit, '(A)') '        </DataArray>'
        WRITE(unit, '(A)') '      </PointData>'
        WRITE(unit, '(A)') '      <Points>'
        WRITE(unit, '(A)') '        <DataArray type="Float32" NumberOfComponents="3">'

        DO i = 1, SIZE(my_obstacles)
            WRITE(unit, '("        ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") my_obstacles(i)%x
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") my_obstacles(i)%y
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") my_obstacles(i)%z
        END DO

        WRITE(unit, '(A)') '        </DataArray>'
        WRITE(unit, '(A)') '      </Points>'
        WRITE(unit, '(A)') '    </Piece>'
        WRITE(unit, '(A)') '  </PolyData>'
        WRITE(unit, '(A)') '</VTKFile>'

        CLOSE(unit)

    END SUBROUTINE write_obstacles

    SUBROUTINE write_grids(grids, n, prefix)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: n
        INTEGER(intk), INTENT(in) :: grids(n)
        CHARACTER(len = 3), INTENT(in) :: prefix

        ! local variables
        INTEGER(intk) :: h, i, j, k, unit, igrid
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        CHARACTER(len = mglet_filename_max) :: filename, my_nobst_char, igrid_char

        ! also write all grids into VTK (one proc per file)
        WRITE(filename, '("Particle_Obstacles/", A, "grids_proc", I0, ".vtp")') TRIM(prefix), myid

        OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

        WRITE(unit, '(A)') '<?xml version="1.0"?>'
        WRITE(unit, '(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
        WRITE(unit, '(A)') '  <PolyData>'

        DO i = 1, n

            igrid = grids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            WRITE(igrid_char, '(I0)') igrid

            WRITE(unit, '(A)') '    <Piece Name="grid'  // TRIM(igrid_char) // &
             '" NumberOfPoints="8" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="6">'
            WRITE(unit, '(A)') '      <Points Name="grid_corners">'
            WRITE(unit, '(A)') '        <DataArray type="Float32" format="ascii" NumberOfComponents="3" Name="grid_corners">'

            ! dirty ...
            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") minx
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") miny
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") minz

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") minx
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") miny
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") maxz

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") minx
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") maxy
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") minz

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") minx
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") maxy
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") maxz

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") maxx
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") miny
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") minz

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") maxx
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") miny
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") maxz

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") maxx
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") maxy
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") minz

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, vtk_float_format, advance="no") maxx
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="no") maxy
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, vtk_float_format, advance="yes") maxz

            WRITE(unit, '(A)') '        </DataArray>'
            WRITE(unit, '(A)') '      </Points>'
            WRITE(unit, '(A)') '      <Polys Name="grid_faces">'
            WRITE(unit, '(A)') '        <DataArray type="Int32" Name="connectivity">'

            WRITE(unit, '("         ")', advance="no")
            WRITE(unit, '(4I2)', advance="yes") 0, 1, 3, 2

            WRITE(unit, '("         ")', advance="no")
            WRITE(unit, '(4I2)', advance="yes") 4, 5, 7, 6

            WRITE(unit, '("         ")', advance="no")
            WRITE(unit, '(4I2)', advance="yes") 0, 1, 5, 4

            WRITE(unit, '("         ")', advance="no")
            WRITE(unit, '(4I2)', advance="yes") 2, 3, 7, 6

            WRITE(unit, '("         ")', advance="no")
            WRITE(unit, '(4I2)', advance="yes") 0, 2, 6, 4

            WRITE(unit, '("         ")', advance="no")
            WRITE(unit, '(4I2)', advance="yes") 1, 3, 7, 5

            WRITE(unit, '(A)') '        </DataArray>'
            WRITE(unit, '(A)') '        <DataArray type="Int32" Name="offsets">'

            WRITE(unit, '("         ")', advance="no")
            WRITE(unit, '(6I3)', advance="yes") 4, 8, 12, 16, 20, 24

            WRITE(unit, '(A)') '        </DataArray>'
            WRITE(unit, '(A)') '      </Polys>'
            WRITE(unit, '(A)') '    </Piece>'

        END DO

        WRITE(unit, '(A)') '  </PolyData>'
        WRITE(unit, '(A)') '</VTKFile>'

        CLOSE(unit)

    END SUBROUTINE write_grids

    LOGICAL FUNCTION in_grid_zone(this, igrid, overlap) result(res)

        ! subroutine arguments
        CLASS(obstacle_t), INTENT(in) :: this
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), OPTIONAL, INTENT(in) :: overlap

        ! local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: ol = 0.0

        IF (PRESENT(overlap)) THEN
            ol = overlap
        END IF

        res = .FALSE.

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        IF(maxx + ol + this%radius + EPSILON(maxx) < this%x) THEN
            RETURN
        END IF

        IF(minx - ol - this%radius - EPSILON(minx) > this%x) THEN
            RETURN
        END IF

        IF(maxy + ol + this%radius + EPSILON(maxy) < this%y) THEN
            RETURN
        END IF

        IF(miny - ol - this%radius - EPSILON(miny) > this%y) THEN
            RETURN
        END IF

        IF(maxz + ol + this%radius + EPSILON(maxz) < this%z) THEN
            RETURN
        END IF

        IF(minz - ol - this%radius - EPSILON(minz) > this%z) THEN
            RETURN
        END IF

        res = .TRUE.

    END FUNCTION in_grid_zone

END MODULE particle_obstacles_mod