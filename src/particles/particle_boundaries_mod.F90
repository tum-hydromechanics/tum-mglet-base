MODULE particle_boundaries_mod
    USE precision_mod, ONLY: realk, intk
    USE core_mod
    USE particle_list_mod

    IMPLICIT NONE

    INTEGER(intk), PARAMETER :: facelist(4,26) = RESHAPE((/ &
        1, 1, 0, 0, &
        1, 2, 0, 0, &
        1, 3, 0, 0, &
        1, 4, 0, 0, &
        1, 5, 0, 0, &
        1, 6, 0, 0, &
        2, 1, 3, 0, &
        2, 1, 4, 0, &
        2, 1, 5, 0, &
        2, 1, 6, 0, &
        2, 2, 3, 0, &
        2, 2, 4, 0, &
        2, 2, 5, 0, &
        2, 2, 6, 0, &
        2, 3, 5, 0, &
        2, 3, 6, 0, &
        2, 4, 5, 0, &
        2, 4, 6, 0, &
        3, 1, 3, 5, &
        3, 1, 3, 6, &
        3, 1, 4, 5, &
        3, 1, 4, 6, &
        3, 2, 3, 5, &
        3, 2, 3, 6, &
        3, 2, 4, 5, &
        3, 2, 4, 6 /), SHAPE(facelist))

    TYPE :: obstacle_t

        REAL(realk) :: radius
        REAL(realk) :: center(3)
        INTEGER(intk) :: particle_boundary_type

    END TYPE obstacle_t

    TYPE :: obstacle_list_t ! every grid has one list! (not an exact analog to the particle list)

        INTEGER(intk) :: igrid
        TYPE(obstacle_t), ALLOCATABLE :: obstacles(:)

    END TYPE obstacle_list_t

    TYPE :: particle_boundaries_t

        INTEGER(intk), ALLOCATABLE :: face_neighbours(:, :)

        REAL(realk), ALLOCATABLE :: face_normals(:, :, :)

        TYPE(obstacle_list_t), ALLOCATABLE :: obstacle_lists(:) ! every grid has one list! (not an exact analog to the particle list)

    END TYPE particle_boundaries_t

    TYPE(particle_boundaries_t) :: particle_boundaries

    CONTAINS

    SUBROUTINE init_particle_boundaries()

        ! local variables
        INTEGER(intk) :: igrid, ibocd, iface, jface, i, j
        INTEGER(intk) :: neighbours(26)
        ! to store which boundary surfaces (faces 1-6) that define a lower order face (7-26) are of connect (or periodic) type:
        INTEGER(intk) :: connect_faces(4)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, n_minx, n_maxx, n_miny, n_maxy, n_minz, n_maxz, magnitude
        CHARACTER(len=8) :: ctyp
        LOGICAL :: found

        ALLOCATE(particle_boundaries%face_normals(3, 26, ngrid))
        particle_boundaries%face_normals = 0.0

        ALLOCATE(particle_boundaries%face_neighbours(26, ngrid))

        ibocd = 3
        DO igrid = 1, ngrid

            CALL get_neighbours(neighbours, igrid)

            DO iface = 1, 6

                CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)

                particle_boundaries%face_neighbours(iface, igrid) = neighbours(iface)

                SELECT CASE(iface)
                    CASE(1)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%face_normals(1, 1, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(2)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%face_normals(1, 2, igrid) = -1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(3)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%face_normals(2, 3, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(4)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%face_normals(2, 4, igrid) = -1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(5)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%face_normals(3, 5, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(6)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%face_normals(3, 6, igrid) = -1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                END SELECT

            END DO

            DO iface = 7, 26

                connect_faces = 0

                DO i = 2, 4
                    IF (facelist(i, iface) == 0) THEN
                        CONTINUE
                    ELSE

                        CALL get_bc_ctyp(ctyp, ibocd, facelist(i, iface), igrid)

                        IF (ctyp == 'CON' .OR. ctyp == 'PER') THEN
                            connect_faces(1) = connect_faces(1) + 1
                            connect_faces(i) = facelist(i, iface)
                        END IF

                        DO j = 1, 3
                            particle_boundaries%face_normals(j, iface, igrid) = &
                            particle_boundaries%face_normals(j, iface, igrid) + particle_boundaries%face_normals(j, facelist(i, iface), igrid)
                        END DO

                    END IF
                END DO

                particle_boundaries%face_neighbours(iface, igrid) = igrid
                found = .FALSE.
                DO jface = 1, 26
                    IF (facelist(1, jface) /= connect_faces(1)) THEN
                        CONTINUE
                    ELSE
                        found = .TRUE.
                        DO i = 2, 1 + facelist(1, jface)
                            IF (facelist(i, jface) /= connect_faces(2) &
                             .AND. facelist(i, jface) /= connect_faces(3) .AND. facelist(i, jface) /= connect_faces(4)) THEN
                                found = .FALSE.
                                EXIT
                            ELSE
                                CONTINUE
                            END IF
                        END DO
                        IF (found .eqv. .TRUE.) THEN
                            particle_boundaries%face_neighbours(iface, igrid) = jface
                        END IF
                    END IF
                END DO

                magnitude = SQRT(particle_boundaries%face_normals(1, iface, igrid)**2 + &
                 particle_boundaries%face_normals(2, iface, igrid)**2 + &
                 particle_boundaries%face_normals(3, iface, igrid)**2)

                DO j = 1, 3
                    particle_boundaries%face_normals(j, iface, igrid) = particle_boundaries%face_normals(j, iface, igrid) / magnitude
                END DO

            END DO
        END DO

        DO igrid = 1, ngrid

            WRITE(*, *) "------ Boundaries ------ grid:   ", igrid
            CALL get_bc_ctyp(ctyp, ibocd, 1, igrid)
            WRITE(*, *) "FRONT:                ", ctyp
            CALL get_bc_ctyp(ctyp, ibocd, 2, igrid)
            WRITE(*, *) "BACK:                 ", ctyp
            CALL get_bc_ctyp(ctyp, ibocd, 3, igrid)
            WRITE(*, *) "RIGHT:                ", ctyp
            CALL get_bc_ctyp(ctyp, ibocd, 4, igrid)
            WRITE(*, *) "LEFT:                 ", ctyp
            CALL get_bc_ctyp(ctyp, ibocd, 5, igrid)
            WRITE(*, *) "BOTTOM:               ", ctyp
            CALL get_bc_ctyp(ctyp, ibocd, 6, igrid)
            WRITE(*, *) "TOP:                  ", ctyp

            WRITE(*, *) "Faces:"
            DO iface = 1, 26
                WRITE(*, *) "Face:                 ", iface
                WRITE(*, *) "Neigbhour grid:       ", particle_boundaries%face_neighbours(iface, igrid)
                WRITE(*, *) "Normal vector (n1):   ", particle_boundaries%face_normals(1, iface, igrid)
                WRITE(*, *) "Normal vector (n2):   ", particle_boundaries%face_normals(2, iface, igrid)
                WRITE(*, *) "Normal vector (n3):   ", particle_boundaries%face_normals(3, iface, igrid)
            END DO

        END DO

    END SUBROUTINE init_particle_boundaries

    !-----------------------------------

    SUBROUTINE finish_particle_boundaries()

        DEALLOCATE(particle_boundaries%face_normals)
        DEALLOCATE(particle_boundaries%obstacle_lists)

    END SUBROUTINE finish_particle_boundaries

    !-----------------------------------

    SUBROUTINE move_particle(particle, dx, dy, dz)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        REAL(realk), INTENT(in) :: dx, dy, dz

        ! local variables
        INTEGER(intk) :: temp_grid, iface, grid_bc, destproc
        INTEGER(intk) :: neighbours(26)
        REAL(realk) :: x, y, z
        REAL(realk) :: dx_from_here, dy_from_here, dz_from_here
        REAL(realk) :: epsilon
        REAL(realk) :: n1, n2, n3

        epsilon = SQRT(dx**(2) + dy**(2) + dz**(2)) / 10**3

        temp_grid = particle%igrid

        x = particle%x
        y = particle%y
        z = particle%z

        dx_from_here = dx
        dy_from_here = dy
        dz_from_here = dz

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                WRITE(*, *) "MOVE PARTICLE: ----------------------START-----------------------"
                WRITE(*, '()')
        END SELECT

        DO WHILE (epsilon < ABS(dx_from_here) .OR. epsilon < ABS(dy_from_here) .OR. epsilon < ABS(dz_from_here))

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, *) "Way to go:"
                    WRITE(*, *) "dx/dy/dz:", dx_from_here, dy_from_here, dz_from_here
                    WRITE(*, '()')
            END SELECT

            CALL move_to_boundary(temp_grid, x, y, z, dx_from_here, dy_from_here, dz_from_here, iface)

            CALL reflect_at_boundary(dx_from_here, dy_from_here, dz_from_here, &
             particle_boundaries%face_normals(1, iface, temp_grid), &
             particle_boundaries%face_normals(2, iface, temp_grid), &
             particle_boundaries%face_normals(3, iface, temp_grid))

            temp_grid = particle_boundaries%face_neighbours(iface, temp_grid)

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, *) "Current Coordinates:"
                    WRITE(*, *) "x/y/z:", x, y, z
                    WRITE(*, '()')
            END SELECT

        END DO

        !IF ( destproc > numprocs .OR. destproc < 0 ) THEN
        !    WRITE(*,*) 'Ill-addressed particle to proc', destproc
        !    CALL errr(__FILE__, __LINE__)
        !END IF

        !IF (destproc /= myid) THEN
        !    particle%state = 4
        !END IF

        particle%x = x
        particle%y = y
        particle%z = z

        CALL update_particle_cell(particle)

    END SUBROUTINE move_particle

    !-----------------------------------

    ! This subroutine only considers grids on the same level
    ! CAUTION: Here, igrid refers to the grid the particle coordinates are currently on and of which the boundaries are relevant.
    ! This might NOT be particle%igrid, which is used to deduce the velocity from.
    SUBROUTINE move_to_boundary(igrid, x, y, z, dx, dy, dz, iface)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        INTEGER(intk), INTENT(out) :: iface

        !local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: lx, ly, lz, rx, ry, rz, dx_to_b, dy_to_b, dz_to_b

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)


        IF (x < minx) THEN
                        SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (px < minx)! Setting px to minx."
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (px < minx)! Setting px to minx."
                    WRITE(*, '()')
            END SELECT
            x = minx
        END IF

        IF (x > maxx) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (px < minx)! Setting px to minx."
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (px > maxx)! Setting px to maxx."
                    WRITE(*, '()')
            END SELECT
            x = maxx
        END IF

        IF (y < miny) THEN
                        SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (py < miny)! Setting py to miny."
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (py < miny)! Setting py to miny."
                    WRITE(*, '()')
            END SELECT
            y = miny
        END IF

        IF (y > maxy) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (py < miny)! Setting py to miny."
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (py > maxy)! Setting py to maxy."
                    WRITE(*, '()')
            END SELECT
            y = maxy
        END IF

        IF (z < minz) THEN
                        SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (pz < minz)! Setting pz to minz."
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (pz < minz)! Setting pz to minz."
                    WRITE(*, '()')
            END SELECT
            z = minz
        END IF

        IF (z > maxz) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (pz < minz)! Setting pz to minz."
                    WRITE(*, '()')
                CASE ("verbose")
                    WRITE(*, *) "WARNING: In move_to_boundary: Particle outside of grid boundaries (pz > maxz)! Setting pz to maxz."
                    WRITE(*, '()')
            END SELECT
            z = maxz
        END IF

        IF (dx < 0) THEN
            lx = (minx - x)
            IF(lx == 0) THEN
                CALL get_current_face(igrid, x, y, z, iface)
                RETURN
            END IF
            rx = dx / lx
        ELSEIF (0 < dx) THEN
            lx = (maxx - x)
            IF(lx == 0) THEN
                CALL get_current_face(igrid, x, y, z, iface)
                RETURN
            END IF
            rx = dx / lx
        ELSE
            rx = 0.0_realk
        END IF

        IF (dy < 0) THEN
            ly = (miny - y)
            IF(ly == 0) THEN
                CALL get_current_face(igrid, x, y, z, iface)
                RETURN
            END IF
            ry = dy / ly
        ELSEIF (0 < dy) THEN
            ly = (maxy - y)
            IF(ly == 0) THEN
                CALL get_current_face(igrid, x, y, z, iface)
                RETURN
            END IF
            ry = dy / ly
        ELSE
            ry = 0.0_realk
        END IF

        IF (dz < 0) THEN
            lz = (minz - z)
            IF(lz == 0) THEN
                CALL get_current_face(igrid, x, y, z, iface)
                RETURN
            END IF
            rz = dz / lz
        ELSEIF (0 < dz) THEN
            lz = (maxz - z)
            IF(lz == 0) THEN
                CALL get_current_face(igrid, x, y, z, iface)
                RETURN
            END IF
            rz = dz / lz
        ELSE
            rz = 0.0_realk
        END IF

        IF (rx < 1.0_realk .AND. ry < 1.0_realk .AND. rz < 1.0_realk) THEN

            x = x + dx
            y = y + dy
            z = z + dz
            dx = 0
            dy = 0
            dz = 0

            RETURN

        END IF

        IF (dx < 0 .AND. ry <= rx .AND. rz <= rx) THEN

            dx_to_b = lx
            dy_to_b = (lx * dy/dx)
            dz_to_b = (lx * dz/dx)
            x = x + dx_to_b
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dx .AND. ry <= rx .AND. rz <= rx) THEN

            dx_to_b = lx
            dy_to_b = (lx * dy/dx)
            dz_to_b = (lx * dz/dx)
            x = x + dx_to_b
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (dy < 0 .AND. rx < ry .AND. rz <= ry) THEN

            dx_to_b = (ly * dx/dy)
            dy_to_b = ly
            dz_to_b = (ly * dz/dy)
            x = x + dx_to_b
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dy .AND. rx < ry .AND. rz <= ry) THEN

            dx_to_b = (ly * dx/dy)
            dy_to_b = ly
            dz_to_b = (ly * dz/dy)
            x = x + dx_to_b
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b


        ELSEIF (dz < 0 .AND. rx < rz .AND. ry < rz) THEN

            dx_to_b = (lz * dx/dz)
            dy_to_b = (lz * dy/dz)
            dz_to_b = lz
            x = x + dx_to_b
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dz .AND. rx < rz .AND. ry < rz) THEN

            dx_to_b = (lz * dx/dz)
            dy_to_b = (lz * dy/dz)
            dz_to_b = lz
            x = x + dx_to_b
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        END IF

        ! get face after x/y/z have (potentially) been altered
        CALL get_current_face(igrid, x, y, z, iface)

    END SUBROUTINE move_to_boundary

    SUBROUTINE get_current_face(igrid, x, y, z, iface)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(in) :: x, y, z
        INTEGER(intk), INTENT(out) :: iface

        ! local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        IF (x == minx) THEN !-------------------------------------------------------- low x
            IF (y == miny) THEN !--------------------------------------------- low y, low x
                IF (z == minz) THEN !---------------------------------- low z, low y, low x
                    iface = 19
                ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, low y, low x
                    iface = 7
                ELSEIF (z == maxz) THEN !----------------------------- high z, low y, low x
                    iface = 20
                END IF
            ELSEIF (miny < y .AND. y < maxy) THEN !--------------------------- mid y, low x
                IF (z == minz) THEN !---------------------------------- low z, mid y, low x
                    iface = 9
                ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, mid y, low x
                    iface = 1
                ELSEIF (z == maxz) THEN !----------------------------- high z, mid y, low x
                    iface = 10
                END IF
            ELSEIF (y == maxy) THEN !---------------------------------------- high y, low x
                IF (z == minz) THEN !--------------------------------- low z, high y, low x
                    iface = 21
                ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, high y, low x
                    iface = 8
                ELSEIF (z == maxz) THEN !---------------------------- high z, high y, low x
                    iface = 22
                END IF
            END IF
        ELSEIF (minx < x .AND. x < maxx) THEN !-------------------------------------- mid x
            IF (y == miny) THEN !--------------------------------------------- low y, mid x
                IF (z == minz) THEN !---------------------------------- low z, low y, mid x
                    iface = 15
                ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, low y, mid x
                    iface = 3
                ELSEIF (z == maxz) THEN !----------------------------- high z, low y, mid x
                    iface = 16
                END IF
            ELSEIF (miny < y .AND. y < maxy) THEN !--------------------------- mid y, mid x
                IF (z == minz) THEN !---------------------------------- low z, mid y, mid x
                    iface = 5
                ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, mid y, mid x
                    iface = 0
                ELSEIF (z == maxz) THEN !----------------------------- high z, mid y, mid x
                    iface = 6
                END IF
            ELSEIF (y == maxy) THEN !---------------------------------------- high y, mid x
                IF (z == minz) THEN !--------------------------------- low z, high y, mid x
                    iface = 17
                ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, high y, mid x
                    iface = 4
                ELSEIF (z == maxz) THEN !---------------------------- high z, high y, mid x
                    iface = 18
                END IF
            END IF
        ELSEIF (x == maxx) THEN !--------------------------------------------------- high x
            IF (y == miny) THEN !-------------------------------------------- low y, high x
                IF (z == minz) THEN !--------------------------------- low z, low y, high x
                    iface = 23
                ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, low y, high x
                    iface = 11
                ELSEIF (z == maxz) THEN !---------------------------- high z, low y, high x
                    iface = 24
                END IF
            ELSEIF (miny < y .AND. y < maxy) THEN !-------------------------- mid y, high x
                IF (z == minz) THEN !--------------------------------- low z, mid y, high x
                    iface = 13
                ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, mid y, high x
                    iface = 2
                ELSEIF (z == maxz) THEN !---------------------------- high z, mid y, high x
                    iface = 14
                END IF
            ELSEIF (y == maxy) THEN !--------------------------------------- high y, high x
                IF (z == minz) THEN !-------------------------------- low z, high y, high x
                    iface = 25
                ELSEIF (minz < z .AND. z < maxz) THEN !-------------- mid z, high y, high x
                    iface = 12
                ELSEIF (z == maxz) THEN !--------------------------- high z, high y, high x
                    iface = 26
                END IF
            END IF
        END IF

    END SUBROUTINE get_current_face

    SUBROUTINE apply_periodic_boundary(x, y, z, oldgrid, newgrid, iface)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: oldgrid, newgrid, iface
        REAL(realk), INTENT(inout) :: x, y, z

        ! local variables
        REAL(realk) :: old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, &
         new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz

        CALL get_bbox(old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, oldgrid)
        CALL get_bbox(new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz, newgrid)

        SELECT CASE(iface)
            CASE(1)
                x = new_maxx + (x - old_minx)
            CASE(2)
                x = new_minx + (x - old_maxx)
            CASE(3)
                y = new_maxy + (y - old_miny)
            CASE(4)
                y = new_miny + (y - old_maxy)
            CASE(5)
                z = new_maxz + (z - old_minz)
            CASE(6)
                z = new_minz + (z - old_maxz)
        END SELECT

    END SUBROUTINE apply_periodic_boundary

    SUBROUTINE reflect_at_boundary(dx, dy, dz, n1, n2, n3)

        ! Presumption: Particle is already exactily on the boundary!

        ! subroutine arguments
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(in) :: n1, n2, n3 ! normal vector components of the surface the particle is reflected from

        ! local variables
        REAL(realk) :: dot_product

        dot_product = n1 * dx + n2 * dy + n3 * dz

        dx = dx - 2 * dot_product * n1
        dy = dy - 2 * dot_product * n2
        dz = dz - 2 * dot_product * n3

    END SUBROUTINE reflect_at_boundary

    SUBROUTINE get_normal_vector(iface, n1, n2, n3)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: iface
        REAL(realk), INTENT(out) :: n1, n2, n3 ! normal vector components

        ! local variables
        LOGICAL :: found

        SELECT CASE(iface)
            CASE(1)
                n1 = 1
                n2 = 0
                n3 = 0
            CASE(2)
                n1 = -1
                n2 = 0
                n3 = 0
            CASE(3)
                n1 = 0
                n2 = 1
                n3 = 0
            CASE(4)
                n1 = 0
                n2 = -1
                n3 = 0
            CASE(5)
                n1 = 0
                n2 = 0
                n3 = 1
            CASE(6)
                n1 = 0
                n2 = 0
                n3 = -1
            CASE(7)
                n1 = 1/SQRT(2.0)
                n2 = 1/SQRT(2.0)
                n3 = 0
            CASE(8)
                n1 = 1/SQRT(2.0)
                n2 = -1/SQRT(2.0)
                n3 = 0
            CASE(9)
                n1 = 1/SQRT(2.0)
                n2 = 0
                n3 = 1/SQRT(2.0)
            CASE(10)
                n1 = 1/SQRT(2.0)
                n2 = 0
                n3 = -1/SQRT(2.0)
            CASE(11)
                n1 = -1/SQRT(2.0)
                n2 = 1/SQRT(2.0)
                n3 = 0
            CASE(12)
                n1 = -1/SQRT(2.0)
                n2 = -1/SQRT(2.0)
                n3 = 0
            CASE(13)
                n1 = -1/SQRT(2.0)
                n2 = 0
                n3 = 1/SQRT(2.0)
            CASE(14)
                n1 = -1/SQRT(2.0)
                n2 = 0
                n3 = -1/SQRT(2.0)
            CASE(15)
                n1 = 0
                n2 = 1/SQRT(2.0)
                n3 = 1/SQRT(2.0)
            CASE(16)
                n1 = 0
                n2 = 1/SQRT(2.0)
                n3 = -1/SQRT(2.0)
            CASE(17)
                n1 = 0
                n2 = -1/SQRT(2.0)
                n3 = 1/SQRT(2.0)
            CASE(18)
                n1 = 0
                n2 = -1/SQRT(2.0)
                n3 = -1/SQRT(2.0)
            CASE(19)
                n1 = 1/SQRT(3.0)
                n2 = 1/SQRT(3.0)
                n3 = 1/SQRT(3.0)
            CASE(20)
                n1 = 1/SQRT(3.0)
                n2 = 1/SQRT(3.0)
                n3 = -1/SQRT(3.0)
            CASE(21)
                n1 = 1/SQRT(3.0)
                n2 = -1/SQRT(3.0)
                n3 = 1/SQRT(3.0)
            CASE(22)
                n1 = 1/SQRT(3.0)
                n2 = -1/SQRT(3.0)
                n3 = -1/SQRT(3.0)
            CASE(23)
                n1 = -1/SQRT(3.0)
                n2 = 1/SQRT(3.0)
                n3 = 1/SQRT(3.0)
            CASE(24)
                n1 = -1/SQRT(3.0)
                n2 = 1/SQRT(3.0)
                n3 = -1/SQRT(3.0)
            CASE(25)
                n1 = -1/SQRT(3.0)
                n2 = -1/SQRT(3.0)
                n3 = 1/SQRT(3.0)
            CASE(26)
                n1 = -1/SQRT(3.0)
                n2 = -1/SQRT(3.0)
                n3 = -1/SQRT(3.0)
        END SELECT

    END SUBROUTINE get_normal_vector

    SUBROUTINE update_coordinates(oldgrid, newgrid, x, y, z)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: oldgrid, newgrid
        REAL(realk), INTENT(inout) :: x, y, z

        ! local variables
        REAL(realk) :: old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, &
         new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz
        LOGICAL :: passed_pb = .FALSE.

        CALL get_bbox(old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, oldgrid)
        CALL get_bbox(new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz, newgrid)

        IF (x < new_minx) THEN
            x = new_maxx - ABS(x - old_minx)
            passed_pb = .TRUE.
        END IF

        IF (new_maxx < x) THEN
            x = new_minx + ABS(x - old_maxx)
            passed_pb = .TRUE.
        END IF

        IF (y < new_miny) THEN
            y = new_maxy - ABS(y - old_miny)
            passed_pb = .TRUE.
        END IF

        IF (new_maxy < y) THEN
            y = new_miny + ABS(y - old_maxy)
            passed_pb = .TRUE.
        END IF

        IF (z < new_minz) THEN
            z = new_maxz - ABS(z - old_minz)
            passed_pb = .TRUE.
        END IF

        IF (new_maxz < z) THEN
            z = new_minz + ABS(z - old_maxz)
            passed_pb = .TRUE.
        END IF

    END SUBROUTINE update_coordinates

END MODULE particle_boundaries_mod