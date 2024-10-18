MODULE particle_boundaries_mod

    USE precision_mod, ONLY: realk, intk
    USE core_mod

    USE particle_list_mod
    USE particle_obstacles_mod
    USE particle_utils_mod

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

    TYPE :: particle_boundaries_t

        INTEGER(intk), ALLOCATABLE :: face_neighbours(:, :)

        REAL(realk), ALLOCATABLE :: face_normals(:, :, :)

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

        CALL start_timer(900)
        CALL start_timer(910)

        ALLOCATE(particle_boundaries%face_normals(3, 26, ngrid))
        particle_boundaries%face_normals = 0.0

        ALLOCATE(particle_boundaries%face_neighbours(26, ngrid))

        ibocd = 2
        DO igrid = 1, ngrid

            CALL get_neighbours(neighbours, igrid)

            DO iface = 1, 6

                CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)

                particle_boundaries%face_neighbours(iface, igrid) = neighbours(iface)

                SELECT CASE(iface)
                    CASE(1)
                        IF (ctyp == 'SWA') THEN
                            particle_boundaries%face_normals(1, 1, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(2)
                        IF (ctyp == 'SWA') THEN
                            particle_boundaries%face_normals(1, 2, igrid) = -1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(3)
                        IF (ctyp == 'SWA') THEN
                            particle_boundaries%face_normals(2, 3, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(4)
                        IF (ctyp == 'SWA') THEN
                            particle_boundaries%face_normals(2, 4, igrid) = -1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(5)
                        IF (ctyp == 'SWA') THEN
                            particle_boundaries%face_normals(3, 5, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(6)
                        IF (ctyp == 'SWA') THEN
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

                        IF (ctyp == "SIO") THEN
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
                    IF (magnitude == 0.0) THEN
                        particle_boundaries%face_normals(j, iface, igrid) = 0.0
                    ELSE
                        particle_boundaries%face_normals(j, iface, igrid) = particle_boundaries%face_normals(j, iface, igrid) / magnitude
                    END IF
                END DO

            END DO
        END DO

        DO igrid = 1, ngrid

            WRITE(*, *) "------ Boundaries, Grid:   ", igrid, "------"
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

        WRITE(*, *) " "

        CALL read_obstacles()

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_boundaries

    !-----------------------------------

    SUBROUTINE finish_particle_boundaries()

        CALL finish_obstacles()

        DEALLOCATE(particle_boundaries%face_neighbours)
        DEALLOCATE(particle_boundaries%face_normals)

    END SUBROUTINE finish_particle_boundaries

    !-----------------------------------

    SUBROUTINE move_particle(particle, dx, dy, dz, dx_eff, dy_eff, dz_eff)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        REAL(realk), INTENT(in) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_eff, dy_eff, dz_eff

        ! local variables
        INTEGER(intk) :: temp_grid, iface, iobst, grid_bc, destproc
        INTEGER(intk) :: counter
        INTEGER(intk) :: neighbours(26)
        REAL(realk) :: x, y, z
        REAL(realk) :: dx_step, dy_step, dz_step
        REAL(realk) :: dx_from_here, dy_from_here, dz_from_here

        REAL(realk) :: epsilon
        REAL(realk) :: n1, n2, n3

        CALL start_timer(924)

        IF (SQRT(dx**(2) + dy**(2) + dz**(2)) == 0) THEN
            CALL stop_timer(924)
            RETURN
        END IF

        epsilon = SQRT(dx**(2) + dy**(2) + dz**(2)) / 10**3

        temp_grid = particle%igrid

        x = particle%x
        y = particle%y
        z = particle%z

        dx_from_here = dx
        dy_from_here = dy
        dz_from_here = dz

        dx_eff = 0.0
        dy_eff = 0.0
        dz_eff = 0.0

        iobst = 0

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                WRITE(*, *) "MOVE PARTICLE: ----------------------START-----------------------"
                WRITE(*, '()')
        END SELECT

        counter = 1
        DO WHILE (epsilon < SQRT(dx_from_here**(2) + dy_from_here**(2) + dz_from_here**(2)) .AND. counter <= 10)

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

            CALL move_to_boundary(temp_grid, x, y, z, &
             dx_from_here, dy_from_here, dz_from_here, dx_step, dy_step, dz_step, iface, iobst)

            dx_eff = dx_eff + dx_step
            dy_eff = dy_eff + dy_step
            dz_eff = dz_eff + dz_step

            IF (0 < iobst) THEN

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*, *) "Particle reflected at Obstacle!"
                END SELECT

                CALL reflect_at_obstacle(iobst, x, y, z, dx_from_here, dy_from_here, dz_from_here)

            ELSEIF (0 < iface) THEN

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*, *) "Particle reflected at Grid face."
                END SELECT

                CALL reflect_at_boundary(dx_from_here, dy_from_here, dz_from_here, &
                particle_boundaries%face_normals(1, iface, temp_grid), &
                particle_boundaries%face_normals(2, iface, temp_grid), &
                particle_boundaries%face_normals(3, iface, temp_grid))

                CALL update_coordinates(temp_grid, particle_boundaries%face_neighbours(iface, temp_grid), iface, x, y, z)

                temp_grid = particle_boundaries%face_neighbours(iface, temp_grid)

            END IF

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

            counter = counter + 1

        END DO

        ! do not update the particle grid here!
        particle%x = particle%x + dx_eff
        particle%y = particle%y + dy_eff
        particle%z = particle%z + dz_eff

        CALL update_particle_cell(particle)

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                WRITE(*, *) "MOVE PARTICLE: -----------------------END------------------------"
                WRITE(*, '()')
        END SELECT

        CALL stop_timer(924)

    END SUBROUTINE move_particle

    !-----------------------------------

    ! This subroutine only considers grids on the same level
    ! CAUTION: Here, temp_grid refers to the grid the particle coordinates are currently on and of which the boundaries are relevant.
    ! This might NOT be particle%igrid, which is used to deduce the velocity from.
    SUBROUTINE move_to_boundary(temp_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, iface, iobst)

        ! subroutine arguments
        INTEGER(intk), INTENT(inout) :: temp_grid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        INTEGER(intk), INTENT(out) :: iface
        INTEGER(intk), INTENT(inout) :: iobst

        !local variables
        INTEGER(intk) :: i, nobst
        REAL(realk) :: dist
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: lx, ly, lz, rx, ry, rz
        REAL(realk) :: s, sa, sb, smin, a, b, c, d, cx, cy, cz, r

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, temp_grid)

        ! STEP 1 - OBSTACLES
        ! find intersection points of the line the particle moves on (straight) and the sphere surface
            ! particle path: X(s) = X + dX * s with s: [0, 1] (X is the vector (x/y/z))
            ! => |X + dX * s - C| =! r (C is the sphere center (cx/cy/cz))
            ! => (x + dx * s -cx)² + (y + dy * s -cy)² + (z + dz * s -cz)² =! r² (r is the sphere radius)
            ! => s1/s2 = sa/sb = (-b +/- sqrt(b² - 4ac)) / 2a (corefficients see code)

        smin = 1.0
        IF (0 < ABS(dx) .AND. ABS(dx) < ABS(dy) .AND. ABS(dx) < ABS(dz)) THEN
            smin =  MIN(smin, ABS(EPSILON(s) / dx))
        ELSEIF (0 < ABS(dy) .AND. ABS(dy) < ABS(dx) .AND. ABS(dy) < ABS(dz)) THEN
            smin =  MIN(smin, ABS(EPSILON(s) / dy))
        ELSEIF (0 < ABS(dz) .AND. ABS(dz) < ABS(dx) .AND. ABS(dz) < ABS(dy)) THEN
            smin =  MIN(smin, ABS(EPSILON(s) / dz))
        END IF

        s = 1.0

        ! first coefficient
        a = (dx**2 + dy**2 + dz**2)

        ! iterate over all obstacles of the grid
        nobst = SIZE(obstacles)

        DO i = 1, nobst

            ! dont check if a particle interacts with the obstacle it has been deflected from in the previous timestep
            IF (i == iobst) THEN
                iobst = 0
                CYCLE
            END IF

            ! for readability
            cx = obstacles(i)%x
            cy = obstacles(i)%y
            cz = obstacles(i)%z
            r = obstacles(i)%radius

            IF (SQRT((x - cx)**2 + (y - cy)**2 + (z - cz)**2) < r) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING: In move_to_boundary: Particle inside Obstacle."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING: In move_to_boundary: Particle Inside Obstacle."
                        WRITE(*, '()')
                END SELECT
            END IF

            ! sphere dependent coefficients
            b = 2*x*dx + 2*y*dy + 2*z*dz - 2*cx*dx - 2*cy*dy - 2*cz*dz
            c = x**2 + y**2 + z**2 + cx**2 + cy**2 + cz**2 - 2*x*cx - 2*y*cy - 2*z*cz - r**2

            d = b**2 - 4*a*c

            IF (d < 0) THEN
                CONTINUE
            END IF

            sa = (-b + SQRT(d)) / 2 / a
            sb = (-b - SQRT(d)) / 2 / a

            ! if sa is outside [0, 1], the intersection point is not on the actual limited path of the particel
            ! => set sa to an arbitrary number higher than 1, so it wont get registered later
            IF (sa < 0 .OR. 1 < sa) THEN
                sa = 1.1
            ! if sa is 0, but the particle moves away from the boundary it is on
            ! => set sa to an arbitrary number higher than 1, so it wont get registered later
            ELSEIF (sa == 0 .AND. r < SQRT((x + smin * dx - cx)**2 + (y + smin * dy - cy)**2 + (z + smin * dz - cz)**2)) THEN
                sa = 1.1
            END IF

            ! if sb is outside [0, 1], the intersection point is not on the actual limited path of the particel
            ! => set sb to an arbitrary number higher than 1, so it wont get registered later
            IF (sb < 0 .OR. 1 < sb) THEN
                sb = 1.1
            ! if sb is 0, but the particle moves away from the boundary it is on
            ! => set sb to an arbitrary number higher than 1, so it wont get registered later
            ELSEIF (sb == 0 .AND. r < SQRT((x + smin * dx - cx)**2 + (y + smin * dy - cy)**2 + (z + smin * dz - cz)**2)) THEN
                sb = 1.1
            END IF

            IF (sa < s) THEN
                s = sa
                iobst = i
            END IF

            IF (sb < s) THEN
                s = sb
                iobst = i
            END IF

        END DO

        ! STEP 2 - GRID BOUNDARIES
        ! now check if any grid boundary is reached before any obstacle is reached

        IF (dx < 0) THEN
            lx = (minx - x)
            IF(lx == 0) THEN
                CALL get_exit_face(temp_grid, x, y, z, dist, iface)
                RETURN
            END IF
            rx = dx * s / lx
        ELSEIF (0 < dx) THEN
            lx = (maxx - x)
            IF(lx == 0) THEN
                CALL get_exit_face(temp_grid, x, y, z, dist, iface)
                RETURN
            END IF
            rx = dx * s / lx
        ELSE
            rx = 0.0_realk
        END IF

        IF (dy < 0) THEN
            ly = (miny - y)
            IF(ly == 0) THEN
                CALL get_exit_face(temp_grid, x, y, z, dist, iface)
                RETURN
            END IF
            ry = dy * s / ly
        ELSEIF (0 < dy) THEN
            ly = (maxy - y)
            IF(ly == 0) THEN
                CALL get_exit_face(temp_grid, x, y, z, dist, iface)
                RETURN
            END IF
            ry = dy * s / ly
        ELSE
            ry = 0.0_realk
        END IF

        IF (dz < 0) THEN
            lz = (minz - z)
            IF(lz == 0) THEN
                CALL get_exit_face(temp_grid, x, y, z, dist, iface)
                RETURN
            END IF
            rz = dz * s / lz
        ELSEIF (0 < dz) THEN
            lz = (maxz - z)
            IF(lz == 0) THEN
                CALL get_exit_face(temp_grid, x, y, z, dist, iface)
                RETURN
            END IF
            rz = dz * s / lz
        ELSE
            rz = 0.0_realk
        END IF

        IF (rx < 1.0_realk .AND. ry < 1.0_realk .AND. rz < 1.0_realk) THEN

            dx_to_b = dx * s
            dy_to_b = dy * s
            dz_to_b = dz * s
            x = x + dx_to_b
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

            iface = 0

            RETURN

        END IF

        ! if the routine did not return yet, no obstacle will be hit before some grid boundary
        iobst = 0

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
        CALL get_exit_face(temp_grid, x, y, z, dist, iface)

    END SUBROUTINE move_to_boundary

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

        ! Presumption: Particle is already exactly on the boundary!

        ! subroutine arguments
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(in) :: n1, n2, n3 ! normal vector components of the surface the particle is reflected from

        ! local variables
        REAL(realk) :: dot_product

        dot_product = n1 * dx + n2 * dy + n3 * dz

        IF (dot_product < 0) THEN
            dx = dx - 2 * dot_product * n1
            dy = dy - 2 * dot_product * n2
            dz = dz - 2 * dot_product * n3
        ELSE
            dx = dx
            dy = dy
            dz = dz
        END IF

    END SUBROUTINE reflect_at_boundary

    ! TODO: make this (partly) an obstacle method
    SUBROUTINE reflect_at_obstacle(iobst, x, y, z, dx, dy, dz)

        ! Presumption 1: Particle is already exactly on the boundary!
        ! Presumption 2: Obstacle is a sphere!

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: iobst
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz

        ! local variables
        REAL(realk) :: n1, n2 , n3, magnitude

        n1 = x - obstacles(iobst)%x
        n2 = y - obstacles(iobst)%y
        n3 = z - obstacles(iobst)%z

        magnitude = SQRT(n1**2 + n2**2 + n3**2)

        n1 = n1 / magnitude
        n2 = n2 / magnitude
        n3 = n3 / magnitude

        CALL reflect_at_boundary(dx, dy, dz, n1, n2, n3)

    END SUBROUTINE reflect_at_obstacle

END MODULE particle_boundaries_mod