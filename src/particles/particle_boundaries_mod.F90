MODULE particle_boundaries_mod

    USE precision_mod, ONLY: realk, intk
    USE core_mod !TODO: specify
    USE utils_mod

    USE particle_obstacles_mod
    USE particle_runtimestat_mod, ONLY: psim_n_replaced_tot, &
     psim_n_bcerr, psim_max_bcerr

    IMPLICIT NONE

    INTEGER(intk), PARAMETER :: facelist_b(4,26) = RESHAPE((/ &
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
        3, 2, 4, 6 /), SHAPE(facelist_b))

    ! TODO: restructure type so that particle_boundaries is an array if size nmygrids
    ! and each attribute has one dimension (that now holds igrid info) less
    TYPE :: particle_boundaries_t

        INTEGER(intk), ALLOCATABLE :: face_neighbours(:, :)

        REAL(realk), ALLOCATABLE :: face_normals(:, :, :)

    END TYPE particle_boundaries_t

    TYPE(particle_boundaries_t) :: particle_boundaries

    CHARACTER(len = 4) :: bc_coupling_mode = "SCAL" ! must be "FLOW", "SCAL" or "PART"

    CONTAINS

    SUBROUTINE init_particle_boundaries()

        ! local variables
        INTEGER(intk) :: igrid, iface, jface, i, j
        INTEGER(intk) :: neighbours(26)
        ! to store which boundary surfaces (faces 1-6) that define a lower order face (7-26) are of connect (or periodic) type:
        INTEGER(intk) :: connect_faces(4)
        REAL(realk) :: magnitude
        CHARACTER(len=3) :: ctyp
        LOGICAL :: found

        CALL start_timer(900)
        CALL start_timer(910)

        ALLOCATE(particle_boundaries%face_normals(3, 26, ngrid))
        particle_boundaries%face_normals = 0.0

        ALLOCATE(particle_boundaries%face_neighbours(26, ngrid))

        ! CON : connective particle boundary (can be periodic)
        ! REF : reflective particle boundary
        DO igrid = 1, ngrid

            ! TODO: possibly cycle if grid is not on particle_level
            !IF (level(igrid) /= particle_level) CYCLE

            CALL get_neighbours(neighbours, igrid)

            DO iface = 1, 6

                CALL get_particle_bc(igrid, iface, bc_coupling_mode, ctyp)

                particle_boundaries%face_neighbours(iface, igrid) = neighbours(iface)

                SELECT CASE(iface)
                    CASE(1)
                        IF (ctyp == "REF") THEN ! ctyp == "SWA"
                            particle_boundaries%face_normals(1, 1, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(2)
                        IF (ctyp == "REF") THEN
                            particle_boundaries%face_normals(1, 2, igrid) = -1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(3)
                        IF (ctyp == "REF") THEN
                            particle_boundaries%face_normals(2, 3, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(4)
                        IF (ctyp == "REF") THEN
                            particle_boundaries%face_normals(2, 4, igrid) = -1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(5)
                        IF (ctyp == "REF") THEN
                            particle_boundaries%face_normals(3, 5, igrid) = 1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                    CASE(6)
                        IF (ctyp == "REF") THEN
                            particle_boundaries%face_normals(3, 6, igrid) = -1.0
                            particle_boundaries%face_neighbours(iface, igrid) = igrid
                        END IF
                END SELECT

            END DO

            DO iface = 7, 26

                connect_faces = 0

                DO i = 2, 4
                    IF (facelist_b(i, iface) == 0) THEN
                        CONTINUE
                    ELSE

                        CALL get_particle_bc(igrid, facelist_b(i, iface), bc_coupling_mode, ctyp)

                        IF (ctyp == "CON") THEN ! ctyp == "SIO"
                            connect_faces(1) = connect_faces(1) + 1
                            connect_faces(i) = facelist_b(i, iface)
                        END IF

                        DO j = 1, 3
                            particle_boundaries%face_normals(j, iface, igrid) = &
                            particle_boundaries%face_normals(j, iface, igrid) + particle_boundaries%face_normals(j, facelist_b(i, iface), igrid)
                        END DO

                    END IF
                END DO

                particle_boundaries%face_neighbours(iface, igrid) = igrid
                found = .FALSE.
                DO jface = 1, 26
                    IF (facelist_b(1, jface) /= connect_faces(1)) THEN
                        CONTINUE
                    ELSE
                        found = .TRUE.
                        DO i = 2, 1 + facelist_b(1, jface)
                            IF (facelist_b(i, jface) /= connect_faces(2) &
                             .AND. facelist_b(i, jface) /= connect_faces(3) .AND. facelist_b(i, jface) /= connect_faces(4)) THEN
                                found = .FALSE.
                                EXIT
                            ELSE
                                CONTINUE
                            END IF
                        END DO
                        IF (found .eqv. .TRUE.) THEN
                            particle_boundaries%face_neighbours(iface, igrid) = neighbours(jface)
                        END IF
                    END IF
                END DO

                magnitude = SQRT(particle_boundaries%face_normals(1, iface, igrid)**2 + &
                    particle_boundaries%face_normals(2, iface, igrid)**2 + &
                    particle_boundaries%face_normals(3, iface, igrid)**2)

                DO j = 1, 3
                    ! EPSILON(magnitude) is an arbitrary value significantely smaller than 1 as
                    ! the shortest valid normal vector up to now should have a magnitude of 1
                    IF (magnitude <= EPSILON(magnitude)) THEN
                        particle_boundaries%face_normals(j, iface, igrid) = 0.0
                    ELSE
                        particle_boundaries%face_normals(j, iface, igrid) = particle_boundaries%face_normals(j, iface, igrid) / magnitude
                    END IF
                END DO

            END DO
        END DO

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "verbose") THEN
                DO igrid = 1, ngrid
                    WRITE(*, *) "------ Boundaries, Grid:   ", igrid, "------"
                    CALL get_particle_bc(igrid, 1, bc_coupling_mode, ctyp)
                    WRITE(*, *) "FRONT:                ", ctyp
                    CALL get_particle_bc(igrid, 2, bc_coupling_mode, ctyp)
                    WRITE(*, *) "BACK:                 ", ctyp
                    CALL get_particle_bc(igrid, 3, bc_coupling_mode, ctyp)
                    WRITE(*, *) "RIGHT:                ", ctyp
                    CALL get_particle_bc(igrid, 4, bc_coupling_mode, ctyp)
                    WRITE(*, *) "LEFT:                 ", ctyp
                    CALL get_particle_bc(igrid, 5, bc_coupling_mode, ctyp)
                    WRITE(*, *) "BOTTOM:               ", ctyp
                    CALL get_particle_bc(igrid, 6, bc_coupling_mode, ctyp)
                    WRITE(*, *) "TOP:                  ", ctyp
                    WRITE(*, *) " "
                    WRITE(*, *) "FACES:"
                    DO iface = 1, 26
                        WRITE(*, *) "Face:                 ", iface
                        WRITE(*, *) "Neigbhour grid:       ", particle_boundaries%face_neighbours(iface, igrid)
                        WRITE(*,'("Normal vector: ")')
                        WRITE(*, *) "n1                    ", particle_boundaries%face_normals(1, iface, igrid)
                        WRITE(*, *) "n2                    ", particle_boundaries%face_normals(2, iface, igrid)
                        WRITE(*, *) "n3                    ", particle_boundaries%face_normals(3, iface, igrid)
                        WRITE(*, *) " "
                    END DO
                END DO
                WRITE(*, '()')
            END IF
        END IF

        ! TODO: remove?
        CALL MPI_Barrier(MPI_COMM_WORLD)

        CALL read_obstacles()

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_boundaries

    !-----------------------------------

    SUBROUTINE finish_particle_boundaries()

        CALL start_timer(900)
        CALL start_timer(910)

        CALL finish_obstacles()

        IF (ALLOCATED(particle_boundaries%face_neighbours)) DEALLOCATE(particle_boundaries%face_neighbours)
        IF (ALLOCATED(particle_boundaries%face_normals)) DEALLOCATE(particle_boundaries%face_normals)

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE finish_particle_boundaries

    !-----------------------------------

    ! TODO: source the following moduled out into particle_motion_mod or so

    SUBROUTINE move_particle(particle, dx, dy, dz, dx_eff, dy_eff, dz_eff, temp_grid_prev)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        REAL(realk), INTENT(in) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_eff, dy_eff, dz_eff
        INTEGER(intk), INTENT(inout), OPTIONAL :: temp_grid_prev

        ! local variables
        INTEGER(intk) :: temp_grid, iface, iobst_local, destgrid, counter
        INTEGER(intk) :: reflect(3)
        REAL(realk) :: x, y, z
        REAL(realk) :: dx_step, dy_step, dz_step
        REAL(realk) :: dx_from_here, dy_from_here, dz_from_here
        REAL(realk) :: eps
        LOGICAL :: dreplace

        dreplace = .FALSE.

        dx_eff = 0.0
        dy_eff = 0.0
        dz_eff = 0.0

        IF (SQRT(dx**(2) + dy**(2) + dz**(2)) <= EPSILON(0.0_realk) ) THEN
            RETURN
        END IF

        ! TODO: make this a reasonable stoping criterion / rethink
        eps = SQRT(dx**(2) + dy**(2) + dz**(2)) / 10.0_realk**3

        IF (PRESENT(temp_grid_prev)) THEN
            temp_grid = temp_grid_prev
        ELSE
            temp_grid = particle%igrid
        END IF

        x = particle%x
        y = particle%y
        z = particle%z

        dx_from_here = dx
        dy_from_here = dy
        dz_from_here = dz

        iobst_local = 0

        IF (TRIM(particle_terminal) == "verbose") THEN
            WRITE(*, *) "---------- MOVE PARTICLE ----------"
            WRITE(*, '()')
        END IF

        counter = 1
        DO WHILE (SQRT(dx_from_here**(2) + dy_from_here**(2) + dz_from_here**(2)) > eps .AND. &
         MAX(ABS(dx_from_here), ABS(dy_from_here), ABS(dz_from_here)) > EPSILON(0.0_realk) .AND. counter <= 10)

            IF (TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, *) "Way to go:"
                WRITE(*, *) "dx/dy/dz:", dx_from_here, dy_from_here, dz_from_here
            END IF

            CALL move_to_boundary(particle%ipart, temp_grid, x, y, z, &
             dx_from_here, dy_from_here, dz_from_here, dx_step, dy_step, dz_step, iface, iobst_local, dreplace)

            ! replace current particle coordinates by a random valid position on the particles curren grid
            IF (dreplace) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, '("WARNING: In move_particle:")')
                    WRITE(*, *) "Particle", particle%ipart, " to be replaced!"
                    WRITE(*, '()')
                END IF
                CALL replace_particle(particle)
                temp_grid = particle%igrid
                x = particle%x
                y = particle%y
                z = particle%z
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    psim_n_replaced_tot = psim_n_replaced_tot + 1
                END IF
                dreplace = .FALSE.
            END IF

            dx_eff = dx_eff + dx_step
            dy_eff = dy_eff + dy_step
            dz_eff = dz_eff + dz_step

            IF (0 < iobst_local) THEN

                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "Particle reflected at Obstacle ", my_obstacles(iobst_local)%iobst, "."
                END IF

                CALL reflect_at_obstacle(iobst_local, x, y, z, dx_from_here, dy_from_here, dz_from_here)

            ELSEIF (0 < iface) THEN

                destgrid = particle_boundaries%face_neighbours(iface, temp_grid)

                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "Particle moved to Grid face ", iface, "with Target Grid ", destgrid, "."
                END IF

                CALL reflect_at_boundary(dx_from_here, dy_from_here, dz_from_here, &
                 particle_boundaries%face_normals(1, iface, temp_grid), &
                 particle_boundaries%face_normals(2, iface, temp_grid), &
                 particle_boundaries%face_normals(3, iface, temp_grid), reflect)

                ! TODO (LONGTERM): restructure boundaries and particle motion such that particle coordinates
                ! do not have to be updated during substeps! (performance)
                CALL update_coordinates(temp_grid, destgrid, iface, x, y, z, reflect)

                temp_grid = destgrid

                ! PROBABLY NOT NEEDED DUE TO ADJUSTMENTS IN MOVE_TO_BOUNDARY
                ! Move particle sligthly to avoid particles getting stuck at a boundary.
                ! The way particle motion works particles would otherwise get stuck on edges and corners,
                ! if they are right on an edge or corner at the beginning of move_particle.
                ! If the following small displacement leads to a particle being outside temp_grid or inside an obstacle,
                ! no problems should occur as the algorithm should be stable in that reguard.
                !x = x + SIGN(EPSILON(x), dx_from_here)
                !y = y + SIGN(EPSILON(y), dy_from_here)
                !z = z + SIGN(EPSILON(z), dz_from_here)

            END IF

            IF (TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, *) "Intermediate Location:"
                WRITE(*, *) "igrid:", temp_grid
                WRITE(*, *) "x:", x
                WRITE(*, *) "y:", y
                WRITE(*, *) "z:", z
                WRITE(*, '()')
            END IF

            counter = counter + 1

        END DO

        ! do not update the particle grid here!
        particle%x = particle%x + dx_eff
        particle%y = particle%y + dy_eff
        particle%z = particle%z + dz_eff

        temp_grid_prev = temp_grid

        CALL update_particle_cell(particle)

    END SUBROUTINE move_particle

    !-----------------------------------

    ! This subroutine only considers grids on the same level
    ! CAUTION: Here, temp_grid refers to the grid the particle coordinates are currently on and of which the boundaries are relevant.
    ! This might NOT be particle%igrid, which is used to deduce the velocity
    SUBROUTINE move_to_boundary(ipart, temp_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, iface, iobst_local, replace)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: ipart
        INTEGER(intk), INTENT(inout) :: temp_grid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        INTEGER(intk), INTENT(out) :: iface
        INTEGER(intk), INTENT(inout) :: iobst_local
        LOGICAL, INTENT(out) :: replace

        !local variables
        INTEGER(intk) :: i, nobst
        REAL(realk) :: dist, dist_to_go, dist_to_center
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: lx, ly, lz, rx, ry, rz
        REAL(realk) :: s, sa, sb, sc, sd, a, b, c, d, cx, cy, cz, r

        replace = .FALSE.

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, temp_grid)

        dx_to_b = 0.0
        dy_to_b = 0.0
        dz_to_b = 0.0

        ! STEP 1 - OBSTACLES
        ! find intersection points of the line the particle moves on (straight) and the sphere surface
            ! particle path: X(s) = X + dX * s with s: [0, 1] (X is the vector (x/y/z))
            ! => |X + dX * s - C| = r (C is the sphere center (cx/cy/cz))
            ! => (x + dx * s -cx)² + (y + dy * s -cy)² + (z + dz * s -cz)² = r² (r is the sphere radius)
            ! => s1/s2 = sa/sb = (-b +/- sqrt(b² - 4ac)) / 2a (corefficients see code)

        dist_to_go = SQRT((dx)**2 + (dy)**2 + (dz)**2)

        s = 1.0

        ! first coefficient
        a = (dx**2 + dy**2 + dz**2)

        ! iterate over all obstacles of the grid
        IF (dread_obstacles_dict) THEN
            nobst = SIZE(my_obstacle_pointers(temp_grid)%grid_obstacles)
        ELSE
            nobst = 0
        END IF

        DO i = 1, nobst

            ! check if a particle interacts with the obstacle it has been deflected from in the previous timestep
            IF (my_obstacle_pointers(temp_grid)%grid_obstacles(i) == iobst_local) THEN
                CYCLE
            END IF

            ! for readability
            cx = my_obstacles(my_obstacle_pointers(temp_grid)%grid_obstacles(i))%x
            cy = my_obstacles(my_obstacle_pointers(temp_grid)%grid_obstacles(i))%y
            cz = my_obstacles(my_obstacle_pointers(temp_grid)%grid_obstacles(i))%z
            r = my_obstacles(my_obstacle_pointers(temp_grid)%grid_obstacles(i))%radius

            ! sphere dependent coefficients
            b = 2*x*dx + 2*y*dy + 2*z*dz - 2*cx*dx - 2*cy*dy - 2*cz*dz
            c = x**2 + y**2 + z**2 + cx**2 + cy**2 + cz**2 - 2*x*cx - 2*y*cy - 2*z*cz - r**2
            d = b**2 - 4*a*c

            IF (d < EPSILON(0.0_realk)) THEN
                CYCLE
            END IF

            IF (a < EPSILON(0.0_realk)) THEN
                CYCLE
            END IF

            sa = (-b + SQRT(d)) / 2 / a
            sb = (-b - SQRT(d)) / 2 / a

            ! if a particle moves towards an obstacle, limit its motion to the closest intersection yet
            IF (sa >= 0.0 .AND. sb >= 0.0) THEN
                sc = MIN(sa, sb)
                IF (sc < s) THEN
                    s = sc
                    iobst_local = my_obstacle_pointers(temp_grid)%grid_obstacles(i)
                END IF
            ! elseif a particle moves away from the current obstacle, cycle
            ELSEIF (sa <= 0.0 .AND. sb <= 0.0) THEN
                CYCLE
            ! else (if sa < 0 and sb > 0 or vice versa) the particle is inside the current obstacle
            ! => replace current particle coordinates by a random valid position on the particles curren grid
            ELSE
                dist_to_center = SQRT((x - cx)**2 + (y - cy)**2 + (z - cz)**2)

                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, '("WARNING: In move_to_boundary:")')
                    WRITE(*, *) "Particle", ipart," inside Obstacle by ", (r - dist_to_center)
                    WRITE(*, '()')
                END IF

                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                    psim_max_bcerr = MAX(psim_max_bcerr, (r - dist_to_center))
                    psim_n_bcerr = psim_n_bcerr + 1
                END IF

                IF ((r - dist_to_center) > aura) THEN
                    IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                        WRITE(*, '("WARNING: In move_to_boundary:")')
                        WRITE(*, *) "Particle", ipart," inside Obstacle by ", (r - dist_to_center)
                        WRITE(*, '()')
                    END IF
                    replace = .TRUE.
                    iface = 0
                    iobst_local = 0
                    RETURN
                END IF

                sc = MIN(sa, sb)
                sd = MAX(sa, sb)

                IF (ABS(sc) < ABS(sd)) THEN
                    s = 0.0
                    iobst_local = my_obstacle_pointers(temp_grid)%grid_obstacles(i)
                    EXIT
                ELSEIF (ABS(sc) >= ABS(sd)) THEN
                    CYCLE
                END IF

            END IF
        END DO

        ! TODO: return if s = 0
        IF (s <= 0.0_realk) THEN
            iface = 0
            RETURN
        END IF

        ! STEP 2 - GRID BOUNDARIES
        ! now check if any grid boundary is reached before any obstacle is reached
        IF (dx < 0) THEN
            lx = (minx - x)
            ! if particle is at boundary in X dir (lx = 0) or particle is outside temp_grid (lx > 0.0),
            ! get exit face and return; so if a particle is incorrectly outside a reflect boundary its
            ! motion vector is reflected towards temp_grid
            IF (lx >= 0.0_realk) THEN
                iobst_local = 0
                ! if a particle is already on a face (esp. edge or corner), its future coordinates have to be
                ! "projected" to assign the right ecit face (otherwise, particles might get stuck on edges or corners)
                CALL get_exit_face(temp_grid, x + dx, y + dy, z + dz, dist, iface)
                RETURN
            END IF
            rx = dx * s / lx
        ELSEIF (0 < dx) THEN
            lx = (maxx - x)
            IF (lx <= 0.0_realk) THEN
                iobst_local = 0
                CALL get_exit_face(temp_grid, x + dx, y + dy, z + dz, dist, iface)
                RETURN
            END IF
            rx = dx * s / lx
        ELSE
            rx = 0.0_realk
        END IF

        IF (dy < 0) THEN
            ly = (miny - y)
            IF (ly >= 0.0_realk) THEN
                iobst_local = 0
                CALL get_exit_face(temp_grid, x + dx, y + dy, z + dz, dist, iface)
                RETURN
            END IF
            ry = dy * s / ly
        ELSEIF (0 < dy) THEN
            ly = (maxy - y)
            IF (ly <= 0.0_realk) THEN
                iobst_local = 0
                CALL get_exit_face(temp_grid, x + dx, y + dy, z + dz, dist, iface)
                RETURN
            END IF
            ry = dy * s / ly
        ELSE
            ry = 0.0_realk
        END IF

        IF (dz < 0) THEN
            lz = (minz - z)
            IF(lz >= 0.0_realk) THEN
                iobst_local = 0
                CALL get_exit_face(temp_grid, x + dx, y + dy, z + dz, dist, iface)
                RETURN
            END IF
            rz = dz * s / lz
        ELSEIF (0 < dz) THEN
            lz = (maxz - z)
            IF(lz <= 0.0_realk) THEN
                iobst_local = 0
                CALL get_exit_face(temp_grid, x + dx, y + dy, z + dz, dist, iface)
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
        iobst_local = 0

        IF (dx < 0 .AND. ry <= rx .AND. rz <= rx) THEN

            dx_to_b = lx
            dy_to_b = (lx * dy/dx)
            dz_to_b = (lx * dz/dx)
            x = minx ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dx .AND. ry <= rx .AND. rz <= rx) THEN

            dx_to_b = lx
            dy_to_b = (lx * dy/dx)
            dz_to_b = (lx * dz/dx)
            x = maxx ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
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
            y = miny ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dy .AND. rx < ry .AND. rz <= ry) THEN

            dx_to_b = (ly * dx/dy)
            dy_to_b = ly
            dz_to_b = (ly * dz/dy)
            x = x + dx_to_b
            y = maxy ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
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
            z = minz ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dz .AND. rx < rz .AND. ry < rz) THEN

            dx_to_b = (lz * dx/dz)
            dy_to_b = (lz * dy/dz)
            dz_to_b = lz
            x = x + dx_to_b
            y = y + dy_to_b
            z = maxz ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        END IF

        ! get face after x/y/z have (potentially) been altered
        CALL get_exit_face(temp_grid, x, y, z, dist, iface)

    END SUBROUTINE move_to_boundary

    SUBROUTINE replace_particle(particle)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle

        ! local variables
        LOGICAL :: valid_location
        INTEGER(intk) :: i, iobst_local, igrid
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, x_new, y_new, z_new, dist_to_center

        ! for readability
        igrid = particle%igrid

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        valid_location = .FALSE.
        DO WHILE (.NOT. valid_location)

            valid_location = .TRUE.

            CALL RANDOM_NUMBER(x_new)
            CALL RANDOM_NUMBER(y_new)
            CALL RANDOM_NUMBER(z_new)

            x_new = minx + x_new * (maxx - minx)
            y_new = miny + y_new * (maxy - miny)
            z_new = minz + z_new * (maxz - minz)

            IF (dread_obstacles_dict) THEN
                DO i = 1, SIZE(my_obstacle_pointers(igrid)%grid_obstacles)

                    iobst_local = my_obstacle_pointers(igrid)%grid_obstacles(i)

                    dist_to_center = SQRT((my_obstacles(iobst_local)%x - x_new)**2 + &
                     (my_obstacles(iobst_local)%y - y_new)**2 + &
                     (my_obstacles(iobst_local)%z - z_new)**2)

                    IF (dist_to_center < my_obstacles(iobst_local)%radius + EPSILON(dist_to_center)) THEN
                        valid_location = .FALSE.
                        EXIT
                    ELSE
                        CONTINUE
                    END IF

                END DO
            END IF

        END DO

        particle%x = x_new
        particle%y = y_new
        particle%z = z_new

        CALL set_particle_cell(particle)

    END SUBROUTINE replace_particle

    ! out of use
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

    SUBROUTINE reflect_at_boundary(dx, dy, dz, n1, n2, n3, reflect)

        ! Presumption: Particle is already exactly on the boundary!

        ! subroutine arguments
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(in) :: n1, n2, n3 ! normal vector components of the surface the particle is reflected from
        INTEGER(intk), INTENT(out), OPTIONAL :: reflect(3)

        ! local variables
        REAL(realk) :: dot_product

        IF (PRESENT(reflect)) THEN
            reflect = 0
            IF (0 < ABS(n1)) reflect(1) = 1
            IF (0 < ABS(n2)) reflect(2) = 1
            IF (0 < ABS(n3)) reflect(3) = 1
        END IF

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
    SUBROUTINE reflect_at_obstacle(iobst_local, x, y, z, dx, dy, dz)

        ! Presumption 1: Particle is already exactly on the boundary!
        ! Presumption 2: Obstacle is a sphere!

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: iobst_local
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz

        ! local variables
        REAL(realk) :: n1, n2 , n3, magnitude

        n1 = x - my_obstacles(iobst_local)%x
        n2 = y - my_obstacles(iobst_local)%y
        n3 = z - my_obstacles(iobst_local)%z

        magnitude = SQRT(n1**2 + n2**2 + n3**2)

        n1 = n1 / magnitude
        n2 = n2 / magnitude
        n3 = n3 / magnitude

        CALL reflect_at_boundary(dx, dy, dz, n1, n2, n3)

    END SUBROUTINE reflect_at_obstacle

    SUBROUTINE get_particle_bc(igrid, iface, coupling_mode, ctyp)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: iface
        CHARACTER(len = 4), INTENT(in) :: coupling_mode
        CHARACTER(len = 3), INTENT(out) :: ctyp

        ! local variables
        INTEGER(intk) :: ibocd, nbocd
        CHARACTER(len = 8) :: ctyp1, ctyp2

        ! TODO: rework particle boundary initialization and storage
        ! SIO = Skalar-RB für Oberflächen der Domain, die durchströmt werden
        ! SWA = Skalar-RB auf Wänden (slip und no-slip)

        SELECT CASE(coupling_mode)
        CASE("PART")

            ibocd = 2
            CALL get_bc_ctyp(ctyp2, ibocd, iface, igrid)

        CASE("SCAL")

            nbocd = nboconds(iface, igrid)
            ibocd = nbocd
            CALL get_bc_ctyp(ctyp2, ibocd, iface, igrid)

            IF (ctyp2 == "SIO" .OR. ctyp2 == "CON") THEN
                ctyp = "CON"
            ELSEIF (ctyp2 == "SWA") THEN
                ctyp = "REF"
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF

        CASE("FLOW")
            ! CAUTION: this part is not complete!
            ibocd = 1
            CALL get_bc_ctyp(ctyp1, ibocd, iface, igrid)
            !ibocd = 2
            !CALL get_bc_ctyp(ctyp2, ibocd, iface, igrid)

            IF (ctyp1 == "CON" .OR. ctyp1 == "PER") THEN !.OR. ctyp1 == "FIX" .OR. ctyp1 == "OP1"
                ctyp = "CON"
            ELSEIF (ctyp1 == "SLI" .OR. ctyp1 == "NOS") THEN
                ctyp = "REF"
            ELSE
                CALL errr(__FILE__, __LINE__)
            END IF

        END SELECT

    END SUBROUTINE get_particle_bc

END MODULE particle_boundaries_mod