MODULE particle_boundaries_mod

    USE precision_mod, ONLY: realk, intk
    USE core_mod !TODO: specify
    USE utils_mod

    USE particle_ofields_mod
    USE particle_obstacles_mod
    USE particle_runtimestat_mod, ONLY: psim_n_replaced_tot, &
     psim_n_bcerr, psim_max_bcerr

    IMPLICIT NONE

    CHARACTER(len = 4) :: bc_coupling_mode = "FLOW" ! must be "FLOW", "SCAL" or "PART"

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

        INTEGER(intk) :: face_neighbours(26)

        REAL(realk) :: face_normals(3, 26) = 0.0

    END TYPE particle_boundaries_t

    TYPE(particle_boundaries_t), ALLOCATABLE :: particle_boundaries(:)

    INTEGER(intk), PARAMETER :: ncorn = 8_intk

    TYPE :: particle_gcorner_boundaries_t

        INTEGER(intk) :: location(3)

        REAL(realk) :: face_coord(3)

        INTEGER(intk) :: face_neighbours(8)

        REAL(realk) :: face_normals(12) = 0.0

    END TYPE particle_gcorner_boundaries_t

    TYPE(particle_gcorner_boundaries_t), ALLOCATABLE :: particle_gcorner_boundaries(:)

#if defined __INTEL_COMPILER
    !$omp declare mapper(particle_boundaries_t :: bnd) map(to: bnd, bnd%face_neighbours, bnd%face_normals)
#endif

    CONTAINS

    SUBROUTINE init_particle_boundaries()

        ! local variables
        INTEGER(intk) :: igrid, icorn, iface, jface, i, j
        INTEGER(intk) :: neighbours(26)
        ! to store which boundary surfaces (faces 1-6) that define a lower order face (7-26) are of connect (or periodic) type:
        INTEGER(intk) :: connect_faces(4)
        REAL(realk) :: magnitude
        CHARACTER(len=3) :: ctyp
        LOGICAL :: found

        CALL start_timer(900)
        CALL start_timer(910)

        ALLOCATE(particle_boundaries(ngrid))

        ! CON : connective particle boundary (can be periodic)
        ! REF : reflective particle boundary
        DO igrid = 1, ngrid

            ! TODO: possibly cycle if grid is not on particle_level
            !IF (level(igrid) /= particle_level) CYCLE

            CALL get_neighbours(neighbours, igrid)

            DO iface = 1, 6

                CALL get_particle_bc(igrid, iface, bc_coupling_mode, ctyp)

                particle_boundaries(igrid)%face_neighbours(iface) = neighbours(iface)

                SELECT CASE(iface)
                    CASE(1)
                        IF (ctyp == "REF") THEN ! ctyp == "SWA"
                            particle_boundaries(igrid)%face_normals(1, 1) = 1.0
                            particle_boundaries(igrid)%face_neighbours(iface) = igrid
                        END IF
                    CASE(2)
                        IF (ctyp == "REF") THEN
                            particle_boundaries(igrid)%face_normals(1, 2) = -1.0
                            particle_boundaries(igrid)%face_neighbours(iface) = igrid
                        END IF
                    CASE(3)
                        IF (ctyp == "REF") THEN
                            particle_boundaries(igrid)%face_normals(2, 3) = 1.0
                            particle_boundaries(igrid)%face_neighbours(iface) = igrid
                        END IF
                    CASE(4)
                        IF (ctyp == "REF") THEN
                            particle_boundaries(igrid)%face_normals(2, 4) = -1.0
                            particle_boundaries(igrid)%face_neighbours(iface) = igrid
                        END IF
                    CASE(5)
                        IF (ctyp == "REF") THEN
                            particle_boundaries(igrid)%face_normals(3, 5) = 1.0
                            particle_boundaries(igrid)%face_neighbours(iface) = igrid
                        END IF
                    CASE(6)
                        IF (ctyp == "REF") THEN
                            particle_boundaries(igrid)%face_normals(3, 6) = -1.0
                            particle_boundaries(igrid)%face_neighbours(iface) = igrid
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
                            particle_boundaries(igrid)%face_normals(j, iface) = &
                            particle_boundaries(igrid)%face_normals(j, iface) + particle_boundaries(igrid)%face_normals(j, facelist_b(i, iface))
                        END DO

                    END IF
                END DO

                particle_boundaries(igrid)%face_neighbours(iface) = igrid
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
                            particle_boundaries(igrid)%face_neighbours(iface) = neighbours(jface)
                        END IF
                    END IF
                END DO

                magnitude = SQRT(particle_boundaries(igrid)%face_normals(1, iface)**2 + &
                    particle_boundaries(igrid)%face_normals(2, iface)**2 + &
                    particle_boundaries(igrid)%face_normals(3, iface)**2)

                DO j = 1, 3
                    ! EPSILON(magnitude) is an arbitrary value significantely smaller than 1 as
                    ! the shortest valid normal vector up to now should have a magnitude of 1
                    IF (magnitude <= EPSILON(magnitude)) THEN
                        particle_boundaries(igrid)%face_normals(j, iface) = 0.0
                    ELSE
                        particle_boundaries(igrid)%face_normals(j, iface) = particle_boundaries(igrid)%face_normals(j, iface) / magnitude
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
                        WRITE(*, *) "Neigbhour grid:       ", particle_boundaries(igrid)%face_neighbours(iface)
                        WRITE(*,'("Normal vector: ")')
                        WRITE(*, *) "n1                    ", particle_boundaries(igrid)%face_normals(1, iface)
                        WRITE(*, *) "n2                    ", particle_boundaries(igrid)%face_normals(2, iface)
                        WRITE(*, *) "n3                    ", particle_boundaries(igrid)%face_normals(3, iface)
                        WRITE(*, *) " "
                    END DO
                END DO
                WRITE(*, '()')
            END IF
        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

        ALLOCATE(particle_gcorner_boundaries(ngrid * 8_intk))

        DO igrid = 1, ngrid
            DO icorn = 1, ncorn
                CALL set_gcorner_boundary(igrid, icorn)
            END DO
        END DO

        !$omp target enter data map(to: ngrid)

#if defined __INTEL_COMPILER
        !$omp target enter data map(mapper(particle_boundaries_t), to: particle_boundaries(1:ngrid))
        !$omp target enter data map(particle_boundaries(1:ngrid)%face_neighbours, &
        !$omp particle_boundaries(1:ngrid)%face_normals)
#else
        !$omp target enter data map(to: particle_boundaries)
#endif
        
        CALL read_obstacles()

        !$omp target enter data map(always, to: my_obstacles_offload)
        !$omp target enter data map(always, to: n_my_obstacles_on_grid)
        !$omp target enter data map(always, to: obstacle_displ)
        !$omp target enter data map(always, to: aura)

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_boundaries

    !-----------------------------------

    SUBROUTINE finish_particle_boundaries()

        CALL start_timer(900)
        CALL start_timer(990)

        !$omp target exit data map(delete: my_obstacles_offload, obstacle_displ, n_my_obstacles_on_grid, aura)

        CALL finish_obstacles()

        !$omp target exit data map(delete: particle_boundaries) 
        DEALLOCATE(particle_boundaries)

        CALL stop_timer(990)
        CALL stop_timer(900)

    END SUBROUTINE finish_particle_boundaries

    !-----------------------------------

    ! TODO: source the following moduled out into particle_motion_mod or so

    SUBROUTINE move_particle(particle, dx, dy, dz, dx_eff, dy_eff, dz_eff, temp_coord_prev, temp_grid_prev)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        REAL(realk), INTENT(in) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_eff, dy_eff, dz_eff
        REAL(realk), INTENT(inout), OPTIONAL :: temp_coord_prev(3)
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

        IF (PRESENT(temp_coord_prev)) THEN
            x = temp_coord_prev(1)
            y = temp_coord_prev(2)
            z = temp_coord_prev(3)
        ELSE
            x = particle%x
            y = particle%y
            z = particle%z
        END IF

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

                CALL reflect_at_obstacle(x, y, z, dx_from_here, dy_from_here, dz_from_here, my_obstacles(iobst_local))

            ELSEIF (0 < iface) THEN

                destgrid = particle_boundaries(temp_grid)%face_neighbours(iface)

                IF (TRIM(particle_terminal) == "verbose") THEN
                    WRITE(*, *) "Particle moved to Grid face ", iface, "with Target Grid ", destgrid, "."
                END IF

                CALL reflect_at_boundary(dx_from_here, dy_from_here, dz_from_here, &
                 particle_boundaries(temp_grid)%face_normals(1, iface), &
                 particle_boundaries(temp_grid)%face_normals(2, iface), &
                 particle_boundaries(temp_grid)%face_normals(3, iface), reflect)

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

        IF (PRESENT(temp_coord_prev)) THEN
            temp_coord_prev(1) = x
            temp_coord_prev(2) = y
            temp_coord_prev(3) = z
        END IF

        IF (PRESENT(temp_grid_prev)) THEN
            temp_grid_prev = temp_grid
        END IF

        ! do not update the particle grid here
        ! and do not apply periodic boundaries here
        particle%x = particle%x + dx_eff
        particle%y = particle%y + dy_eff
        particle%z = particle%z + dz_eff

        !particle%xyz_abs(1) = particle%xyz_abs(1) + dx_eff
        !particle%xyz_abs(2) = particle%xyz_abs(2) + dy_eff
        !particle%xyz_abs(3) = particle%xyz_abs(3) + dz_eff

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
        REAL(realk) :: dist, dist_to_center
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: lx, ly, lz, rx, ry, rz
        REAL(realk) :: s, sa, sb, sc, sd, a, b0, b, c0, c, d, cx, cy, cz, r

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

        s = 1.0

        ! first coefficient
        a = (dx**2 + dy**2 + dz**2)

        IF (a < EPSILON(0.0_realk)) THEN

            ! iterate over all obstacles of the grid
            IF (dread_obstacles_dict) THEN
                nobst = SIZE(my_obstacle_pointers(temp_grid)%grid_obstacles)
            ELSE
                nobst = 0
            END IF

            b0 = 2*x*dx + 2*y*dy + 2*z*dz
            c0 = x**2 + y**2 + z**2

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
                b = b0 - 2*cx*dx - 2*cy*dy - 2*cz*dz
                c = c0 + cx**2 + cy**2 + cz**2 - 2*x*cx - 2*y*cy - 2*z*cz - r**2
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

        END IF

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


    SUBROUTINE move_particle_target(particle, dx, dy, dz, dx_eff, dy_eff, dz_eff, temp_x, temp_y, temp_z, temp_grid_prev, boundaries, obstacles, dreplace)

        !$omp declare target

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        REAL(realk), INTENT(in) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_eff, dy_eff, dz_eff
        REAL(realk), INTENT(inout) :: temp_x, temp_y, temp_z
        INTEGER(intk), INTENT(inout) :: temp_grid_prev
        TYPE(particle_boundaries_t), INTENT(in) :: boundaries(ngrid)
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles
        LOGICAL, INTENT(inout) :: dreplace 
        
        ! local variables
        INTEGER(intk) :: temp_grid, iface, iobst_local, destgrid, i
        INTEGER(intk) :: reflect(3)
        REAL(realk) :: x, y, z
        REAL(realk) :: dx_step, dy_step, dz_step
        REAL(realk) :: dx_from_here, dy_from_here, dz_from_here
        REAL(realk) :: bbox(6)

        dreplace = .FALSE.

        dx_eff = 0.0
        dy_eff = 0.0
        dz_eff = 0.0

        temp_grid = temp_grid_prev

        x = temp_x
        y = temp_y
        z = temp_z

        dx_from_here = dx
        dy_from_here = dy
        dz_from_here = dz

        iobst_local = 0

        ! to avoid branch divergence here, just iterate to the max. number of iterations that would be a stoping criterion anyways
        DO i = 1, 10

            CALL get_bbox_target(bbox(1), bbox(2), bbox(3), bbox(4), bbox(5), bbox(6), temp_grid)

            CALL move_to_boundary_target(particle%igrid, temp_grid, x, y, z, &
             dx_from_here, dy_from_here, dz_from_here, dx_step, dy_step, dz_step, iface, iobst_local, obstacles, bbox)

            ! TODO: reactivate particle replacement (?)
            ! replace current particle coordinates by a random valid position on the particles curren grid
            IF (dreplace) THEN
                RETURN
            END IF

            dx_eff = dx_eff + dx_step
            dy_eff = dy_eff + dy_step
            dz_eff = dz_eff + dz_step

            IF (0 < iobst_local) THEN

                CALL reflect_at_obstacle(x, y, z, dx_from_here, dy_from_here, dz_from_here, obstacles(iobst_local))

            ELSEIF (0 < iface) THEN

                destgrid = boundaries(temp_grid)%face_neighbours(iface)

                CALL reflect_at_boundary(dx_from_here, dy_from_here, dz_from_here, &
                 boundaries(temp_grid)%face_normals(1, iface), &
                 boundaries(temp_grid)%face_normals(2, iface), &
                 boundaries(temp_grid)%face_normals(3, iface), reflect)

                CALL update_coordinates_target(temp_grid, destgrid, iface, x, y, z, bbox, reflect)

                temp_grid = destgrid

            END IF

        END DO

        temp_x = x
        temp_y = y
        temp_z = z

        temp_grid_prev = temp_grid

        ! do not update the particle grid here
        ! and do not apply periodic boundaries here
        particle%x = particle%x + dx_eff
        particle%y = particle%y + dy_eff
        particle%z = particle%z + dz_eff

        !particle%xyz_abs(1) = particle%xyz_abs(1) + dx_eff
        !particle%xyz_abs(2) = particle%xyz_abs(2) + dy_eff
        !particle%xyz_abs(3) = particle%xyz_abs(3) + dz_eff

    END SUBROUTINE move_particle_target

    !-----------------------------------

   ! This subroutine only considers grids on the same level
    ! CAUTION: Here, temp_grid refers to the grid the particle coordinates are currently on and of which the boundaries are relevant.
    ! This might NOT be particle%igrid, which is used to deduce the velocity
    SUBROUTINE move_to_boundary_target(old_grid, temp_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, iface, iobst_local, obstacles, bbox)

        !$omp declare target

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: old_grid
        INTEGER(intk), INTENT(in) :: temp_grid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        INTEGER(intk), INTENT(out) :: iface
        INTEGER(intk), INTENT(inout) :: iobst_local
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles
        REAL(realk), INTENT(in) :: bbox(6)

        !local variables
        REAL(realk) :: s, dist
        LOGICAL :: dget_exit_face

        dx_to_b = 0.0
        dy_to_b = 0.0
        dz_to_b = 0.0
        
        CALL s_to_obstacle(old_grid, x, y, z, dx, dy, dz, iobst_local, s, obstacles)

        IF (s <= 0.0_realk) THEN
            iface = 0
            RETURN
        END IF

        CALL to_grid_boundary(temp_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, s, dget_exit_face, bbox)

        IF (dget_exit_face) THEN
            ! get face after x/y/z have (potentially) been altered
            iobst_local = 0
            CALL get_exit_face_target(bbox, x, y, z, dist, iface)
        ELSE
            iface = 0
        END IF

    END SUBROUTINE move_to_boundary_target


    SUBROUTINE s_to_obstacle(old_grid, x, y, z, dx, dy, dz, iobst_local, s, obstacles)

        !$omp declare target

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: old_grid
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(in) :: dx, dy, dz
        INTEGER(intk), INTENT(inout) :: iobst_local
        REAl(realk), INTENT(inout) :: s
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles

        !local variables
        INTEGER(intk) :: i, nobst
        REAL(realk) :: sa, sb, sc, sd, a, b, b0, c, c0, d, r

        ! STEP 1 - OBSTACLES
        ! find intersection points of the line the particle moves on (straight) and the sphere surface
            ! particle path: X(s) = X + dX * s with s: [0, 1] (X is the vector (x/y/z))
            ! => |X + dX * s - C| = r (C is the sphere center (cx/cy/cz))
            ! => (x + dx * s -cx)² + (y + dy * s -cy)² + (z + dz * s -cz)² = r² (r is the sphere radius)
            ! => s1/s2 = sa/sb = (-b +/- sqrt(b² - 4ac)) / 2a (corefficients see code)

        s = 1.0

        ! first coefficient
        a = (dx**2 + dy**2 + dz**2)

        b0 = 2*x*dx + 2*y*dy + 2*z*dz
        c0 = x**2 + y**2 + z**2

        ! iterate over all obstacles of the grid
        nobst = n_my_obstacles_on_grid(old_grid) * a_greater_b(a, 0.0_realk)

        DO i = 1, nobst

            ! check if a particle interacts with the obstacle it has been deflected from in the previous timestep
            IF (i == iobst_local .OR. obstacles(i)%iobst < 0) THEN
                CYCLE
            END IF

            b = b0 - &
                2*obstacles(i)%x*dx - &
                2*obstacles(i)%y*dy - &
                2*obstacles(i)%z*dz
            c = c0 + &
                obstacles(i)%x**2 + &
                obstacles(i)%y**2 + &
                obstacles(i)%z**2 - &
                obstacles(i)%x - &
                obstacles(i)%y - &
                obstacles(i)%z - &
                obstacles(i)%radius**2
            d = b**2 - 4*a*c

            IF (d < EPSILON(0.0_realk)) THEN
                CYCLE
            END IF

            sa = (-b + SQRT(d)) / 2 / a
            sb = (-b - SQRT(d)) / 2 / a

            ! if a particle moves towards an obstacle, limit its motion to the closest intersection yet
            IF (sa >= 0.0 .AND. sb >= 0.0) THEN
                sc = MIN(sa, sb)
                IF (sc < s) THEN
                    s = sc
                    iobst_local = i
                END IF
            ! elseif a particle moves away from the current obstacle, cycle
            ELSEIF (sa <= 0.0 .AND. sb <= 0.0) THEN
                CYCLE
            ! else (if sa < 0 and sb > 0 or vice versa) the particle is inside the current obstacle
            ! => replace current particle coordinates by a random valid position on the particles curren grid
            ELSE
                !dist_to_center = SQRT((x - cx)**2 + (y - cy)**2 + (z - cz)**2)

                !IF ((r - dist_to_center) > aura) THEN
                !    replace = .TRUE.
                !    iface = 0
                !    iobst_local = 0
                !    RETURN
                !END IF

                sc = MIN(sa, sb)
                sd = MAX(sa, sb)

                IF (ABS(sc) < ABS(sd)) THEN
                    s = 0.0
                    iobst_local = i
                    EXIT
                ELSEIF (ABS(sc) >= ABS(sd)) THEN
                    CYCLE
                END IF

            END IF
        END DO

    END SUBROUTINE s_to_obstacle


    SUBROUTINE to_grid_boundary(temp_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, s, dget_exit_face, bbox)

        !$omp declare target

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: temp_grid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        REAL(realk), INTENT(inout) :: s
        LOGICAL, INTENT(out) :: dget_exit_face
        REAL(realk), INTENT(in) :: bbox(6)

        !local variables
        REAL(realk) :: lx, ly, lz, rx, ry, rz

        ! STEP 2 - GRID BOUNDARIES
        ! now check if any grid boundary is reached before any obstacle is reached
        IF (dx < 0) THEN
            lx = (bbox(1) - x)
            ! if particle is at boundary in X dir (lx = 0) or particle is outside temp_grid (lx > 0.0),
            ! get exit face and return; so if a particle is incorrectly outside a reflect boundary its
            ! motion vector is reflected towards temp_grid
            IF (lx >= 0.0_realk) THEN
                ! if a particle is already on a face (esp. edge or corner), its future coordinates have to be
                ! "projected" to assign the right ecit face (otherwise, particles might get stuck on edges or corners)
                dget_exit_face = .TRUE.
                RETURN
            END IF
            rx = dx * s / lx
        ELSEIF (0 < dx) THEN
            lx = (bbox(2) - x)
            IF (lx <= 0.0_realk) THEN
                dget_exit_face = .TRUE.
                RETURN
            END IF
            rx = dx * s / lx
        ELSE
            rx = 0.0_realk
        END IF

        IF (dy < 0) THEN
            ly = (bbox(3) - y)
            IF (ly >= 0.0_realk) THEN
                dget_exit_face = .TRUE.
                RETURN
            END IF
            ry = dy * s / ly
        ELSEIF (0 < dy) THEN
            ly = (bbox(4) - y)
            IF (ly <= 0.0_realk) THEN
                dget_exit_face = .TRUE.
                RETURN
            END IF
            ry = dy * s / ly
        ELSE
            ry = 0.0_realk
        END IF

        IF (dz < 0) THEN
            lz = (bbox(5) - z)
            IF(lz >= 0.0_realk) THEN
                dget_exit_face = .TRUE.
                RETURN
            END IF
            rz = dz * s / lz
        ELSEIF (0 < dz) THEN
            lz = (bbox(6) - z)
            IF(lz <= 0.0_realk) THEN
                dget_exit_face = .TRUE.
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

            dget_exit_face = .FALSE.
            RETURN
        END IF

        dget_exit_face = .TRUE.

        IF (dx < 0 .AND. ry <= rx .AND. rz <= rx) THEN

            dx_to_b = lx
            dy_to_b = (lx * dy/dx)
            dz_to_b = (lx * dz/dx)
            x = bbox(1) ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
            y = y + dy_to_b
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dx .AND. ry <= rx .AND. rz <= rx) THEN

            dx_to_b = lx
            dy_to_b = (lx * dy/dx)
            dz_to_b = (lx * dz/dx)
            x = bbox(2) ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
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
            y = bbox(3) ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
            z = z + dz_to_b
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dy .AND. rx < ry .AND. rz <= ry) THEN

            dx_to_b = (ly * dx/dy)
            dy_to_b = ly
            dz_to_b = (ly * dz/dy)
            x = x + dx_to_b
            y = bbox(4) ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
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
            z = bbox(5) ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dz .AND. rx < rz .AND. ry < rz) THEN

            dx_to_b = (lz * dx/dz)
            dy_to_b = (lz * dy/dz)
            dz_to_b = lz
            x = x + dx_to_b
            y = y + dy_to_b
            z = bbox(6) ! keep this expression so no floating point errors occur and the particle is EXACTLY at the boundary
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        END IF

    END SUBROUTINE to_grid_boundary


    SUBROUTINE move_to_boundary_target2(temp_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, iface, iobst_local, replace, obstacles, bbox)

        !$omp declare target

        ! subroutine arguments
        INTEGER(intk), INTENT(inout) :: temp_grid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        INTEGER(intk), INTENT(out) :: iface
        INTEGER(intk), INTENT(inout) :: iobst_local
        LOGICAL, INTENT(out) :: replace
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles
        REAL(realk), INTENT(in) :: bbox(6)

        !local variables
        INTEGER(intk) :: closestbx, closestby, closestbz, sum
        REAL(realk) :: dist, s
        REAL(realk) :: lx, ly, lz
        REAL(realk) :: newcoord(3)

        replace = .FALSE.
        s = 1.0_realk

        dx_to_b = 0.0
        dy_to_b = 0.0
        dz_to_b = 0.0

        CALL s_to_obstacle2(temp_grid, x, y, z, dx, dy, dz, iobst_local, s, obstacles)

        CALL to_grid_boundary2(temp_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, s, iobst_local, bbox)

        CALL get_exit_face_target(bbox, x, y, z, dist, iface)

    END SUBROUTINE move_to_boundary_target2


    SUBROUTINE s_to_obstacle2(temp_grid, x, y, z, dx, dy, dz, iobst_local, s, obstacles)

        !$omp declare target

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: temp_grid
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(in) :: dx, dy, dz
        INTEGER(intk), INTENT(inout) :: iobst_local
        REAl(realk), INTENT(inout) :: s
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles

        !local variables
        INTEGER(intk) :: i, nobst
        REAL(realk) :: sa, sb, a, b, b0, c, c0, d, r
        REAL(realk) :: cond0, cond1, cond2

        ! STEP 1 - OBSTACLES

        ! find intersection points of the line the particle moves on (straight) and the sphere surface
            ! particle path: X(s) = X + dX * s with s: [0, 1] (X is the vector (x/y/z))
            ! => |X + dX * s - C| = r (C is the sphere center (cx/cy/cz))
            ! => (x + dx * s -cx)² + (y + dy * s -cy)² + (z + dz * s -cz)² = r² (r is the sphere radius)
            ! => s1/s2 = sa/sb = (-b +/- sqrt(b² - 4ac)) / 2a (corefficients see code)

        ! first coefficient
        a = (dx**2 + dy**2 + dz**2)
        b0 = 2*x*dx + 2*y*dy + 2*z*dz
        c0 = x**2 + y**2 + z**2

        ! iterate over all obstacles of the grid
        nobst = n_my_obstacles_on_grid(temp_grid)
        DO i = 1, nobst

            ! check if a particle interacts with the obstacle it has been deflected from in the previous timestep
            cond0 = a_greater_b(i, iobst_local) &
                    + a_greater_b(iobst_local, i)

            ! sphere dependent coefficients
            b = b0 - &
                2*obstacles(i)%x*dx - &
                2*obstacles(i)%y*dy - &
                2*obstacles(i)%z*dz
            c = c0 + &
                obstacles(i)%x**2 + &
                obstacles(i)%y**2 + &
                obstacles(i)%z**2 - &
                2*x*obstacles(i)%x - &
                2*y*obstacles(i)%y - &
                2*z*obstacles(i)%z - &
                obstacles(i)%radius**2
            d = b**2 - 4*a*c

            sa = (-b + SQRT(MAX(0.0_realk, d))) / 2 / MAX(EPSILON(a), a)
            sb = (-b - SQRT(MAX(0.0_realk, d))) / 2 / MAX(EPSILON(a), a)

            cond1 = a_greaterequal_b(sa, 0.0_realk) * a_greaterequal_b(sb, 0.0_realk) &
                    * a_greater_b(a, 0.0_realk) * a_greater_b(d, 0.0_realk)

            cond2 = a_greater_b(ABS(MAX(sa, sb)), ABS(MIN(sa, sb))) &
                    * a_greater_b(a, 0.0_realk) * a_greater_b(d, 0.0_realk)

            s = s - (s - MIN(sa, sb)) * cond1 - s * cond2

            iobst_local = iobst_local * NINT(1.0_realk - cond1) &
                 + i * NINT(cond1)

            iobst_local = iobst_local * NINT(1.0_realk - cond2) &
                 + i * NINT(cond2)
        
        END DO

    END SUBROUTINE s_to_obstacle2


    SUBROUTINE to_grid_boundary2(temp_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, s, iobst_local, bbox)

        !$omp declare target

        ! subroutine arguments
        INTEGER(intk), INTENT(inout) :: temp_grid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        REAL(realk), INTENT(inout) :: s
        INTEGER(intk), INTENT(inout) :: iobst_local
        REAL(realk), INTENT(in) :: bbox(6)

        !local variables
        INTEGER(intk) :: closestbx, closestby, closestbz, sum
        REAL(realk) :: lx, ly, lz, dist
        REAL(realk) :: newcoord(3)

        ! STEP 2 - GRID BOUNDARIES
        ! signed distance of particle to grid boundaries
        lx = 0.0_realk + (bbox(1) - x) * a_greater_b(0.0_realk, dx) + (bbox(2) - x) * a_greater_b(dx, 0.0_realk) ! = lx
        ly = 0.0_realk + (bbox(3) - y) * a_greater_b(0.0_realk, dy) + (bbox(4) - y) * a_greater_b(dy, 0.0_realk) ! = ly
        lz = 0.0_realk + (bbox(5) - z) * a_greater_b(0.0_realk, dz) + (bbox(6) - z) * a_greater_b(dz, 0.0_realk) ! = lz
        
        closestbx = 0
        closestby = 0
        closestbz = 0

        ! check which face is hit first
        IF (ABS(dx) > 0.0_realk) THEN
            closestbx = MAX(0_intk, CEILING(s - MAX(0.0_realk, lx/dx)) * INT(SIGN(1.0_realk, lx/dx)))
            s = MIN(s, MAX(0.0_realk, lx/dx))
        END IF

        IF (ABS(dy) > 0.0_realk) THEN
            closestby = MAX(0_intk, CEILING(s - MAX(0.0_realk, ly/dy)) * INT(SIGN(1.0_realk, ly/dy)))
            closestbx = MAX(0_intk, closestbx - closestby)
            s = MIN(s, MAX(0.0_realk, ly/dy))
        END IF

        IF (ABS(dz) > 0.0_realk) THEN
            closestbz = MAX(0_intk, CEILING(s - MAX(0.0_realk, lz/dz)) * INT(SIGN(1.0_realk, lz/dz)))
            closestby = MAX(0_intk, closestby - closestbz)
            closestbx = MAX(0_intk, closestbx - closestbz)
            s = MIN(s, MAX(0.0_realk, lz/dz))
        END IF

        dx_to_b = 0.0_realk + dx * s
        dy_to_b = 0.0_realk + dy * s
        dz_to_b = 0.0_realk + dz * s

        dx = dx - dx_to_b
        dy = dy - dy_to_b
        dz = dz - dz_to_b

        newcoord(1) = bbox(1)
        newcoord(2)  = x + dx_to_b
        newcoord(3)  = bbox(2)
        x = newcoord(2 + closestbx * INT(SIGN(1.0_realk, dx)))
        
        newcoord(1) = bbox(3)
        newcoord(2)  = y + dy_to_b
        newcoord(3)  = bbox(4)
        y = newcoord(2 + closestby * INT(SIGN(1.0_realk, dy)))
        
        newcoord(1) = bbox(5)
        newcoord(2) = z + dz_to_b
        newcoord(3) = bbox(6)
        z = newcoord(2 + closestbz * INT(SIGN(1.0_realk, dz)))

        sum = closestbx + closestby + closestbz
        IF (sum > 0) THEN
            iobst_local = 0_intk
        END IF

    END SUBROUTINE to_grid_boundary2


    SUBROUTINE move_particle_target3(particle, pstag, dx, dy, dz, dx_eff, dy_eff, dz_eff, gcorner_boundary, obstacles, dreplace)

        !$omp declare target

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(inout) :: pstag(3)
        REAL(realk), INTENT(in) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_eff, dy_eff, dz_eff
        TYPE(particle_gcorner_boundaries_t), INTENT(in) :: gcorner_boundary 
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles
        LOGICAL, INTENT(inout) :: dreplace 
        
        ! local variables
        INTEGER(intk) :: temp_grid, idir, iobst_local, destgrid, i
        INTEGER(intk) :: reflect(3), pstag_counter(3)
        REAL(realk) :: n(3)
        REAL(realk) :: x, y, z
        REAL(realk) :: dx_step, dy_step, dz_step
        REAL(realk) :: dx_from_here, dy_from_here, dz_from_here

        dreplace = .FALSE.

        dx_eff = 0.0
        dy_eff = 0.0
        dz_eff = 0.0

        x = particle%x
        y = particle%y
        z = particle%z

        dx_from_here = dx
        dy_from_here = dy
        dz_from_here = dz

        iobst_local = 0
        idir = 0
        pstag_counter = 0

        ! to avoid branch divergence here, just iterate to the max. number of iterations that would be a stoping criterion anyways
        DO i = 1, 10

            CALL move_to_boundary_target3(gcorner_boundary, pstag, particle%igrid, x, y, z, &
             dx_from_here, dy_from_here, dz_from_here, dx_step, dy_step, dz_step, idir, iobst_local, obstacles)

            dx_eff = dx_eff + dx_step
            dy_eff = dy_eff + dy_step
            dz_eff = dz_eff + dz_step

            IF (0 < iobst_local) THEN

                CALL reflect_at_obstacle(x, y, z, dx_from_here, dy_from_here, dz_from_here, obstacles(iobst_local))

            ELSEIF (0 < idir) THEN
                
                CALL get_gcorner_normal(gcorner_boundary, pstag, idir, n)

                n(idir) = n(idir) * ((-1_intk) ** pstag(idir))

                CALL reflect_at_boundary(dx_from_here, dy_from_here, dz_from_here, &
                 n(1), n(2), n(3))

                !update pstag
                pstag(idir) = pstag(idir) + (SIGN(1_intk - INT(ABS(n(idir))), -1_intk)) ** pstag_counter(idir)
                pstag_counter(idir) = pstag_counter(idir) + 1_intk
            END IF

        END DO

        ! do not update the particle grid here
        ! and do not apply periodic boundaries here
        particle%x = particle%x + dx_eff
        particle%y = particle%y + dy_eff
        particle%z = particle%z + dz_eff

        !particle%xyz_abs(1) = particle%xyz_abs(1) + dx_eff
        !particle%xyz_abs(2) = particle%xyz_abs(2) + dy_eff
        !particle%xyz_abs(3) = particle%xyz_abs(3) + dz_eff

    END SUBROUTINE move_particle_target3


    SUBROUTINE move_to_boundary_target3(gcorner_boundary, pstag, old_grid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, idir, iobst_local, obstacles)

        !$omp declare target

        ! subroutine arguments
        TYPE(particle_gcorner_boundaries_t), INTENT(in) :: gcorner_boundary 
        INTEGER(intk), INTENT(in) :: pstag(3)
        INTEGER(intk), INTENT(in) :: old_grid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        INTEGER(intk), INTENT(out) :: idir
        INTEGER(intk), INTENT(inout) :: iobst_local
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles

        !local variables
        REAL(realk) :: s

        dx_to_b = 0.0
        dy_to_b = 0.0
        dz_to_b = 0.0
        
        CALL s_to_obstacle(old_grid, x, y, z, dx, dy, dz, iobst_local, s, obstacles)

        IF (s <= 0.0_realk) THEN
            idir = 0
            RETURN
        END IF

        CALL to_grid_boundary3(gcorner_boundary, pstag, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, s, idir, iobst_local)

    END SUBROUTINE move_to_boundary_target3


    SUBROUTINE to_grid_boundary3(gcorner_boundary, pstag, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, s, idir, iobst_local)

        !$omp declare target

        ! subroutine arguments
        TYPE(particle_gcorner_boundaries_t), INTENT(in) :: gcorner_boundary
        INTEGER(intk), INTENT(in) :: pstag(3)
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        REAL(realk), INTENT(inout) :: s
        INTEGER(intk), INTENT(out) :: idir
        INTEGER(intk), INTENT(inout) :: iobst_local

        !local variables
        REAL(realk) :: lx_a, ly_a, lz_a, rx, ry, rz

        ! STEP 2 - GRID BOUNDARIES

        idir = 0_intk

        ! abs distance of particle to grid boundaries
        lx_a = ABS(x - gcorner_boundary%face_coord(1))
        ly_a = ABS(y - gcorner_boundary%face_coord(2))
        lz_a = ABS(z - gcorner_boundary%face_coord(3))
        
        ! relative distance to boundaries (negative value indicates particle is moving away from boundary)
        rx = gcorner_boundary%location(1) * ((-1_intk) ** pstag(1)) * (s * dx / lx_a)
        ry = gcorner_boundary%location(2) * ((-1_intk) ** pstag(2)) * (s * dy / ly_a)
        rz = gcorner_boundary%location(3) * ((-1_intk) ** pstag(3)) * (s * dz / lz_a)

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
            RETURN
        END IF

        iobst_local = 0_intk

        IF (ry <= rx .AND. rz <= rx) THEN

            idir = 1_intk

            dx_to_b = gcorner_boundary%location(1) * lx_a
            dy_to_b = (lx_a / ABS(dx) * dy)
            dz_to_b = (lx_a / ABS(dx) * dz)
            ! avoid floating point errors and put the particle EXACTLY at the boundary
            x = gcorner_boundary%face_coord(1) 
            ! avoid particles to hit corners or edges => - EPSILON(x_i) * gcorner_boundary%location(i)
            y = y + dy_to_b - EPSILON(y) * gcorner_boundary%location(2)
            z = z + dz_to_b - EPSILON(z) * gcorner_boundary%location(3)
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (rx < ry .AND. rz <= ry) THEN

            idir = 2_intk

            dx_to_b = (ly_a / ABS(dy) * dx)
            dy_to_b = gcorner_boundary%location(2) * ly_a
            dz_to_b = (ly_a / ABS(dy) * dz)
            ! avoid particles to hit corners or edges => - EPSILON(x_i) * gcorner_boundary%location(i)
            x = x + dx_to_b - EPSILON(x) * gcorner_boundary%location(1)
            ! avoid floating point errors and put the particle EXACTLY at the boundary
            y = gcorner_boundary%face_coord(2) 
            z = z + dz_to_b - EPSILON(z) * gcorner_boundary%location(3)
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (rx < rz .AND. ry < rz) THEN

            idir = 3_intk

            dx_to_b = (lz_a / ABS(dz) * dx)
            dy_to_b = (lz_a / ABS(dz) * dy)
            dz_to_b = gcorner_boundary%location(3) * lz_a
            ! avoid particles to hit corners or edges => - EPSILON(x_i) * gcorner_boundary%location(i)
            x = x + dx_to_b - EPSILON(x) * gcorner_boundary%location(1)
            y = y + dy_to_b - EPSILON(y) * gcorner_boundary%location(2)
            ! avoid floating point errors and put the particle EXACTLY at the boundary
            z = gcorner_boundary%face_coord(3) 
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        END IF

    END SUBROUTINE to_grid_boundary3


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

    SUBROUTINE replace_particle_target(particle, obstacles, kk, jj, ii, x, y, z, dx, dy, dz)

        !$omp declare target

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        TYPE(obstacle_t), POINTER, CONTIGUOUS, DIMENSION(:), INTENT(in) :: obstacles
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk), dx(ii), dy(jj), dz(kk)

        ! local variables
        LOGICAL :: valid_location
        INTEGER(intk) :: i, igrid
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, x_new, y_new, z_new, dist_to_center

        ! for readability
        igrid = particle%igrid

        CALL get_bbox_target(minx, maxx, miny, maxy, minz, maxz, igrid)

        valid_location = .FALSE.
        DO WHILE (.NOT. valid_location)

            valid_location = .TRUE.

#if defined __GFORTRAN__
            CALL RANDOM_NUMBER(x_new)
            CALL RANDOM_NUMBER(y_new)
            CALL RANDOM_NUMBER(z_new)
#else
            CALL lcg(particle%seed, x_new)
            CALL lcg(particle%seed, y_new)
            CALL lcg(particle%seed, z_new)
#endif

            x_new = minx + x_new * (maxx - minx)
            y_new = miny + y_new * (maxy - miny)
            z_new = minz + z_new * (maxz - minz)

            DO i = 1, n_my_obstacles_on_grid(igrid)

                dist_to_center = SQRT((obstacles(i)%x - x_new)**2 + &
                 (obstacles(i)%y - y_new)**2 + &
                 (obstacles(i)%z - z_new)**2)

                IF (dist_to_center < obstacles(i)%radius + EPSILON(dist_to_center)) THEN
                    valid_location = .FALSE.
                    EXIT
                ELSE
                    CONTINUE
                END IF

            END DO

        END DO

        particle%x = x_new
        particle%y = y_new
        particle%z = z_new

        CALL set_particle_cell_target(particle, kk, jj, ii, x, y, z, dx, dy, dz)

    END SUBROUTINE replace_particle_target


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
        
        !$omp declare target
        
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
    SUBROUTINE reflect_at_obstacle(x, y, z, dx, dy, dz, obstacle)

        !$omp declare target

        ! Presumption 1: Particle is already exactly on the boundary!
        ! Presumption 2: Obstacle is a sphere!

        ! subroutine arguments
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        TYPE(obstacle_t) :: obstacle

        ! local variables
        REAL(realk) :: n1, n2 , n3, magnitude

        n1 = x - obstacle%x
        n2 = y - obstacle%y
        n3 = z - obstacle%z

        magnitude = SQRT(n1**2 + n2**2 + n3**2)

        n1 = n1 / magnitude
        n2 = n2 / magnitude
        n3 = n3 / magnitude

        CALL reflect_at_boundary(dx, dy, dz, n1, n2, n3)

    END SUBROUTINE reflect_at_obstacle


    SUBROUTINE get_gcorner_neighbour(gcorner_boundary, stag, neighbour)

        ! subroutine arguments
        TYPE(particle_gcorner_boundaries_t), INTENT(in) :: gcorner_boundary
        INTEGER(intk), INTENT(in) :: stag(3)
        INTEGER(intk), INTENT(out) :: neighbour

        ! local variables
        INTEGER(intk) :: ineighbour  

        ineighbour = stag(1) * 4_intk + stag(2) * 2_intk + stag(3) * 1_intk + 1_intk

        neighbour = gcorner_boundary%face_neighbours(ineighbour)

    END SUBROUTINE get_gcorner_neighbour


    SUBROUTINE get_gcorner_normal(gcorner_boundary, stag, ndir, n)

        ! subroutine arguments
        TYPE(particle_gcorner_boundaries_t), INTENT(in) :: gcorner_boundary
        INTEGER(intk), INTENT(in) :: stag(3), ndir
        REAL(realk), INTENT(out) :: n(3)

        ! local variables 
        INTEGER(intk) :: inormal
        REAL(realk) :: nfactor

        n = 0.0_realk

        inormal = 4 * ndir + stag(MOD(ndir + 2_intk, 3_intk)) * 2 + stag(MOD(ndir + 1_intk, 3_intk)) * 1 + 1

        nfactor = SIGN(1.0_realk, (0.5_realk - 1.0_realk * REAL(stag(ndir), realk)))

        n(ndir) = nfactor * gcorner_boundary%face_normals(inormal)

    END SUBROUTINE get_gcorner_normal


    SUBROUTINE set_gcorner_boundary(igrid, icorn) 

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, icorn

        ! local variables 
        INTEGER(intk) :: jgrid, iface, idir, iface_of_corner
        INTEGER(intk) :: neighbours(26)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        iface_of_corner = 18_intk + icorn

        CALL get_neighbours(neighbours, igrid)

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        ! TODO: unify convention for neighbours an normals!

        ! mapping the grid wise particle_boundaries onto the particle_gcorner_boundaries
        SELECT CASE(iface_of_corner)
            CASE(19_intk)

                ! LOCATION
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(1) = -1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(2) = -1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(3) = -1_intk

                ! FACE COORDINATES
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(1) = minx
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(2) = miny
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(3) = minz

                ! NEIGBHOURS
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(1) = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(2) = particle_boundaries(igrid)%face_neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(3) = particle_boundaries(igrid)%face_neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(4) = particle_boundaries(igrid)%face_neighbours(15)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(5) = particle_boundaries(igrid)%face_neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(6) = particle_boundaries(igrid)%face_neighbours(9)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(7) = particle_boundaries(igrid)%face_neighbours(7)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(8) = particle_boundaries(igrid)%face_neighbours(19)

                ! FACE NORMALS
                ! x-direction (1-4)
                idir = 1_intk
                iface = 1_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 1) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 2) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 3) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(15)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 4) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! y-direction (5-8) 
                idir = 2_intk
                iface = 3_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 5) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 6) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 7) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(9)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 8) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! z-direction (9-12) 
                idir = 3_intk
                iface = 5_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 9) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(10) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(11) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(7)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(12) = particle_boundaries(jgrid)%face_normals(idir, iface)

            CASE(20_intk)

                ! LOCATION
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(1) = -1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(2) = -1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(3) =  1_intk

                ! FACE COORDINATES
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(1) = minx
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(2) = miny
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(3) = maxz

                ! NEIGBHOURS
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(1) = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(2) = particle_boundaries(igrid)%face_neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(3) = particle_boundaries(igrid)%face_neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(4) = particle_boundaries(igrid)%face_neighbours(16)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(5) = particle_boundaries(igrid)%face_neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(6) = particle_boundaries(igrid)%face_neighbours(10)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(7) = particle_boundaries(igrid)%face_neighbours(7)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(8) = particle_boundaries(igrid)%face_neighbours(20)

                ! FACE NORMALS
                ! x-direction (1-4) 
                idir = 1_intk
                iface = 1_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 1) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 2) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 3) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(16)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 4) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! y-direction (5-8) 
                idir = 2_intk
                iface = 3_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 5) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 6) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 7) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(10)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 8) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! z-direction (9-12)
                idir = 3_intk 
                iface = 6_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 9) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(10) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(11) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(7)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(12) = particle_boundaries(jgrid)%face_normals(idir, iface)

            CASE(21_intk)

                ! LOCATION
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(1) = -1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(2) =  1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(3) = -1_intk

                ! FACE COORDINATES
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(1) = minx
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(2) = maxy
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(3) = minz

                ! NEIGBHOURS
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(1) = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(2) = particle_boundaries(igrid)%face_neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(3) = particle_boundaries(igrid)%face_neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(4) = particle_boundaries(igrid)%face_neighbours(17)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(5) = particle_boundaries(igrid)%face_neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(6) = particle_boundaries(igrid)%face_neighbours(9)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(7) = particle_boundaries(igrid)%face_neighbours(8)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(8) = particle_boundaries(igrid)%face_neighbours(21)

                ! FACE NORMALS
                ! x-direction (1-4) 
                idir = 1_intk
                iface = 1_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 1) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 2) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 3) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(17)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 4) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! y-direction (5-8) 
                idir = 2_intk
                iface = 4_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 5) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 6) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 7) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(9)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 8) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! z-direction (9-12) 
                idir = 3_intk
                iface = 5_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 9) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(10) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(11) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(8)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(12) = particle_boundaries(jgrid)%face_normals(idir, iface)

            CASE(22_intk)

                ! LOCATION
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(1) = -1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(2) =  1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(3) =  1_intk

                ! FACE COORDINATES
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(1) = minx
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(2) = maxy
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(3) = maxz

                ! NEIGBHOURS
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(1) = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(2) = particle_boundaries(igrid)%face_neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(3) = particle_boundaries(igrid)%face_neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(4) = particle_boundaries(igrid)%face_neighbours(18)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(5) = particle_boundaries(igrid)%face_neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(6) = particle_boundaries(igrid)%face_neighbours(10)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(7) = particle_boundaries(igrid)%face_neighbours(8)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(8) = particle_boundaries(igrid)%face_neighbours(22)

                ! FACE NORMALS
                ! x-direction (1-4) 
                idir = 1_intk
                iface = 1_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 1) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 2) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 3) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(18)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 4) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! y-direction (5-8) 
                idir = 2_intk
                iface = 4_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 5) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 6) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 7) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(10)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 8) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! z-direction (9-12) 
                idir = 3_intk
                iface = 6_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 9) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(1)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(10) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(11) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(8)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(12) = particle_boundaries(jgrid)%face_normals(idir, iface)

            CASE(23_intk)

                ! LOCATION
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(1) =  1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(2) = -1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(3) = -1_intk

                ! FACE COORDINATES
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(1) = maxx
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(2) = miny
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(3) = minz

                ! NEIGBHOURS
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(1) = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(2) = particle_boundaries(igrid)%face_neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(3) = particle_boundaries(igrid)%face_neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(4) = particle_boundaries(igrid)%face_neighbours(15)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(5) = particle_boundaries(igrid)%face_neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(6) = particle_boundaries(igrid)%face_neighbours(13)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(7) = particle_boundaries(igrid)%face_neighbours(11)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(8) = particle_boundaries(igrid)%face_neighbours(23)

                ! FACE NORMALS
                ! x-direction (1-4) 
                idir = 1_intk
                iface = 2_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 1) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 2) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 3) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(15)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 4) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! y-direction (5-8) 
                idir = 2_intk
                iface = 3_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 5) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 6) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 7) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(13)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 8) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! z-direction (9-12) 
                idir = 3_intk
                iface = 5_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 9) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(10) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(11) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(11)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(12) = particle_boundaries(jgrid)%face_normals(idir, iface)

            CASE(24_intk)

                ! LOCATION
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(1) =  1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(2) = -1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(3) =  1_intk

                ! FACE COORDINATES
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(1) = maxx
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(2) = miny
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(3) = maxz

                ! NEIGBHOURS
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(1) = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(2) = particle_boundaries(igrid)%face_neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(3) = particle_boundaries(igrid)%face_neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(4) = particle_boundaries(igrid)%face_neighbours(16)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(5) = particle_boundaries(igrid)%face_neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(6) = particle_boundaries(igrid)%face_neighbours(14)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(7) = particle_boundaries(igrid)%face_neighbours(11)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(8) = particle_boundaries(igrid)%face_neighbours(24)

                ! FACE NORMALS
                ! x-direction (1-4) 
                idir = 1_intk
                iface = 2_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 1) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 2) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 3) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(16)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 4) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! y-direction (5-8) 
                idir = 2_intk
                iface = 3_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 5) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 6) = particle_boundaries(jgrid)%face_normals(idir, iface)


                jgrid = neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 7) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(14)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 8) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! z-direction (9-12)
                idir = 3_intk 
                iface = 6_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 9) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(10) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(3)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(11) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(11)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(12) = particle_boundaries(jgrid)%face_normals(idir, iface)

            CASE(25_intk)

                ! LOCATION
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(1) =  1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(2) =  1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(3) = -1_intk

                ! FACE COORDINATES
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(1) = maxx
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(2) = maxy
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(3) = minz

                ! NEIGBHOURS
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(1) = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(2) = particle_boundaries(igrid)%face_neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(3) = particle_boundaries(igrid)%face_neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(4) = particle_boundaries(igrid)%face_neighbours(17)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(5) = particle_boundaries(igrid)%face_neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(6) = particle_boundaries(igrid)%face_neighbours(13)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(7) = particle_boundaries(igrid)%face_neighbours(12)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(8) = particle_boundaries(igrid)%face_neighbours(25)

                ! FACE NORMALS
                ! x-direction (1-4) 
                idir = 1_intk
                iface = 2_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 1) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 2) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 3) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(17)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 4) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! y-direction (5-8)
                idir = 2_intk 
                iface = 4_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 5) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(5)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 6) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 7) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(13)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 8) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! z-direction (9-12) 
                idir = 3_intk
                iface = 5_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 9) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(10) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(11) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(12)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(12) = particle_boundaries(jgrid)%face_normals(idir, iface)

            CASE(26_intk)   

                ! LOCATION
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(1) =  1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(2) =  1_intk
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%location(3) =  1_intk

                ! FACE COORDINATES
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(1) = maxx
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(2) = maxy
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_coord(3) = maxz

                ! NEIGBHOURS
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(1) = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(2) = particle_boundaries(igrid)%face_neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(3) = particle_boundaries(igrid)%face_neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(4) = particle_boundaries(igrid)%face_neighbours(18)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(5) = particle_boundaries(igrid)%face_neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(6) = particle_boundaries(igrid)%face_neighbours(14)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(7) = particle_boundaries(igrid)%face_neighbours(12)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_neighbours(8) = particle_boundaries(igrid)%face_neighbours(26)

                ! FACE NORMALS
                ! x-direction (1-4) 
                idir = 1_intk
                iface = 2_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 1) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 2) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 3) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(18)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 4) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! y-direction (5-8) 
                idir = 2_intk
                iface = 4_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 5) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(6)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 6) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 7) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(14)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 8) = particle_boundaries(jgrid)%face_normals(idir, iface)

                ! z-direction (9-12) 
                idir = 3_intk
                iface = 6_intk
                jgrid = igrid
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals( 9) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(2)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(10) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(4)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(11) = particle_boundaries(jgrid)%face_normals(idir, iface)

                jgrid = neighbours(12)
                particle_gcorner_boundaries((igrid - 1_intk) * 8_intk + icorn)%face_normals(12) = particle_boundaries(jgrid)%face_normals(idir, iface)

        END SELECT

    END SUBROUTINE set_gcorner_boundary

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