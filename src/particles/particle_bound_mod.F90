MODULE particle_bound_mod
    USE precision_mod, ONLY: realk, intk
    USE core_mod
    USE particle_list_mod

    ! particle boundary types integer identifiers
    ! connect  => 1
    ! periodic => 2
    ! reflect  => 3


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

        INTEGER(intk), ALLOCATABLE :: front(:)
        INTEGER(intk), ALLOCATABLE :: back(:)
        INTEGER(intk), ALLOCATABLE :: right(:)
        INTEGER(intk), ALLOCATABLE :: left(:)
        INTEGER(intk), ALLOCATABLE :: bottom(:)
        INTEGER(intk), ALLOCATABLE :: top(:)

        TYPE(obstacle_list_t), ALLOCATABLE :: obstacle_lists(:) ! every grid has one list! (not an exact analog to the particle list)

    END TYPE particle_boundaries_t

    TYPE(particle_boundaries_t) :: particle_boundaries

    CONTAINS

    SUBROUTINE init_particle_boundaries()

        ! local variables
        INTEGER(intk) :: i, igrid, iface, ibocd, nbocd
        INTEGER(intk) :: neighbours(26)

        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, n_minx, n_maxx, n_miny, n_maxy, n_minz, n_maxz

        CHARACTER(len=8) :: ctyp

        ALLOCATE(particle_boundaries%front(ngrid))
        ALLOCATE(particle_boundaries%back(ngrid))
        ALLOCATE(particle_boundaries%right(ngrid))
        ALLOCATE(particle_boundaries%left(ngrid))
        ALLOCATE(particle_boundaries%bottom(ngrid))
        ALLOCATE(particle_boundaries%top(ngrid))
        ALLOCATE(particle_boundaries%obstacle_lists(ngrid))

        particle_boundaries%front = 0
        particle_boundaries%back = 0
        particle_boundaries%right = 0
        particle_boundaries%left = 0
        particle_boundaries%bottom = 0
        particle_boundaries%top = 0



        DO i = 1, ngrid
            igrid = i ! (?) <------------------------------------------
            DO iface = 1, 6
                ibocd = 3
                    CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                    WRITE(*, *) "CTYP: ", ctyp
                    SELECT CASE(iface)
                    CASE(1)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%front(igrid) = 3
                        ELSEIF (ctyp == 'PER') THEN
                            particle_boundaries%front(igrid) = 2
                        ELSEIF (ctyp == 'CON') THEN
                            CALL get_neighbours(neighbours, igrid)
                            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
                            CALL get_bbox(n_minx, n_maxx, n_miny, n_maxy, n_minz, n_maxz, neighbours(iface))
                            IF (maxx <= n_minx .OR. igrid == neighbours(iface)) THEN
                                particle_boundaries%front(igrid) = 2
                            ELSE
                                particle_boundaries%front(igrid) = 1
                            END IF
                        END IF
                    CASE(2)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%back(igrid) = 3
                        ELSEIF (ctyp == 'PER') THEN
                            particle_boundaries%back(igrid) = 2
                        ELSEIF (ctyp == 'CON') THEN
                            CALL get_neighbours(neighbours, igrid)
                            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
                            CALL get_bbox(n_minx, n_maxx, n_miny, n_maxy, n_minz, n_maxz, neighbours(iface))
                            IF (n_maxx <= minx .OR. igrid == neighbours(iface)) THEN
                                particle_boundaries%back(igrid) = 2
                            ELSE
                                particle_boundaries%back(igrid) = 1
                            END IF
                        END IF
                    CASE(3)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%right(igrid) = 3
                        ELSEIF (ctyp == 'PER') THEN
                            particle_boundaries%right(igrid) = 2
                        ELSEIF (ctyp == 'CON') THEN
                            CALL get_neighbours(neighbours, igrid)
                            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
                            CALL get_bbox(n_minx, n_maxx, n_miny, n_maxy, n_minz, n_maxz, neighbours(iface))
                            IF (maxy <= n_miny .OR. igrid == neighbours(iface)) THEN
                                particle_boundaries%right(igrid) = 2
                            ELSE
                                particle_boundaries%right(igrid) = 1
                            END IF
                        END IF
                    CASE(4)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%left(igrid) = 3
                        ELSEIF (ctyp == 'PER') THEN
                            particle_boundaries%left(igrid) = 2
                        ELSEIF (ctyp == 'CON') THEN
                            CALL get_neighbours(neighbours, igrid)
                            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
                            CALL get_bbox(n_minx, n_maxx, n_miny, n_maxy, n_minz, n_maxz, neighbours(iface))
                            IF (n_maxy <= miny .OR. igrid == neighbours(iface)) THEN
                                particle_boundaries%left(igrid) = 2
                            ELSE
                                particle_boundaries%left(igrid) = 1
                            END IF
                        END IF
                    CASE(5)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%bottom(igrid) = 3
                        ELSEIF (ctyp == 'PER') THEN
                            particle_boundaries%bottom(igrid) = 2
                        ELSEIF (ctyp == 'CON') THEN
                            CALL get_neighbours(neighbours, igrid)
                            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
                            CALL get_bbox(n_minx, n_maxx, n_miny, n_maxy, n_minz, n_maxz, neighbours(iface))
                            IF (maxz <= n_minz .OR. igrid == neighbours(iface)) THEN
                                particle_boundaries%bottom(igrid) = 2
                            ELSE
                                particle_boundaries%bottom(igrid) = 1
                            END IF
                        END IF
                    CASE(6)
                        IF (ctyp == 'REF') THEN
                            particle_boundaries%top(igrid) = 3
                        ELSEIF (ctyp == 'PER') THEN
                            particle_boundaries%top(igrid) = 2
                        ELSEIF (ctyp == 'CON') THEN
                            CALL get_neighbours(neighbours, igrid)
                            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
                            CALL get_bbox(n_minx, n_maxx, n_miny, n_maxy, n_minz, n_maxz, neighbours(iface))
                            IF (n_maxz <= minz .OR. igrid == neighbours(iface)) THEN
                                particle_boundaries%top(igrid) = 2
                            ELSE
                                particle_boundaries%top(igrid) = 1
                            END IF
                        END IF
                    END SELECT
            END DO
        END DO

        DO i = 1, ngrid

            WRITE(*, *) "Boundaries of Grid:   ", i
            WRITE(*, *) "FRONT:                ", particle_boundaries%front(i)
            WRITE(*, *) "BACK:                 ", particle_boundaries%back(i)
            WRITE(*, *) "RIGHT:                ", particle_boundaries%right(i)
            WRITE(*, *) "LEFT:                 ", particle_boundaries%left(i)
            WRITE(*, *) "BOTTOM:               ", particle_boundaries%bottom(i)
            WRITE(*, *) "TOP:                  ", particle_boundaries%top(i)

        END DO

    END SUBROUTINE init_particle_boundaries

    !-----------------------------------

    SUBROUTINE finish_particle_boundaries()

        DEALLOCATE(particle_boundaries%front)
        DEALLOCATE(particle_boundaries%back)
        DEALLOCATE(particle_boundaries%right)
        DEALLOCATE(particle_boundaries%left)
        DEALLOCATE(particle_boundaries%bottom)
        DEALLOCATE(particle_boundaries%top)
        DEALLOCATE(particle_boundaries%obstacle_lists)

    END SUBROUTINE finish_particle_boundaries

    !-----------------------------------

    SUBROUTINE move_particle(particle, dx, dy, dz)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(inout) :: particle
        REAL(realk), INTENT(inout) :: dx, dy, dz

        !local variables
        INTEGER(intk) :: newgrid, oldgrid, iface, grid_bc, destproc
        INTEGER(intk) :: neighbours(26)
        REAL(realk) :: x, y, z
        REAL(realk) :: dx_to_here, dy_to_here, dz_to_here, dx_from_here, dy_from_here, dz_from_here
        REAL(realk) :: epsilon
        REAL(realk) :: n1, n2, n3

        epsilon = SQRT(dx**(2) + dy**(2) + dz**(2)) / 10**3

        oldgrid = particle%igrid
        newgrid = particle%igrid
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
                CALL print_particle_status(particle)
        END SELECT

        DO WHILE (epsilon < ABS(dx_from_here) .OR. epsilon < ABS(dy_from_here) .OR. epsilon < ABS(dz_from_here))

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, *) ""
                    WRITE(*, *) "Way to go:"
                    WRITE(*, *) "dx/dy/dz:", dx_from_here, dy_from_here, dz_from_here
            END SELECT

            CALL move_to_boundary(oldgrid, x, y, z, &
            dx_from_here, dy_from_here, dz_from_here, dx_to_here, dy_to_here, dz_to_here, iface)

            SELECT CASE(iface)
                CASE(0)
                    grid_bc = 0
                CASE(1)
                    grid_bc = particle_boundaries%front(oldgrid)
                CASE(2)
                    grid_bc = particle_boundaries%back(oldgrid)
                CASE(3)
                    grid_bc = particle_boundaries%right(oldgrid)
                CASE(4)
                    grid_bc = particle_boundaries%left(oldgrid)
                CASE(5)
                    grid_bc = particle_boundaries%bottom(oldgrid)
                CASE(6)
                    grid_bc = particle_boundaries%top(oldgrid)
            END SELECT


            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, *) ""
                    WRITE(*, *) "Boundary:"
                    WRITE(*, *) "iface:", iface
            END SELECT

            IF (grid_bc == 0) THEN


                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*, *) "Same Grid."
                        WRITE(*, *) ""
                END SELECT

                newgrid = oldgrid

                particle%state = MAX(particle%state, 2)

            ELSEIF (grid_bc == 1) THEN

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*, *) "Connect."
                        WRITE(*, *) ""
                END SELECT

                CALL get_neighbours(neighbours, oldgrid)
                newgrid = neighbours(iface)

                particle%state = MAX(particle%state, 3)

            ELSEIF (grid_bc == 2) THEN

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*, *) "Periodic."
                        WRITE(*, *) ""
                END SELECT

                CALL get_neighbours(neighbours, oldgrid)
                newgrid = neighbours(iface)

                CALL apply_periodic_boundary(x, y, z, oldgrid, newgrid, iface)

                particle%state = MAX(particle%state, 3)

            ELSEIF (grid_bc == 3) THEN

                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        CONTINUE
                    CASE ("verbose")
                        WRITE(*, *) "Reflect."
                END SELECT

                newgrid = oldgrid

                CALL get_normal_vector(iface, n1, n2, n3)
                CALL apply_reflect_boundary(dx_from_here, dy_from_here, dz_from_here, n1, n2, n3)

                particle%state = MAX(particle%state, 2)

            ELSE

                ! TODO: ERROR/WARNING

            END IF

            oldgrid = newgrid

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, *) "New igrid:", newgrid
                    WRITE(*, *) "x/y/z:", x, y, z
            END SELECT

        END DO

        destproc = idprocofgrd(newgrid)

        IF ( destproc > numprocs .OR. destproc < 0 ) THEN
            WRITE(*,*) 'Ill-addressed particle to proc', destproc
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (destproc /= myid) THEN
            particle%state = 4
        END IF

        particle%iproc = destproc
        particle%igrid = newgrid
        particle%x = x
        particle%y = y
        particle%z = z

        IF (particle%state == 2) THEN
            CALL update_particle_cell(particle)
            particle%state = 1
        ELSEIF (particle%state == 3) THEN
            CALL set_particle_cell(particle)
            particle%state = 1
        END IF

        dx = 0
        dy = 0
        dz = 0

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                CALL print_particle_status(particle)
                WRITE(*, *) "MOVE PARTICLE: ----------------------END-----------------------"
        END SELECT

    END SUBROUTINE move_particle

    !-----------------------------------

    ! this subroutine only considers grids on the same level
    SUBROUTINE move_to_boundary(igrid, x, y, z, dx, dy, dz, dx_to_b, dy_to_b, dz_to_b, iface)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(inout) :: dx, dy, dz
        REAL(realk), INTENT(out) :: dx_to_b, dy_to_b, dz_to_b
        INTEGER(intk), INTENT(out) :: iface

        !local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: lx, ly, lz, rx, ry, rz

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        IF (dx < 0) THEN
            lx = (minx - x)
            rx = dx / lx
        ELSEIF (0 < dx) THEN
            lx = (maxx - x)
            rx = dx / lx
        ELSE
            rx = 0.0_realk
        END IF

        IF (dy < 0) THEN
            ly = (miny - y)
            ry = dy / ly
        ELSEIF (0 < dy) THEN
            ly = (maxy - y)
            ry = dy / ly
        ELSE
            ry = 0.0_realk
        END IF

        IF (dz < 0) THEN
            lz = (minz - z)
            rz = dz / lz
        ELSEIF (0 < dz) THEN
            lz = (maxz - z)
            rz = dz / lz
        ELSE
            rz = 0.0_realk
        END IF

        IF (rx < 1.0_realk .AND. ry < 1.0_realk .AND. rz < 1.0_realk) THEN

            iface = 0

            x = x + dx
            y = y + dy
            z = z + dz
            dx = 0
            dy = 0
            dz = 0

            RETURN

        END IF

        IF (dx < 0 .AND. ry <= rx .AND. rz <= rx) THEN

            iface = 1

            x = x + lx
            y = y + (lx * dy/dx)
            z = z + (lx * dz/dx)
            dx_to_b = lx
            dy_to_b = (lx * dy/dx)
            dz_to_b = (lx * dz/dx)
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dx .AND. ry <= rx .AND. rz <= rx) THEN

            iface = 2

            x = x + lx
            y = y + (lx * dy/dx)
            z = z + (lx * dz/dx)
            dx_to_b = lx
            dy_to_b = (lx * dy/dx)
            dz_to_b = (lx * dz/dx)
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (dy < 0 .AND. rx < ry .AND. rz <= ry) THEN

            iface = 3

            x = x + (ly * dx/dy)
            y = y + ly
            z = z + (ly * dz/dy)
            dx_to_b = (ly * dx/dy)
            dy_to_b = ly
            dz_to_b = (ly * dz/dy)
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dy .AND. rx < ry .AND. rz <= ry) THEN

            iface = 4

            x = x + (ly * dx/dy)
            y = y + ly
            z = z + (lz * dz/dy)
            dx_to_b = (ly * dx/dy)
            dy_to_b = ly
            dz_to_b = (ly * dz/dy)
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b


        ELSEIF (dz < 0 .AND. rx < rz .AND. ry < rz) THEN

            iface = 5

            x = x + (lz * dx/dz)
            y = y + (lz * dy/dz)
            z = z + lz
            dx_to_b = (lz * dx/dz)
            dy_to_b = (lz * dy/dz)
            dz_to_b = lz
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        ELSEIF (0 < dz .AND. rx < rz .AND. ry < rz) THEN

            iface = 6

            x = x + (lz * dx/dz)
            y = y + (lz * dy/dz)
            z = z + lz
            dx_to_b = (lz * dx/dz)
            dy_to_b = (lz * dy/dz)
            dz_to_b = lz
            dx = dx - dx_to_b
            dy = dy - dy_to_b
            dz = dz - dz_to_b

        END IF

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

    SUBROUTINE apply_reflect_boundary(dx, dy, dz, n1, n2, n3)

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

    END SUBROUTINE apply_reflect_boundary

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

END MODULE particle_bound_mod