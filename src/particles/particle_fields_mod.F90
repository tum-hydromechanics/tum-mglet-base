MODULE particle_fields_mod

    USE grids_mod
    USE field_mod
    USE fields_mod

    USE particle_list_mod

    IMPLICIT NONE

CONTAINS

    SUBROUTINE init_particle_field()

        ! local variables
        INTEGER(intk), PARAMETER :: units_part(7) = [0, 0, 0, 0, 0, 0, 0]

        ! TODO: TIMER

        IF (.NOT. dwrite_npcfield) THEN
            RETURN
        END IF

        ! TODO: set dwrite to .FALSE.

        ! number per cell (NPC) field to count particles on cell
        ! with buffer so particle%ijkcell point to the right cell directly (neccessary?)
        CALL set_field("NPC", units=units_part, dwrite=.TRUE., buffers=.TRUE.)

        CALL update_particle_field()

    END SUBROUTINE init_particle_field

    SUBROUTINE update_particle_field()

        ! local_variables
        INTEGER(intk) :: igrid, ipart, ii, jj, kk, i_p, j_p, k_p
        REAL(realk) :: x_p, y_p, z_p
        TYPE(field_t), POINTER :: x_f, y_f, z_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz
        TYPE(field_t), POINTER :: npc_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: npc

        IF (.NOT. dwrite_npcfield) THEN
            RETURN
        END IF

        ! TODO: TIMER

        CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(npc_f, "NPC")

        npc_f%arr = 0.0
        !DO i = 1, nmygrids
        !    igrid = mygrids(i)
        !    CALL npc_f%get_ptr(npc, igrid)
        !    npc = 0.0
        !END DO

        DO ipart = 1, my_particle_list%ifinal

            IF (my_particle_list%particles(ipart)%state <= 0) CYCLE

            igrid = my_particle_list%particles(ipart)%igrid

            ! sanity check
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL x_f%get_ptr(x, igrid)
            CALL y_f%get_ptr(y, igrid)
            CALL z_f%get_ptr(z, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL npc_f%get_ptr(npc, igrid)

            i_p = my_particle_list%particles(ipart)%ijkcell(1)
            j_p = my_particle_list%particles(ipart)%ijkcell(2)
            k_p = my_particle_list%particles(ipart)%ijkcell(3)

            x_p = my_particle_list%particles(ipart)%x
            y_p = my_particle_list%particles(ipart)%y
            z_p = my_particle_list%particles(ipart)%z

            ! sanity checks
            IF (i_p < 1 .OR. i_p > ii .OR. j_p < 1 .OR. j_p > jj .OR. k_p < 1 .OR. k_p > kk) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING in update_particle_field: Particle ijkcell outside valid Range!"
                    CASE ("verbose")
                        WRITE(*, *) "WARNING in update_particle_field: Particle ijkcell outside valid Range!"
                    END SELECT
            END IF

            IF (x_p + EPSILON(x_p) < x(i_p) - dx(MAX(1, i_p - 1)) / 2 .OR. &
             x_p - EPSILON(x_p) > x(i_p) + dx(i_p) / 2 .OR. &
             y_p + EPSILON(y_p) < y(j_p) - dy(MAX(1, j_p - 1)) / 2 .OR. &
             y_p - EPSILON(y_p) > y(j_p) + dy(j_p) / 2 .OR. &
             z_p + EPSILON(z_p) < z(k_p) - dz(MAX(1, k_p - 1)) / 2 .OR. &
             z_p - EPSILON(z_p) > z(k_p) + dz(k_p) / 2) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) "WARNING in update_particle_field: Particle ijkcell does not coincide with Particle Coordinates."
                        WRITE(*, '()')
                    CASE ("verbose")
                        WRITE(*, *) "WARNING in update_particle_field: Particle ijkcell does not coincide with Particle Coordinates."
                        WRITE(*, '()')
                END SELECT
            END IF

            ! actual field computation
            npc(k_p, j_p, i_p) = npc(k_p, j_p, i_p) + 1

            SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CONTINUE
            CASE ("verbose")
                WRITE(*, *) "Particle registered on NPC field: igrid = ", igrid, "; cell: i=", i_p, " /j=", j_p, " /k=", k_p, "."
                WRITE(*, '()')
            END SELECT

        END DO

    END SUBROUTINE update_particle_field

END MODULE particle_fields_mod