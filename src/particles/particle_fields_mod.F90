MODULE particle_fields_mod

    USE grids_mod
    USE field_mod
    USE fields_mod

    USE particle_list_mod

    IMPLICIT NONE

CONTAINS

    SUBROUTINE init_particle_field()

        ! local variables
        INTEGER(intk), PARAMETER :: units_part(7) = [0, -3, 0, 0, 0, 0, 0]

        ! TODO: switch to using intfield_t

        ! number per cell (NPC) field to count particles on cell
        ! with buffer so particle%ijkcell point to the right cell directly (neccessary?)
        CALL set_field("NPC", units=units_part, buffers=.TRUE.)

        CALL update_particle_field()

    END SUBROUTINE init_particle_field

    SUBROUTINE update_particle_field()

        ! local_variables
        INTEGER(intk) :: igrid, ipart, i, j, k
        TYPE(field_t), POINTER :: npc_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: npc

        CALL get_field(npc_f, "NPC")

        npc_f%arr = 0.0

        DO ipart = 1, my_particle_list%ifinal

            igrid = my_particle_list%particles(ipart)%igrid

            CALL npc_f%get_ptr(npc, igrid)

            i = my_particle_list%particles(ipart)%ijkcell(1)
            j = my_particle_list%particles(ipart)%ijkcell(2)
            k = my_particle_list%particles(ipart)%ijkcell(3)

            npc(k, j, i) = npc(k, j, i) + 1

        END DO

    END SUBROUTINE update_backup_fields

END MODULE particle_fields_mod