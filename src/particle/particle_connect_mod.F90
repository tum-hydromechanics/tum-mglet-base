MODULE particle_connect_mod

USE core_mod
USE particlecore_mod

IMPLICIT NONE

CONTAINS

    SUBROUTINE init_pconnect()

    END SUBROUTINE init_pconnect

    SUBROUTINE pconnect(particle, nbrgrid, nbrproc)

        CLASS()baseparticle_t :: particle
        INTEGER :: nbrgrid, nbrproc

        CALL get_target_grid(my_particle_list%particles(i), nbrgrid)

        nbrproc = idprocofgrd(nbrgrid)

    END SUBROUTINE pconnect

    SUBROUTINE get_target_grid(particle, nbrgrid)

        !subroutine arguments
        CLASS(baseparticle_t), INTENT(in) :: particle
        INTEGER(intk), INTENT(out) :: nbrgrid

        !local variables
        INTEGER(intk) :: iface, neighbours(26)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz


        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, particle%igrid)

        IF (particle%x < minx) THEN !-------------------------------------------------------- low x

            IF (particle%y < miny) THEN !--------------------------------------------- low y, low x

                IF (particle%z < minz) THEN !---------------------------------- low z, low y, low x

                    iface = 19

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !---- mid z, low y, low x

                    iface = 7

                ELSEIF (maxz < particle%z) THEN !----------------------------- high z, low y, low x

                    iface = 20

                END IF

            ELSEIF (miny <= particle%y .AND. particle%y <= maxy) THEN !--------------- mid y, low x

                IF (particle%z < minz) THEN !---------------------------------- low z, mid y, low x

                    iface = 9

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !---- mid z, mid y, low x

                    iface = 1

                ELSEIF (maxz < particle%z) THEN !----------------------------- high z, mid y, low x

                    iface = 10

                END IF

            ELSEIF (maxy < particle%y) THEN !---------------------------------------- high y, low x

                IF (particle%z < minz) THEN !--------------------------------- low z, high y, low x

                    iface = 21

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !--- mid z, high y, low x

                    iface = 8

                ELSEIF (maxz < particle%z) THEN !---------------------------- high z, high y, low x

                    iface = 22

                END IF

            END IF

        ELSEIF (minx <= particle%x .AND. particle%x <= maxx) THEN !-------------------------- mid x

            IF (particle%y < miny) THEN !--------------------------------------------- low y, mid x

                IF (particle%z < minz) THEN !---------------------------------- low z, low y, mid x

                    iface = 15

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !---- mid z, low y, mid x

                    iface = 3

                ELSEIF (maxz < particle%z) THEN !----------------------------- high z, low y, mid x

                    iface = 16

                END IF

            ELSEIF (miny <= particle%y .AND. particle%y <= maxy) THEN !--------------- mid y, mid x

                IF (particle%z < minz) THEN !---------------------------------- low z, mid y, mid x

                    iface = 5

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !---- mid z, mid y, mid x

                    iface = 0

                ELSEIF (maxz < particle%z) THEN !----------------------------- high z, mid y, mid x

                    iface = 6

                END IF

            ELSEIF (maxy < particle%y) THEN !---------------------------------------- high y, mid x

                IF (particle%z < minz) THEN !--------------------------------- low z, high y, mid x

                    iface = 17

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !--- mid z, high y, mid x

                    iface = 4

                ELSEIF (maxz < particle%z) THEN !---------------------------- high z, high y, mid x

                    iface = 18

                END IF

            END IF

        ELSEIF (maxx < particle%x) THEN !--------------------------------------------------- high x

            IF (particle%y < miny) THEN !-------------------------------------------- low y, high x

                IF (particle%z < minz) THEN !--------------------------------- low z, low y, high x

                    iface = 23

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !--- mid z, low y, high x

                    iface = 11

                ELSEIF (maxz < particle%z) THEN !---------------------------- high z, low y, high x

                    iface = 24

                END IF

            ELSEIF (miny <= particle%y .AND. particle%y <= maxy) THEN !-------------- mid y, high x

                IF (particle%z < minz) THEN !--------------------------------- low z, mid y, high x

                    iface = 13

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !--- mid z, mid y, high x

                    iface = 2

                ELSEIF (maxz < particle%z) THEN !---------------------------- high z, mid y, high x

                    iface = 14

                END IF

            ELSEIF (maxy < particle%y) THEN !--------------------------------------- high y, high x

                IF (particle%z < minz) THEN !-------------------------------- low z, high y, high x

                    iface = 25

                ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !-- mid z, high y, high x

                    iface = 12

                ELSEIF (maxz < particle%z) THEN !--------------------------- high z, high y, high x

                    iface = 26

                END IF

            END IF

        END IF !-------------------------------------------------------------


        !itypbc = itypboconds(1, iface, particle%igrid)

        CALL get_neighbours(neighbours, particle%igrid)

        nbrgrid = neighbours(iface)

        !CALL get_nbrs(iface, neighbours, nbrgrid, nbrface)

    END SUBROUTINE get_target_grid

    SUBROUTINE get_exit_face(particle, pdx, pdy, pdz)

        ! subroutine arguments
        CLASS(baseparticle_t), INTENT(inout) :: particle
        REAL(realk), INTENT(in) :: pdx, pdy, pdz
        !INTEGER(intk), INTENT(out) :: sface_arr(3)

        !local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: lx, ly, lz, rx, ry, rz

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, particle%igrid)

        IF (pdx < 0) THEN

            lx = (minx - particle%x)
            rx = pdx / lx

        ELSEIF (0 < pdx) THEN

            lx = (maxx - particle%x)
            rx = pdx / lx

        ELSE

            rx = 0.0_realk

        END IF

        IF (pdy < 0) THEN

            ly = (miny - particle%y)
            ry = pdy / ly

        ELSEIF (0 < pdy) THEN

            ly = (maxy - particle%y)
            ry = pdy / ly

        ELSE

            ry = 0.0_realk

        END IF

        IF (pdz < 0) THEN

            lz = (minz - particle%z)
            rz = pdz / lz

        ELSEIF (0 < pdz) THEN

            lz = (maxz - particle%z)
            rz = pdz / lz

        ELSE

            rz = 0.0_realk

        END IF

        IF (rx <= 1.0_realk .AND. ry <= 1.0_realk .AND. rz <= 1.0_realk) THEN

            particle%facepath = 0

            RETURN

        END IF

        IF (pdx < 0 .AND. ry <= rx .AND. rz <= rx) THEN

            particle%facepath(1) = 1

            IF (pdy < 0 .AND. rz <= ry) THEN

                particle%facepath(2) = 3

                IF (pdz < 0) THEN

                    particle%facepath(3) = 5

                ELSEIF (0 < pdz) THEN

                    particle%facepath(3) = 6

                END IF

            ELSEIF (0 < pdy .AND. rz <= ry) THEN

                particle%facepath(2) = 4

                IF (pdz < 0) THEN

                    particle%facepath(3) = 5

                ELSEIF (0 < pdz) THEN

                    particle%facepath(3) = 6

                END IF

            ELSEIF (pdz < 0) THEN

                particle%facepath(2) = 5

                IF (pdy < 0) THEN

                    particle%facepath(3) = 3

                ELSEIF (0 < pdy) THEN

                    particle%facepath(3) = 4

                END IF

            ELSEIF (0 < pdz) THEN

                particle%facepath(2) = 6

                IF (pdy < 0) THEN

                    particle%facepath(3) = 3

                ELSEIF (0 < pdy) THEN

                    particle%facepath(3) = 4

                END IF

            END IF

        ELSEIF (0 < pdx .AND. ry <= rx .AND. rz <= rx) THEN

            particle%facepath(1) = 2

            IF (pdy < 0 .AND. rz <= ry) THEN

                particle%facepath(2) = 3

                IF (pdz < 0) THEN

                    particle%facepath(3) = 5

                ELSEIF (0 < pdz) THEN

                    particle%facepath(3) = 6

                END IF

            ELSEIF (0 < pdy .AND. rz <= ry) THEN

                particle%facepath(2) = 4

                IF (pdz < 0) THEN

                    particle%facepath(3) = 5

                ELSEIF (0 < pdz) THEN

                    particle%facepath(3) = 6

                END IF

            ELSEIF (pdz < 0) THEN

                particle%facepath(2) = 5

                IF (pdy < 0) THEN

                    particle%facepath(3) = 3

                ELSEIF (0 < pdy) THEN

                    particle%facepath(3) = 4

                END IF

            ELSEIF (0 < pdz) THEN

                particle%facepath(2) = 6

                IF (pdy < 0) THEN

                    particle%facepath(3) = 3

                ELSEIF (0 < pdy) THEN

                    particle%facepath(3) = 4

                END IF

            END IF

        ELSEIF (pdy < 0 .AND. rx <= ry .AND. rz <= ry) THEN

            particle%facepath(1) = 3

            IF (pdx < 0 .AND. rz <= rx) THEN

                particle%facepath(2) = 1

                IF (pdz < 0) THEN

                    particle%facepath(3) = 5

                ELSEIF (0 < pdz) THEN

                    particle%facepath(3) = 6

                END IF

            ELSEIF (0 < pdx .AND. rz <= rx) THEN

                particle%facepath(2) = 2

                IF (pdz < 0) THEN

                    particle%facepath(3) = 5

                ELSEIF (0 < pdz) THEN

                    particle%facepath(3) = 6

                END IF

            ELSEIF (pdz < 0) THEN

                particle%facepath(2) = 5

                IF (pdx < 0) THEN

                    particle%facepath(3) = 1

                ELSEIF (0 < pdx) THEN

                    particle%facepath(3) = 2

                END IF

            ELSEIF (0 < pdz) THEN

                particle%facepath(2) = 6

                IF (pdx < 0) THEN

                    particle%facepath(3) = 1

                ELSEIF (0 < pdx) THEN

                    particle%facepath(3) = 2

                END IF

            END IF

        ELSEIF (0 < pdy .AND. rx <= ry .AND. rz <= ry) THEN

            particle%facepath(1) = 4

            IF (pdx < 0 .AND. rz <= rx) THEN

                particle%facepath(2) = 1

                IF (pdz < 0) THEN

                    particle%facepath(3) = 5

                ELSEIF (0 < pdz) THEN

                    particle%facepath(3) = 6

                END IF

            ELSEIF (0 < pdx .AND. rz <= rx) THEN

                particle%facepath(2) = 2

                IF (pdz < 0) THEN

                    particle%facepath(3) = 5

                ELSEIF (0 < pdz) THEN

                    particle%facepath(3) = 6

                END IF

            ELSEIF (pdz < 0) THEN

                particle%facepath(2) = 5

                IF (pdx < 0) THEN

                    particle%facepath(3) = 1

                ELSEIF (0 < pdx) THEN

                    particle%facepath(3) = 2

                END IF

            ELSEIF (0 < pdz) THEN

                particle%facepath(2) = 6

                IF (pdx < 0) THEN

                    particle%facepath(3) = 1

                ELSEIF (0 < pdx) THEN

                    particle%facepath(3) = 2

                END IF

            END IF

        ELSEIF (pdz < 0 .AND. rx <= rz .AND. ry <= rz) THEN

            particle%facepath(1) = 5

            IF (pdx < 0 .AND. ry <= rx) THEN

                particle%facepath(2) = 1

                IF (pdy < 0) THEN

                    particle%facepath(3) = 3

                ELSEIF (0 < pdy) THEN

                    particle%facepath(3) = 4

                END IF

            ELSEIF (0 < pdx .AND. ry <= rx) THEN

                particle%facepath(2) = 2

                IF (pdy < 0) THEN

                    particle%facepath(3) = 3

                ELSEIF (0 < pdy) THEN

                    particle%facepath(3) = 4

                END IF

            ELSEIF (pdy < 0) THEN

                particle%facepath(2) = 3

                IF (pdx < 0) THEN

                    particle%facepath(3) = 1

                ELSEIF (0 < pdx) THEN

                    particle%facepath(3) = 2

                END IF

            ELSEIF (0 < pdy) THEN

                particle%facepath(2) = 4

                IF (pdx < 0) THEN

                    particle%facepath(3) = 1

                ELSEIF (0 < pdx) THEN

                    particle%facepath(3) = 2

                END IF

            END IF

        ELSEIF (0 < pdz .AND. rx <= rz .AND. ry <= rz) THEN

            particle%facepath(1) = 6

            IF (pdx < 0 .AND. ry <= rx) THEN

                particle%facepath(2) = 1

                IF (pdy < 0) THEN

                    particle%facepath(3) = 3

                ELSEIF (0 < pdy) THEN

                    particle%facepath(3) = 4

                END IF

            ELSEIF (0 < pdx .AND. ry <= rx) THEN

                particle%facepath(2) = 2

                IF (pdy < 0) THEN

                    particle%facepath(3) = 3

                ELSEIF (0 < pdy) THEN

                    particle%facepath(3) = 4

                END IF

            ELSEIF (pdy < 0) THEN

                particle%facepath(2) = 3

                IF (pdx < 0) THEN

                    particle%facepath(3) = 1

                ELSEIF (0 < pdx) THEN

                    particle%facepath(3) = 2

                END IF

            ELSEIF (0 < pdy) THEN

                particle%facepath(2) = 4

                IF (pdx < 0) THEN

                    particle%facepath(3) = 1

                ELSEIF (0 < pdx) THEN

                    particle%facepath(3) = 2

                END IF

            END IF

        END IF

    END SUBROUTINE get_exit_face

    SUBROUTINE is_true_neigbour()



    END SUBROUTINE is_true_neigbour

END MODULE particle_connect_mod