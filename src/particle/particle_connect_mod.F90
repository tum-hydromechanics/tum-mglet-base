MODULE particle_connect_mod

USE core_mod
USE particlecore_mod

IMPLICIT NONE

CONTAINS

    SUBROUTINE get_target_grid(igrid, x_new, y_new, z_new, iface, nbrgrid, nbrface)

        !subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(in) :: x_new, y_new, z_new
        INTEGER(intk), INTENT(out) :: iface, nbrgrid, nbrface

        !local variables
        INTEGER(intk) :: neighbours(26)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz


        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        IF (x_new < minx) THEN !--------------------------------------- low x

            IF (y_new < miny) THEN !---------------------------- low y, low x

                IF (z_new < minz) THEN !----------------- low z, low y, low x

                    iface = 19

                ELSEIF (minz <= z_new <= maxz) THEN !---- mid z, low y, low x

                    iface = 7

                ELSEIF (maxz < z_new) THEN !------------ high z, low y, low x

                    iface = 20

                END IF

            ELSEIF (miny <= y_new <= maxy) THEN !--------------- mid y, low x

                IF (z_new < minz) THEN !----------------- low z, mid y, low x

                    iface = 9

                ELSEIF (minz <= z_new <= maxz) THEN !---- mid z, mid y, low x

                    iface = 1

                ELSEIF (maxz < z_new) THEN !------------ high z, mid y, low x

                    iface = 10

                END IF

            ELSEIF (maxy < y_new) THEN !----------------------- high y, low x

                IF (z_new < minz) THEN !---------------- low z, high y, low x

                    iface = 21

                ELSEIF (minz <= z_new <= maxz) THEN !--- mid z, high y, low x

                    iface = 8

                ELSEIF (maxz < z_new) THEN !----------- high z, high y, low x

                    iface = 22

                END IF

            END IF

        ELSEIF (minx <= x_new <= maxx) THEN !-------------------------- mid x

            IF (y_new < miny) THEN !---------------------------- low y, mid x

                IF (z_new < minz) THEN !----------------- low z, low y, mid x

                    iface = 15

                ELSEIF (minz <= z_new <= maxz) THEN !---- mid z, low y, mid x

                    iface = 3

                ELSEIF (maxz < z_new) THEN !------------ high z, low y, mid x

                    iface = 16

                END IF

            ELSEIF (miny <= y_new <= maxy) THEN !--------------- mid y, mid x

                IF (z_new < minz) THEN !----------------- low z, mid y, mid x

                    iface = 5

                ELSEIF (minz <= z_new <= maxz) THEN !---- mid z, mid y, mid x

                    iface = 0

                ELSEIF (maxz < z_new) THEN !------------ high z, mid y, mid x

                    iface = 6

                END IF

            ELSEIF (maxy < y_new) THEN !----------------------- high y, mid x

                IF (z_new < minz) THEN !---------------- low z, high y, mid x

                    iface = 17

                ELSEIF (minz <= z_new <= maxz) THEN !--- mid z, high y, mid x

                    iface = 4

                ELSEIF (maxz < z_new) THEN !----------- high z, high y, mid x

                    iface = 18

                END IF

            END IF

        ELSEIF (maxx < x_new) THEN !---------------------------------- high x

            IF (y_new < miny) THEN !--------------------------- low y, high x

                IF (z_new < minz) THEN !---------------  low z, low y, high x

                    iface = 23

                ELSEIF (minz <= z_new <= maxz) THEN !--- mid z, low y, high x

                    iface = 11

                ELSEIF (maxz < z_new) THEN !----------- high z, low y, high x

                    iface = 24

                END IF

            ELSEIF (miny <= y_new <= maxy) THEN !-------------- mid y, high x

                IF (z_new < minz) THEN !---------------- low z, mid y, high x

                    iface = 13

                ELSEIF (minz <= z_new <= maxz) THEN !--- mid z, mid y, high x

                    iface = 2

                ELSEIF (maxz < z_new) THEN !----------- high z, mid y, high x

                    iface = 14

                END IF

            ELSEIF (maxy < y_new) THEN !---------------------- high y, high x

                IF (z_new < minz) THEN !--------------- low z, high y, high x

                    iface = 25

                ELSEIF (minz <= z_new <= maxz) THEN !-- mid z, high y, high x

                    iface = 12

                ELSEIF (maxz < z_new) THEN !---------- high z, high y, high x

                    iface = 26

                END IF

            END IF

        END IF !-------------------------------------------------------------


        !itypbc = itypboconds(1, iface, igrid)

        CALL get_neighbours(neighbours, igrid)

        nbrgrid = neighbours(iface)

        !CALL get_nbrs(iface, neighbours, nbrgrid, nbrface)

    END SUBROUTINE get_target_grid

    SUBROUTINE get_exit_face(pdx, pdy, pdz, lx, ly, lz)



    END SUBROUTINE get_exit_face

END MODULE particle_connect_mod