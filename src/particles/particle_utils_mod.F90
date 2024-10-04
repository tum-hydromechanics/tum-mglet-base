MODULE particle_utils_mod

    USE particle_core_mod

    IMPLICIT NONE

    INTERFACE get_exit_face
        MODULE PROCEDURE :: get_exit_face_c
        MODULE PROCEDURE :: get_exit_face_p
    END INTERFACE get_exit_face

    INTERFACE update_coordinates2
        MODULE PROCEDURE :: update_coordinates_c
        MODULE PROCEDURE :: update_coordinates_p
    END INTERFACE update_coordinates2

    CONTAINS

    SUBROUTINE get_exit_face_p(particle, dist, iface)

        TYPE(baseparticle_t), INTENT(in) :: particle
        REAL(realk), INTENT(out) :: dist
        INTEGER(intk), INTENT(out) :: iface

        CALL get_exit_face_c(particle%igrid, particle%x, particle%y, particle%z, dist, iface)

    END SUBROUTINE get_exit_face_p

    SUBROUTINE get_exit_face_c(igrid, x, y, z, dist, iface)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(out) :: dist
        INTEGER(intk), INTENT(out) :: iface

        ! local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: dist_x, dist_y, dist_z

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        IF (minx < x .AND. x < maxx .AND. &
            miny < y .AND. y < maxy .AND. &
            minz < z .AND. z < maxz) THEN
                iface = 0
                dist = 0.0
        ELSE
            ! checking the geometrical relation
            IF (x <= minx) THEN !-------------------------------------------------------- low x
                IF (y <= miny) THEN !--------------------------------------------- low y, low x
                    IF (z <= minz) THEN !---------------------------------- low z, low y, low x
                        iface = 19
                        dist_x = ABS(x - minx)
                        dist_y = ABS(y - miny)
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, low y, low x
                        iface = 7
                        dist_x = ABS(x - minx)
                        dist_y = ABS(y - miny)
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !----------------------------- high z, low y, low x
                        iface = 20
                        dist_x = ABS(x - minx)
                        dist_y = ABS(y - miny)
                        dist_z = ABS(z - maxz)
                    END IF
                ELSEIF (miny < y .AND. y < maxy) THEN !--------------------------- mid y, low x
                    IF (z <= minz) THEN !---------------------------------- low z, mid y, low x
                        iface = 9
                        dist_x = ABS(x - minx)
                        dist_y = 0.0
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, mid y, low x
                        iface = 1
                        dist_x = ABS(x - minx)
                        dist_y = 0.0
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !----------------------------- high z, mid y, low x
                        iface = 10
                        dist_x = ABS(x - minx)
                        dist_y = 0.0
                        dist_z = ABS(z - maxz)
                    END IF
                ELSEIF (maxy <= y) THEN !---------------------------------------- high y, low x
                    IF (z <= minz) THEN !--------------------------------- low z, high y, low x
                        iface = 21
                        dist_x = ABS(x - minx)
                        dist_y = ABS(y - maxy)
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, high y, low x
                        iface = 8
                        dist_x = ABS(x - minx)
                        dist_y = ABS(y - maxy)
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !---------------------------- high z, high y, low x
                        iface = 22
                        dist_x = ABS(x - minx)
                        dist_y = ABS(y - maxy)
                        dist_z = ABS(z - maxz)
                    END IF
                END IF
            ELSEIF (minx < x .AND. x < maxx) THEN !-------------------------------------- mid x
                IF (y <= miny) THEN !--------------------------------------------- low y, mid x
                    IF (z <= minz) THEN !---------------------------------- low z, low y, mid x
                        iface = 15
                        dist_x = 0.0
                        dist_y = ABS(y - miny)
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, low y, mid x
                        iface = 3
                        dist_x = 0.0
                        dist_y = ABS(y - miny)
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !----------------------------- high z, low y, mid x
                        iface = 16
                        dist_x = 0.0
                        dist_y = ABS(y - miny)
                        dist_z = ABS(z - maxz)
                    END IF
                ELSEIF (miny < y .AND. y < maxy) THEN !--------------------------- mid y, mid x
                    IF (z <= minz) THEN !---------------------------------- low z, mid y, mid x
                        iface = 5
                        dist_x = 0.0
                        dist_y = 0.0
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, mid y, mid x
                        iface = 0
                        dist_x = 0.0
                        dist_y = 0.0
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !----------------------------- high z, mid y, mid x
                        iface = 6
                        dist_x = 0.0
                        dist_y = 0.0
                        dist_z = ABS(z - maxz)
                    END IF
                ELSEIF (maxy <= y) THEN !---------------------------------------- high y, mid x
                    IF (z <= minz) THEN !--------------------------------- low z, high y, mid x
                        iface = 17
                        dist_x = 0.0
                        dist_y = ABS(y - maxy)
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, high y, mid x
                        iface = 4
                        dist_x = 0.0
                        dist_y = ABS(y - maxy)
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !---------------------------- high z, high y, mid x
                        iface = 18
                        dist_x = 0.0
                        dist_y = ABS(y - maxy)
                        dist_z = ABS(z - maxz)
                    END IF
                END IF
            ELSEIF (maxx <= x) THEN !--------------------------------------------------- high x
                IF (y <= miny) THEN !-------------------------------------------- low y, high x
                    IF (z <= minz) THEN !--------------------------------- low z, low y, high x
                        iface = 23
                        dist_x = ABS(x - maxx)
                        dist_y = ABS(y - miny)
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, low y, high x
                        iface = 11
                        dist_x = ABS(x - maxx)
                        dist_y = ABS(y - miny)
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !---------------------------- high z, low y, high x
                        iface = 24
                        dist_x = ABS(x - maxx)
                        dist_y = ABS(y - miny)
                        dist_z = ABS(z - maxz)
                    END IF
                ELSEIF (miny < y .AND. y < maxy) THEN !-------------------------- mid y, high x
                    IF (z <= minz) THEN !--------------------------------- low z, mid y, high x
                        iface = 13
                        dist_x = ABS(x - maxx)
                        dist_y = 0.0
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, mid y, high x
                        iface = 2
                        dist_x = ABS(x - maxx)
                        dist_y = 0.0
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !---------------------------- high z, mid y, high x
                        iface = 14
                        dist_x = ABS(x - maxx)
                        dist_y = 0.0
                        dist_z = ABS(z - maxz)
                    END IF
                ELSEIF (maxy <= y) THEN !--------------------------------------- high y, high x
                    IF (z <= minz) THEN !-------------------------------- low z, high y, high x
                        iface = 25
                        dist_x = ABS(x - maxx)
                        dist_y = ABS(y - maxy)
                        dist_z = ABS(z - minz)
                    ELSEIF (minz < z .AND. z < maxz) THEN !-------------- mid z, high y, high x
                        iface = 12
                        dist_x = ABS(x - maxx)
                        dist_y = ABS(y - maxy)
                        dist_z = 0.0
                    ELSEIF (maxz <= z) THEN !--------------------------- high z, high y, high x
                        iface = 26
                        dist_x = ABS(x - maxx)
                        dist_y = ABS(y - maxy)
                        dist_z = ABS(z - maxz)
                    END IF
                END IF
            END IF !-------------------------------------------------------------

            dist = SQRT(dist_x**2 + dist_y**2 + dist_z**2)

        END IF

    END SUBROUTINE get_exit_face_c

    SUBROUTINE update_coordinates_p(particle, destgrid, iface)

        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: destgrid, iface

        CALL update_coordinates_c(particle%igrid, destgrid, iface, particle%x, particle%y, particle%z)

    END SUBROUTINE update_coordinates_p

    SUBROUTINE update_coordinates_c(igrid, destgrid, iface, x, y, z)

        ! subroutine arguments

        INTEGER(intk), INTENT(in) :: igrid, destgrid, iface
        REAL(realk), INTENT(inout) :: x, y, z

        ! local variables
        REAL(realk) :: old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, &
         new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz
        LOGICAL :: passed_pb

        passed_pb = .FALSE.

        CALL get_bbox(old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, igrid)
        CALL get_bbox(new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz, destgrid)

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

    END SUBROUTINE update_coordinates_c


END MODULE particle_utils_mod