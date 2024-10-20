MODULE particle_utils_mod

    USE particle_core_mod

    IMPLICIT NONE

    INTERFACE get_exit_face
        MODULE PROCEDURE :: get_exit_face_c
        MODULE PROCEDURE :: get_exit_face_p
    END INTERFACE get_exit_face

    INTERFACE update_coordinates
        MODULE PROCEDURE :: update_coordinates_c
        MODULE PROCEDURE :: update_coordinates_p
    END INTERFACE update_coordinates

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

    SUBROUTINE update_coordinates_c(igrid, destgrid, iface, x, y, z, reflect)

        ! subroutine arguments

        INTEGER(intk), INTENT(in) :: igrid, destgrid, iface
        REAL(realk), INTENT(inout) :: x, y, z
        INTEGER(intk), INTENT(in), OPTIONAL :: reflect(3)

        ! local variables
        REAL(realk) :: old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, &
         new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz
        LOGICAL :: passed_pb

        passed_pb = .FALSE.

        IF (iface == 0) THEN
            RETURN
        END IF

        CALL get_bbox(old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, igrid)
        CALL get_bbox(new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz, destgrid)

        IF (PRESENT(reflect)) THEN ! this case is for the particle boundaries module

            IF (igrid == destgrid) THEN

                IF (reflect(1) == 0) THEN
                    IF (x <= new_minx) THEN
                        x = new_maxx - ABS(x - old_minx)
                        passed_pb = .TRUE.
                    END IF

                    IF (new_maxx <= x) THEN
                        x = new_minx + ABS(x - old_maxx)
                        passed_pb = .TRUE.
                    END IF
                END IF

                IF (reflect(2) == 0) THEN
                    IF (y <= new_miny) THEN
                        y = new_maxy - ABS(y - old_miny)
                        passed_pb = .TRUE.
                    END IF

                    IF (new_maxy <= y) THEN
                        y = new_miny + ABS(y - old_maxy)
                        passed_pb = .TRUE.
                    END IF
                END IF

                IF (reflect(3) == 0) THEN
                    IF (z <= new_minz) THEN
                        z = new_maxz - ABS(z - old_minz)
                        passed_pb = .TRUE.
                    END IF

                    IF (new_maxz <= z) THEN
                        z = new_minz + ABS(z - old_maxz)
                        passed_pb = .TRUE.
                    END IF
                END IF

            ELSE

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

            END IF

        ELSE ! this case is for the particle exchange module

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

        END IF

    END SUBROUTINE update_coordinates_c

    FUNCTION is_inside_grid(igrid, x, y, z) result(res)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(in) :: x, y, z

        ! local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        ! return value
        LOGICAL :: res

        res = .TRUE.

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        IF (x < minx) THEN
            res = .FALSE.
        END IF

        IF (x > maxx) THEN
            res = .FALSE.
        END IF

        IF (y < miny) THEN
            res = .FALSE.
        END IF

        IF (y > maxy) THEN
            res = .FALSE.
        END IF

        IF (z < minz) THEN
            res = .FALSE.
        END IF

        IF (z > maxz) THEN
            res = .FALSE.
        END IF

    END FUNCTION is_inside_grid

        !IF (minx <= particle%x .AND. particle%x <= maxx .AND. &
        !    miny <= particle%y .AND. particle%y <= maxy .AND. &
        !    minz <= particle%z .AND. particle%z <= maxz) THEN
        !            iface = 0
        !ELSE
        !    ! checking the geometrical relation
        !    IF (particle%x < minx) THEN !-------------------------------------------------------- low x
        !        IF (particle%y < miny) THEN !--------------------------------------------- low y, low x
        !            IF (particle%z < minz) THEN !---------------------------------- low z, low y, low x
        !                iface = 19
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !---- mid z, low y, low x
        !                iface = 7
        !            ELSEIF (maxz < particle%z) THEN !----------------------------- high z, low y, low x
        !                iface = 20
        !            END IF
        !        ELSEIF (miny <= particle%y .AND. particle%y <= maxy) THEN !--------------- mid y, low x
        !            IF (particle%z < minz) THEN !---------------------------------- low z, mid y, low x
        !                iface = 9
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !---- mid z, mid y, low x
        !                iface = 1
        !            ELSEIF (maxz < particle%z) THEN !----------------------------- high z, mid y, low x
        !                iface = 10
        !            END IF
        !        ELSEIF (maxy < particle%y) THEN !---------------------------------------- high y, low x
        !            IF (particle%z < minz) THEN !--------------------------------- low z, high y, low x
        !                iface = 21
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !--- mid z, high y, low x
        !                iface = 8
        !            ELSEIF (maxz < particle%z) THEN !---------------------------- high z, high y, low x
        !                iface = 22
        !            END IF
        !        END IF
        !    ELSEIF (minx <= particle%x .AND. particle%x <= maxx) THEN !-------------------------- mid x
        !        IF (particle%y < miny) THEN !--------------------------------------------- low y, mid x
        !            IF (particle%z < minz) THEN !---------------------------------- low z, low y, mid x
        !                iface = 15
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !---- mid z, low y, mid x
        !                iface = 3
        !            ELSEIF (maxz < particle%z) THEN !----------------------------- high z, low y, mid x
        !                iface = 16
        !            END IF
        !        ELSEIF (miny <= particle%y .AND. particle%y <= maxy) THEN !--------------- mid y, mid x
        !            IF (particle%z < minz) THEN !---------------------------------- low z, mid y, mid x
        !                iface = 5
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !---- mid z, mid y, mid x
        !                iface = 0
        !            ELSEIF (maxz < particle%z) THEN !----------------------------- high z, mid y, mid x
        !                iface = 6
        !            END IF
        !        ELSEIF (maxy < particle%y) THEN !---------------------------------------- high y, mid x
        !            IF (particle%z < minz) THEN !--------------------------------- low z, high y, mid x
        !                iface = 17
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !--- mid z, high y, mid x
        !                iface = 4
        !            ELSEIF (maxz < particle%z) THEN !---------------------------- high z, high y, mid x
        !                iface = 18
        !            END IF
        !        END IF
        !    ELSEIF (maxx < particle%x) THEN !--------------------------------------------------- high x
        !        IF (particle%y < miny) THEN !-------------------------------------------- low y, high x
        !            IF (particle%z < minz) THEN !--------------------------------- low z, low y, high x
        !                iface = 23
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !--- mid z, low y, high x
        !                iface = 11
        !            ELSEIF (maxz < particle%z) THEN !---------------------------- high z, low y, high x
        !                iface = 24
        !            END IF
        !        ELSEIF (miny <= particle%y .AND. particle%y <= maxy) THEN !-------------- mid y, high x
        !            IF (particle%z < minz) THEN !--------------------------------- low z, mid y, high x
        !                iface = 13
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !--- mid z, mid y, high x
        !                iface = 2
        !            ELSEIF (maxz < particle%z) THEN !---------------------------- high z, mid y, high x
        !                iface = 14
        !            END IF
        !        ELSEIF (maxy < particle%y) THEN !--------------------------------------- high y, high x
        !            IF (particle%z < minz) THEN !-------------------------------- low z, high y, high x
        !                iface = 25
        !            ELSEIF (minz <= particle%z .AND. particle%z <= maxz) THEN !-- mid z, high y, high x
        !                iface = 12
        !            ELSEIF (maxz < particle%z) THEN !--------------------------- high z, high y, high x
        !                iface = 26
        !            END IF
        !        END IF
        !    END IF !-------------------------------------------------------------
        !END IF

    !SUBROUTINE get_current_face(igrid, x, y, z, iface)
!
    !    ! subroutine arguments
    !    INTEGER(intk), INTENT(in) :: igrid
    !    REAL(realk), INTENT(in) :: x, y, z
    !    INTEGER(intk), INTENT(out) :: iface
!
    !    ! local variables
    !    REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
!
    !    CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
!
    !    IF (x == minx) THEN !-------------------------------------------------------- low x
    !        IF (y == miny) THEN !--------------------------------------------- low y, low x
    !            IF (z == minz) THEN !---------------------------------- low z, low y, low x
    !                iface = 19
    !            ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, low y, low x
    !                iface = 7
    !            ELSEIF (z == maxz) THEN !----------------------------- high z, low y, low x
    !                iface = 20
    !            END IF
    !        ELSEIF (miny < y .AND. y < maxy) THEN !--------------------------- mid y, low x
    !            IF (z == minz) THEN !---------------------------------- low z, mid y, low x
    !                iface = 9
    !            ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, mid y, low x
    !                iface = 1
    !            ELSEIF (z == maxz) THEN !----------------------------- high z, mid y, low x
    !                iface = 10
    !            END IF
    !        ELSEIF (y == maxy) THEN !---------------------------------------- high y, low x
    !            IF (z == minz) THEN !--------------------------------- low z, high y, low x
    !                iface = 21
    !            ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, high y, low x
    !                iface = 8
    !            ELSEIF (z == maxz) THEN !---------------------------- high z, high y, low x
    !                iface = 22
    !            END IF
    !        END IF
    !    ELSEIF (minx < x .AND. x < maxx) THEN !-------------------------------------- mid x
    !        IF (y == miny) THEN !--------------------------------------------- low y, mid x
    !            IF (z == minz) THEN !---------------------------------- low z, low y, mid x
    !                iface = 15
    !            ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, low y, mid x
    !                iface = 3
    !            ELSEIF (z == maxz) THEN !----------------------------- high z, low y, mid x
    !                iface = 16
    !            END IF
    !        ELSEIF (miny < y .AND. y < maxy) THEN !--------------------------- mid y, mid x
    !            IF (z == minz) THEN !---------------------------------- low z, mid y, mid x
    !                iface = 5
    !            ELSEIF (minz < z .AND. z < maxz) THEN !---------------- mid z, mid y, mid x
    !                iface = 0
    !            ELSEIF (z == maxz) THEN !----------------------------- high z, mid y, mid x
    !                iface = 6
    !            END IF
    !        ELSEIF (y == maxy) THEN !---------------------------------------- high y, mid x
    !            IF (z == minz) THEN !--------------------------------- low z, high y, mid x
    !                iface = 17
    !            ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, high y, mid x
    !                iface = 4
    !            ELSEIF (z == maxz) THEN !---------------------------- high z, high y, mid x
    !                iface = 18
    !            END IF
    !        END IF
    !    ELSEIF (x == maxx) THEN !--------------------------------------------------- high x
    !        IF (y == miny) THEN !-------------------------------------------- low y, high x
    !            IF (z == minz) THEN !--------------------------------- low z, low y, high x
    !                iface = 23
    !            ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, low y, high x
    !                iface = 11
    !            ELSEIF (z == maxz) THEN !---------------------------- high z, low y, high x
    !                iface = 24
    !            END IF
    !        ELSEIF (miny < y .AND. y < maxy) THEN !-------------------------- mid y, high x
    !            IF (z == minz) THEN !--------------------------------- low z, mid y, high x
    !                iface = 13
    !            ELSEIF (minz < z .AND. z < maxz) THEN !--------------- mid z, mid y, high x
    !                iface = 2
    !            ELSEIF (z == maxz) THEN !---------------------------- high z, mid y, high x
    !                iface = 14
    !            END IF
    !        ELSEIF (y == maxy) THEN !--------------------------------------- high y, high x
    !            IF (z == minz) THEN !-------------------------------- low z, high y, high x
    !                iface = 25
    !            ELSEIF (minz < z .AND. z < maxz) THEN !-------------- mid z, high y, high x
    !                iface = 12
    !            ELSEIF (z == maxz) THEN !--------------------------- high z, high y, high x
    !                iface = 26
    !            END IF
    !        END IF
    !    END IF
!
    !END SUBROUTINE get_current_face

END MODULE particle_utils_mod