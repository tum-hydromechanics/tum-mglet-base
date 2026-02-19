MODULE particle_utils_mod

    USE particle_ofields_mod
    USE particle_basetype_mod

    IMPLICIT NONE

    INTEGER(intk), ALLOCATABLE :: facelist_utils(:)

    INTERFACE get_exit_face
        MODULE PROCEDURE :: get_exit_face_c
        MODULE PROCEDURE :: get_exit_face_p
    END INTERFACE get_exit_face

    INTERFACE update_coordinates
        MODULE PROCEDURE :: update_coordinates_c
        MODULE PROCEDURE :: update_coordinates_p
    END INTERFACE update_coordinates

    INTERFACE get_exit_face_target
        MODULE PROCEDURE :: get_exit_face_c_target
        MODULE PROCEDURE :: get_exit_face_p_target
    END INTERFACE get_exit_face_target

    INTERFACE update_coordinates_target
        MODULE PROCEDURE :: update_coordinates_c_target
        MODULE PROCEDURE :: update_coordinates_p_target
    END INTERFACE update_coordinates_target

    INTERFACE a_greater_b
        MODULE PROCEDURE i_greater_i
        MODULE PROCEDURE r_greater_r
    END INTERFACE a_greater_b

    INTERFACE a_greaterequal_b
        MODULE PROCEDURE i_greaterequal_i
        MODULE PROCEDURE r_greaterequal_r
    END INTERFACE a_greaterequal_b

#ifdef __GFORTRAN__
    !$omp declare target (facelist_utils)
#endif

    CONTAINS

    SUBROUTINE init_particle_utils()

        ALLOCATE(facelist_utils(27))

        facelist_utils(1)  = 19   ! low z,   low y,   low x  
        facelist_utils(2)  = 7    ! mid z,   low y,   low x
        facelist_utils(3)  = 20   ! high z,  low y,   low x 
        facelist_utils(4)  = 9    ! low z,   mid y,   low x
        facelist_utils(5)  = 1    ! mid z,   mid y,   low x
        facelist_utils(6)  = 10   ! high z,  mid y,   low x 
        facelist_utils(7)  = 21   ! low z,   high y,  low x 
        facelist_utils(8)  = 8    ! mid z,   high y,  low x
        facelist_utils(9)  = 22   ! high z,  high y,  low x 
        facelist_utils(10) = 15   ! low z,   low y,   mid x 
        facelist_utils(11) = 3    ! mid z,   low y,   mid x
        facelist_utils(12) = 16   ! high z,  low y,   mid x 
        facelist_utils(13) = 5    ! low z,   mid y,   mid x
        facelist_utils(14) = 0    ! mid z,   mid y,   mid x
        facelist_utils(15) = 6    ! high z,  mid y,   mid x
        facelist_utils(16) = 17   ! low z,   high y,  mid x 
        facelist_utils(17) = 4    ! mid z,   high y,  mid x
        facelist_utils(18) = 18   ! high z,  high y,  mid x 
        facelist_utils(19) = 23   ! low z,   low y,   high x 
        facelist_utils(20) = 11   ! mid z,   low y,   high x 
        facelist_utils(21) = 24   ! high z,  low y,   high x 
        facelist_utils(22) = 13   ! low z,   mid y,   high x 
        facelist_utils(23) = 2    ! mid z,   mid y,   high x
        facelist_utils(24) = 14   ! high z,  mid y,   high x 
        facelist_utils(25) = 25   ! low z,   high y,  high x 
        facelist_utils(26) = 12   ! mid z,   high y,  high x 
        facelist_utils(27) = 26   ! high z,  high y,  high x 

        !$omp target enter data map(to: facelist_utils)

    END SUBROUTINE init_particle_utils


    SUBROUTINE finish_particle_utils()

        !$omp target exit data map(delete: facelist_utils)

    END SUBROUTINE finish_particle_utils


    ! get the face (iface) that indicates the correct neigbouring grid (of the outdated particle%igrid)
    ! that the particle is actually on after its displacement
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
                RETURN
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

    ! get the face (iface) that indicates the correct neigbouring grid (of the outdated particle%igrid)
    ! that the particle is actually on after its displacement
    SUBROUTINE get_exit_face_p_target(bbox, particle, dist, iface)

        !$omp declare target
        
        REAL(realk), INTENT(in) :: bbox(6)
        TYPE(baseparticle_t), INTENT(in) :: particle
        REAL(realk), INTENT(out) :: dist
        INTEGER(intk), INTENT(out) :: iface

        CALL get_exit_face_c_target(bbox, particle%x, particle%y, particle%z, dist, iface)

    END SUBROUTINE get_exit_face_p_target

SUBROUTINE get_exit_face_c_target(bbox, x, y, z, dist, iface)

        !$omp declare target

        ! subroutine arguments
        REAL(realk), INTENT(in) :: bbox(6)
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(out) :: dist
        INTEGER(intk), INTENT(out) :: iface

        ! local variables
        REAL(realk) :: dist_x, dist_y, dist_z

        IF (bbox(1) < x .AND. x < bbox(2) .AND. &
            bbox(3) < y .AND. y < bbox(4) .AND. &
            bbox(5) < z .AND. z < bbox(6)) THEN
                iface = 0
                dist = 0.0
                RETURN
        ELSE
            ! checking the geometrical relation
            IF (x <= bbox(1)) THEN !-------------------------------------------------------- low x
                IF (y <= bbox(3)) THEN !--------------------------------------------- low y, low x
                    IF (z <= bbox(5)) THEN !---------------------------------- low z, low y, low x
                        iface = 19
                        dist_x = ABS(x - bbox(1))
                        dist_y = ABS(y - bbox(3))
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !---------------- mid z, low y, low x
                        iface = 7
                        dist_x = ABS(x - bbox(1))
                        dist_y = ABS(y - bbox(3))
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !----------------------------- high z, low y, low x
                        iface = 20
                        dist_x = ABS(x - bbox(1))
                        dist_y = ABS(y - bbox(3))
                        dist_z = ABS(z - bbox(6))
                    END IF
                ELSEIF (bbox(3) < y .AND. y < bbox(4)) THEN !--------------------------- mid y, low x
                    IF (z <= bbox(5)) THEN !---------------------------------- low z, mid y, low x
                        iface = 9
                        dist_x = ABS(x - bbox(1))
                        dist_y = 0.0
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !---------------- mid z, mid y, low x
                        iface = 1
                        dist_x = ABS(x - bbox(1))
                        dist_y = 0.0
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !----------------------------- high z, mid y, low x
                        iface = 10
                        dist_x = ABS(x - bbox(1))
                        dist_y = 0.0
                        dist_z = ABS(z - bbox(6))
                    END IF
                ELSEIF (bbox(4) <= y) THEN !---------------------------------------- high y, low x
                    IF (z <= bbox(5)) THEN !--------------------------------- low z, high y, low x
                        iface = 21
                        dist_x = ABS(x - bbox(1))
                        dist_y = ABS(y - bbox(4))
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !--------------- mid z, high y, low x
                        iface = 8
                        dist_x = ABS(x - bbox(1))
                        dist_y = ABS(y - bbox(4))
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !---------------------------- high z, high y, low x
                        iface = 22
                        dist_x = ABS(x - bbox(1))
                        dist_y = ABS(y - bbox(4))
                        dist_z = ABS(z - bbox(6))
                    END IF
                END IF
            ELSEIF (bbox(1) < x .AND. x < bbox(2)) THEN !-------------------------------------- mid x
                IF (y <= bbox(3)) THEN !--------------------------------------------- low y, mid x
                    IF (z <= bbox(5)) THEN !---------------------------------- low z, low y, mid x
                        iface = 15
                        dist_x = 0.0
                        dist_y = ABS(y - bbox(3))
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !---------------- mid z, low y, mid x
                        iface = 3
                        dist_x = 0.0
                        dist_y = ABS(y - bbox(3))
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !----------------------------- high z, low y, mid x
                        iface = 16
                        dist_x = 0.0
                        dist_y = ABS(y - bbox(3))
                        dist_z = ABS(z - bbox(6))
                    END IF
                ELSEIF (bbox(3) < y .AND. y < bbox(4)) THEN !--------------------------- mid y, mid x
                    IF (z <= bbox(5)) THEN !---------------------------------- low z, mid y, mid x
                        iface = 5
                        dist_x = 0.0
                        dist_y = 0.0
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !---------------- mid z, mid y, mid x
                        iface = 0
                        dist_x = 0.0
                        dist_y = 0.0
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !----------------------------- high z, mid y, mid x
                        iface = 6
                        dist_x = 0.0
                        dist_y = 0.0
                        dist_z = ABS(z - bbox(6))
                    END IF
                ELSEIF (bbox(4) <= y) THEN !---------------------------------------- high y, mid x
                    IF (z <= bbox(5)) THEN !--------------------------------- low z, high y, mid x
                        iface = 17
                        dist_x = 0.0
                        dist_y = ABS(y - bbox(4))
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !--------------- mid z, high y, mid x
                        iface = 4
                        dist_x = 0.0
                        dist_y = ABS(y - bbox(4))
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !---------------------------- high z, high y, mid x
                        iface = 18
                        dist_x = 0.0
                        dist_y = ABS(y - bbox(4))
                        dist_z = ABS(z - bbox(6))
                    END IF
                END IF
            ELSEIF (bbox(2) <= x) THEN !--------------------------------------------------- high x
                IF (y <= bbox(3)) THEN !-------------------------------------------- low y, high x
                    IF (z <= bbox(5)) THEN !--------------------------------- low z, low y, high x
                        iface = 23
                        dist_x = ABS(x - bbox(2))
                        dist_y = ABS(y - bbox(3))
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !--------------- mid z, low y, high x
                        iface = 11
                        dist_x = ABS(x - bbox(2))
                        dist_y = ABS(y - bbox(3))
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !---------------------------- high z, low y, high x
                        iface = 24
                        dist_x = ABS(x - bbox(2))
                        dist_y = ABS(y - bbox(3))
                        dist_z = ABS(z - bbox(6))
                    END IF
                ELSEIF (bbox(3) < y .AND. y < bbox(4)) THEN !-------------------------- mid y, high x
                    IF (z <= bbox(5)) THEN !--------------------------------- low z, mid y, high x
                        iface = 13
                        dist_x = ABS(x - bbox(2))
                        dist_y = 0.0
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !--------------- mid z, mid y, high x
                        iface = 2
                        dist_x = ABS(x - bbox(2))
                        dist_y = 0.0
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !---------------------------- high z, mid y, high x
                        iface = 14
                        dist_x = ABS(x - bbox(2))
                        dist_y = 0.0
                        dist_z = ABS(z - bbox(6))
                    END IF
                ELSEIF (bbox(4) <= y) THEN !--------------------------------------- high y, high x
                    IF (z <= bbox(5)) THEN !-------------------------------- low z, high y, high x
                        iface = 25
                        dist_x = ABS(x - bbox(2))
                        dist_y = ABS(y - bbox(4))
                        dist_z = ABS(z - bbox(5))
                    ELSEIF (bbox(5) < z .AND. z < bbox(6)) THEN !-------------- mid z, high y, high x
                        iface = 12
                        dist_x = ABS(x - bbox(2))
                        dist_y = ABS(y - bbox(4))
                        dist_z = 0.0
                    ELSEIF (bbox(6) <= z) THEN !--------------------------- high z, high y, high x
                        iface = 26
                        dist_x = ABS(x - bbox(2))
                        dist_y = ABS(y - bbox(4))
                        dist_z = ABS(z - bbox(6))
                    END IF
                END IF
            END IF !-------------------------------------------------------------

            dist = SQRT(dist_x**2 + dist_y**2 + dist_z**2)

        END IF

    END SUBROUTINE get_exit_face_c_target

    
    SUBROUTINE get_exit_face_c_target2(bbox, x, y, z, dist, iface)

        !$omp declare target

        ! subroutine arguments
        REAL(realk), INTENT(in) :: bbox(6)
        REAL(realk), INTENT(in) :: x, y, z
        REAL(realk), INTENT(out) :: dist
        INTEGER(intk), INTENT(out) :: iface

        ! local variables
        INTEGER(intk) :: i, j, k

        i = 1 - NINT(a_greaterequal_b(bbox(1), x)) + NINT(a_greaterequal_b(x, bbox(2)))
        j = 1 - NINT(a_greaterequal_b(bbox(3), y)) + NINT(a_greaterequal_b(y, bbox(4)))
        k = 1 - NINT(a_greaterequal_b(bbox(5), z)) + NINT(a_greaterequal_b(z, bbox(6)))

        IF (i == 1 .AND. j == 1 .AND. k == 1) THEN
            iface = 0
            dist  = 0.0_realk
            RETURN
        END IF

        iface = facelist_utils(i*9 + j*3 + k + 1)
        dist = SQRT(a_greaterequal_b(bbox(1), x) * (x - bbox(1))**2 + a_greaterequal_b(x, bbox(2)) * (x - bbox(2))**2 + &
                    a_greaterequal_b(bbox(3), y) * (y - bbox(3))**2 + a_greaterequal_b(y, bbox(4)) * (y - bbox(4))**2 + &
                    a_greaterequal_b(bbox(5), z) * (z - bbox(5))**2 + a_greaterequal_b(z, bbox(6)) * (z - bbox(6))**2 )


    END SUBROUTINE get_exit_face_c_target2


    ! update particle coordinates such if particle crossed a periodic boundary
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
                    ELSEIF (new_maxx <= x) THEN
                        x = new_minx + ABS(x - old_maxx)
                        passed_pb = .TRUE.
                    END IF
                END IF

                IF (reflect(2) == 0) THEN
                    IF (y <= new_miny) THEN
                        y = new_maxy - ABS(y - old_miny)
                        passed_pb = .TRUE.
                    ELSEIF (new_maxy <= y) THEN
                        y = new_miny + ABS(y - old_maxy)
                        passed_pb = .TRUE.
                    END IF
                END IF

                IF (reflect(3) == 0) THEN
                    IF (z <= new_minz) THEN
                        z = new_maxz - ABS(z - old_minz)
                        passed_pb = .TRUE.
                    ELSEIF (new_maxz <= z) THEN
                        z = new_minz + ABS(z - old_maxz)
                        passed_pb = .TRUE.
                    END IF
                END IF

            ELSE

                IF (x < new_minx) THEN
                    x = new_maxx - ABS(x - old_minx)
                    passed_pb = .TRUE.
                ELSEIF (new_maxx < x) THEN
                    x = new_minx + ABS(x - old_maxx)
                    passed_pb = .TRUE.
                END IF

                IF (y < new_miny) THEN
                    y = new_maxy - ABS(y - old_miny)
                    passed_pb = .TRUE.
                ELSEIF (new_maxy < y) THEN
                    y = new_miny + ABS(y - old_maxy)
                    passed_pb = .TRUE.
                END IF

                IF (z < new_minz) THEN
                    z = new_maxz - ABS(z - old_minz)
                    passed_pb = .TRUE.
                ELSEIF (new_maxz < z) THEN
                    z = new_minz + ABS(z - old_maxz)
                    passed_pb = .TRUE.
                END IF

            END IF

        ELSE ! this case is for the particle exchange module

            IF (x < new_minx) THEN
                x = new_maxx - ABS(x - old_minx)
                passed_pb = .TRUE.
            ELSEIF (new_maxx < x) THEN
                x = new_minx + ABS(x - old_maxx)
                passed_pb = .TRUE.
            END IF

            IF (y < new_miny) THEN
                y = new_maxy - ABS(y - old_miny)
                passed_pb = .TRUE.
            ELSEIF (new_maxy < y) THEN
                y = new_miny + ABS(y - old_maxy)
                passed_pb = .TRUE.
            END IF

            IF (z < new_minz) THEN
                z = new_maxz - ABS(z - old_minz)
                passed_pb = .TRUE.
            ELSEIF (new_maxz < z) THEN
                z = new_minz + ABS(z - old_maxz)
                passed_pb = .TRUE.
            END IF

        END IF

    END SUBROUTINE update_coordinates_c


    ! update particle coordinates such if particle crossed a periodic boundary
    SUBROUTINE update_coordinates_p_target(particle, destgrid, iface, old_bbox)

        !$omp declare target

        TYPE(baseparticle_t), INTENT(inout) :: particle
        INTEGER(intk), INTENT(in) :: destgrid, iface
        REAL(realk), INTENT(in) :: old_bbox(6)

        CALL update_coordinates_c_target(particle%igrid, destgrid, iface, particle%x, particle%y, particle%z, old_bbox)

    END SUBROUTINE update_coordinates_p_target


    SUBROUTINE update_coordinates_c_target(igrid, destgrid, iface, x, y, z, old_bbox, reflect)

        !$omp declare target

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, destgrid, iface
        REAL(realk), INTENT(inout) :: x, y, z
        REAL(realk), INTENT(in) :: old_bbox(6)
        INTEGER(intk), INTENT(in), OPTIONAL :: reflect(3)

        ! local variables
        REAL(realk) :: new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz
         
        LOGICAL :: passed_pb

        passed_pb = .FALSE.

        IF (iface == 0) THEN
            RETURN
        END IF

        CALL get_bbox_target(new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz, destgrid)

        IF (PRESENT(reflect)) THEN ! this case is for the particle boundaries module

            IF (igrid == destgrid) THEN

                IF (reflect(1) == 0) THEN
                    IF (x <= new_minx) THEN
                        x = new_maxx - ABS(x - old_bbox(1))
                        passed_pb = .TRUE.
                    ELSEIF (new_maxx <= x) THEN
                        x = new_minx + ABS(x - old_bbox(2))
                        passed_pb = .TRUE.
                    END IF
                END IF

                IF (reflect(2) == 0) THEN
                    IF (y <= new_miny) THEN
                        y = new_maxy - ABS(y - old_bbox(3))
                        passed_pb = .TRUE.
                    ELSEIF (new_maxy <= y) THEN
                        y = new_miny + ABS(y - old_bbox(4))
                        passed_pb = .TRUE.
                    END IF
                END IF

                IF (reflect(3) == 0) THEN
                    IF (z <= new_minz) THEN
                        z = new_maxz - ABS(z - old_bbox(5))
                        passed_pb = .TRUE.
                    ELSEIF (new_maxz <= z) THEN
                        z = new_minz + ABS(z - old_bbox(6))
                        passed_pb = .TRUE.
                    END IF
                END IF

            ELSE

                IF (x < new_minx) THEN
                    x = new_maxx - ABS(x - old_bbox(1))
                    passed_pb = .TRUE.
                ELSEIF (new_maxx < x) THEN
                    x = new_minx + ABS(x - old_bbox(2))
                    passed_pb = .TRUE.
                END IF

                IF (y < new_miny) THEN
                    y = new_maxy - ABS(y - old_bbox(3))
                    passed_pb = .TRUE.
                ELSEIF (new_maxy < y) THEN
                    y = new_miny + ABS(y - old_bbox(4))
                    passed_pb = .TRUE.
                END IF

                IF (z < new_minz) THEN
                    z = new_maxz - ABS(z - old_bbox(5))
                    passed_pb = .TRUE.
                ELSEIF (new_maxz < z) THEN
                    z = new_minz + ABS(z - old_bbox(6))
                    passed_pb = .TRUE.
                END IF

            END IF

        ELSE ! this case is for the particle exchange module

            IF (x < new_minx) THEN
                x = new_maxx - ABS(x - old_bbox(1))
                passed_pb = .TRUE.
            ELSEIF (new_maxx < x) THEN
                x = new_minx + ABS(x - old_bbox(2))
                passed_pb = .TRUE.
            END IF

            IF (y < new_miny) THEN
                y = new_maxy - ABS(y - old_bbox(3))
                passed_pb = .TRUE.
            ELSEIF (new_maxy < y) THEN
                y = new_miny + ABS(y - old_bbox(4))
                passed_pb = .TRUE.
            END IF

            IF (z < new_minz) THEN
                z = new_maxz - ABS(z - old_bbox(5))
                passed_pb = .TRUE.
            ELSEIF (new_maxz < z) THEN
                z = new_minz + ABS(z - old_bbox(6))
                passed_pb = .TRUE.
            END IF

        END IF

    END SUBROUTINE update_coordinates_c_target


    ! function that returns true if given coordinates lie within given grid and false otherwise
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

    SUBROUTINE sort_conns_unique(list, sort_idx, check_redundancy, forbidden_val)
        ! Input array to be sorted
        INTEGER(int32), INTENT(inout) :: list(:,:)
        INTEGER(intk), INTENT(in) :: sort_idx
        LOGICAL, INTENT(in) :: check_redundancy 
        INTEGER(intk), INTENT(in), OPTIONAL :: forbidden_val 

        INTEGER(intk) :: i, j

        ! Temporary storage
        INTEGER(int32), ALLOCATABLE :: temp(:)

        ALLOCATE(temp(SIZE(list, 1)))

        ! Sort by sending processor number (field 2) (rising order)
        DO i = 2, SIZE(list, 2)
            j = i - 1
            temp(:) = list(:,i)
            DO WHILE (j >= 1)
                IF (list(sort_idx, j) > temp(sort_idx)) THEN
                    list(:,j+1) = list(:,j)
                    j = j - 1
                ELSE
                    EXIT
                END IF
            END DO
            list(:,j+1) = temp(:)
        END DO

        ! Check for redundant entries
        IF (check_redundancy .AND. PRESENT(forbidden_val)) THEN
            DO i = 2, SIZE(list, 2)
                IF ( list(sort_idx, i) == list(sort_idx, i-1) ) THEN
                    WRITE(*,*) 'Redundant listing: ', list(sort_idx, i)
                    CALL errr(__FILE__, __LINE__)
                END IF
                IF (list(sort_idx, i) == forbidden_val) THEN
                    WRITE(*,*) 'Forbidden Value listed:  ', list(sort_idx, i)
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO
        ELSEIF (check_redundancy) THEN
            DO i = 2, SIZE(list, 2)
                IF ( list(sort_idx, i) == list(sort_idx, i-1) ) THEN
                    WRITE(*,*) 'Redundant listing: ', list(sort_idx, i)
                    CALL errr(__FILE__, __LINE__)
                END IF
            END DO
        END IF

    END SUBROUTINE sort_conns_unique

    PURE SUBROUTINE conditional_update_ii(old, new, large_to_old, large_to_new)

        !$omp declare target 
        
        !subroutine arguments
        INTEGER(intk), INTENT(inout) :: old
        INTEGER(intk), INTENT(in) :: new
        INTEGER(intk), INTENT(in) :: large_to_old, large_to_new

        old = old * MAX(0_intk, SIGN(1_intk, large_to_old - large_to_new)) &
            + new * MAX(0_intk, SIGN(1_intk, large_to_new - large_to_old))

    END SUBROUTINE conditional_update_ii


    PURE SUBROUTINE conditional_update_ir(old, new, large_to_old, large_to_new)

        !$omp declare target 
        
        !subroutine arguments
        INTEGER(intk), INTENT(inout) :: old
        INTEGER(intk), INTENT(in) :: new
        REAL(realk), INTENT(in) :: large_to_old, large_to_new

        old = old * INT(MAX(0.0_realk, SIGN(1.0_realk, large_to_old - large_to_new))) &
            + new * INT(MAX(0.0_realk, SIGN(1.0_realk, large_to_new - large_to_old)))

    END SUBROUTINE conditional_update_ir


    PURE SUBROUTINE conditional_update_ri(old, new, large_to_old, large_to_new)

        !$omp declare target 
        
        !subroutine arguments
        REAL(realk), INTENT(inout) :: old
        REAL(realk), INTENT(in) :: new
        INTEGER(realk), INTENT(in) :: large_to_old, large_to_new

        old = old * REAL(MAX(0_intk, SIGN(1_intk, large_to_old - large_to_new))) &
            + new * REAL(MAX(0_intk, SIGN(1_intk, large_to_new - large_to_old)))

    END SUBROUTINE conditional_update_ri


    SUBROUTINE conditional_update_rr(old, new, large_to_old, large_to_new)

        !$omp declare target 
        
        !subroutine arguments
        REAL(realk), INTENT(inout) :: old
        REAL(realk), INTENT(in) :: new
        REAL(realk), INTENT(in) :: large_to_old, large_to_new

        old = old * r_greater_r(large_to_old, large_to_new) &
            + new * r_greater_r(large_to_new, large_to_old)

    END SUBROUTINE conditional_update_rr


    REAL(realk) FUNCTION r_greater_r(a, b) result(factor)

        !$omp declare target 

        REAL(realk), INTENT(in) :: a, b

        factor = MAX(0.0_realk, SIGN(1.0_realk, a - b))

    END FUNCTION r_greater_r


    REAL(realk) FUNCTION r_greaterequal_r(a, b) result(factor)

        !$omp declare target 

        REAL(realk), INTENT(in) :: a, b

        factor = MAX(0.0_realk, SIGN(1.0_realk, a - b + EPSILON(a)))

    END FUNCTION r_greaterequal_r


    REAL(realk) FUNCTION i_greater_i(a, b) result(factor)

        !$omp declare target 

        INTEGER(intk), INTENT(in) :: a, b

        factor = REAL(MAX(0_intk, SIGN(1_intk, a - b)), realk)

    END FUNCTION i_greater_i


    REAL(realk) FUNCTION i_greaterequal_i(a, b) result(factor)

        !$omp declare target 

        INTEGER(intk), INTENT(in) :: a, b

        factor = REAL(MAX(0_intk, SIGN(1_intk, a - b + 1_intk)), realk)

    END FUNCTION i_greaterequal_i



END MODULE particle_utils_mod