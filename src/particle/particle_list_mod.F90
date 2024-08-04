! this file is for particle lists

MODULE particle_list_mod

    !===================================

    USE particlecore_mod
    USE particle_dict_mod

    IMPLICIT NONE

    !-----------------------------------

    TYPE :: particle_list_t

        INTEGER(intk) :: iproc           ! REMOVE (obsolete) ?

        INTEGER(intk) :: max_np          ! max number of particles of this process/list
        INTEGER(intk) :: active_np = 0   ! number of active particles of this process/list
        INTEGER(intk) :: ifinal          ! index of last entry of the list which holds an active particle

        TYPE(baseparticle_t), ALLOCATABLE :: particles(:)
        !LOGICAL, ALLOCATABLE :: particle_stored(:)
        ! each logical value reflects whether a particle is stored in the list
        ! at the respective index. Is this a feasable and good way to keep track
        ! of particle storage (especially as is_active in particle_t carries the same information?

        CONTAINS

        PROCEDURE :: defragment, resize

    END TYPE particle_list_t

    !-----------------------------------

    ! LOGICAL :: dread_particles = .FALSE.

    TYPE(particle_list_t) :: my_particle_list ! rather declare in init_particle_list? NO

    PUBLIC :: my_particle_list

    !===================================

CONTAINS

    SUBROUTINE defragment( this )

        ! SIMON: Here just as an idea...

        ! Subroutine arguments
        CLASS(particle_list_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(intk) :: i, j, ifin, n
        LOGICAL :: cont

        ! local copy
        ifin = this%ifinal
        cont = .TRUE.

        IF ( this%active_np == ifin ) THEN
            ! all slots are occupied
            RETURN
        END IF

        DO i = 1, ifin

            ! search for empty slot
            IF ( this%particles(i)%is_active /= 1 ) THEN

                ! search from the end of list and find particle to fill in
                DO j = this%ifinal, 1, -1
                    ! finished if positions before "i" are considered
                    IF ( j < (i+1) ) THEN
                        cont = .FALSE.
                        EXIT
                    END IF
                    ! fill empty slot with last valid particle
                    IF ( this%particles(j)%is_active == 1 ) THEN
                        this%particles(i) = this%particles(j)
                        this%particles(j)%is_active = 0
                        this%ifinal = j - 1
                        EXIT
                    END IF
                END DO

                ! empty slot could not be filled
                IF ( .NOT. cont ) THEN
                    this%ifinal = i - 1
                    EXIT
                END IF

            END IF

        END DO

        ! Debug check
        DO i = 1, this%ifinal
            IF ( this%particles(i)%is_active /= 1 ) THEN
                WRITE(*,*) "defragmented list contains inactive particle at ", i
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        ! Some output (to be silenced later)
        IF ( ifin /= this%ifinal ) THEN
            n = ifin - this%ifinal
            WRITE(*,*) " - defragmentation: ", n, " slots cleared"
        END IF

    END SUBROUTINE defragment



    SUBROUTINE resize( this, n )

        ! SIMON: Here just as an idea...

        ! Subroutine arguments
        CLASS(particle_list_t), INTENT(inout) :: this
        INTEGER(intk), INTENT(in) :: n

        ! Local variables
        TYPE(baseparticle_t), ALLOCATABLE :: new_particles(:)

        IF ( n < this%ifinal ) THEN
            WRITE(*,*) "new list size is too small"
            CALL errr(__FILE__, __LINE__)
        END IF

        ALLOCATE( new_particles(n) )
        new_particles(1:this%ifinal) = this%particles(1:this%ifinal)

        CALL MOVE_ALLOC( new_particles, this%particles )

        ! checking if temporary array is gone
        IF ( ALLOCATED(new_particles) ) THEN
            WRITE(*,*) "move_alloc failed to deallocate"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! checking if desired array is there
        IF ( ALLOCATED(this%particles) ) THEN
            IF ( SIZE(this%particles, dim=1, kind=intk) /= n ) THEN
                WRITE(*,*) "move_alloc allocated wrong size"
                CALL errr(__FILE__, __LINE__)
            END IF
        ELSE
            WRITE(*,*) "move_alloc deallocated wrong array"
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE resize



    SUBROUTINE init_particle_list()

        ! local variables
        INTEGER(intk) :: i, j, global_np
        INTEGER(intk), ALLOCATABLE :: ipart_arr(:), p_igrid_arr(:)
        REAL(realk), ALLOCATABLE :: x(:), y(:), z(:)

        my_particle_list%iproc = myid
        my_particle_list%max_np = plist_len

        ALLOCATE(my_particle_list%particles(my_particle_list%max_np))

        IF (dread_particles) THEN

            CALL read_particles(dread_particles, my_particle_list%max_np, global_np, ipart_arr, p_igrid_arr, x, y, z)

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*,*) ' '
                        WRITE(*, '("Initializing ", I0, " Particle(s):")') global_np
                    CASE ("verbose")
                        WRITE(*,*) ' '
                        WRITE(*, '("Initializing ", I0, " Particle(s):")') global_np
                END SELECT
            END IF

            DO i = 1, global_np

                DO j = 1, nmygrids

                    IF (mygrids(j) == p_igrid_arr(i)) THEN

                        my_particle_list%active_np = my_particle_list%active_np + 1_intk

                        ! CALL my_particle_list%particles(my_particle_list%active_np)%init(ipart = ipart_arr(i), &
                        !  x = x(i), y = y(i), z = z(i), igrid = p_igrid_arr(i))

                        ! CALL my_particle_list%particles(i)%print_status()

                        CALL set_particle( my_particle_list%particles(my_particle_list%active_np), &
                        ipart = ipart_arr(i), x = x(i), y = y(i), z = z(i), &
                        igrid = p_igrid_arr(i))

                        CALL print_particle_status( my_particle_list%particles(i) )


                    END IF

                END DO

            END DO

            my_particle_list%ifinal = my_particle_list%active_np

        END IF

        IF (.NOT. dread_particles) THEN

            my_particle_list%active_np = init_npart

            ALLOCATE(ipart_arr(my_particle_list%active_np))
            ALLOCATE(p_igrid_arr(my_particle_list%active_np))
            ALLOCATE(x(my_particle_list%active_np))
            ALLOCATE(y(my_particle_list%active_np))
            ALLOCATE(z(my_particle_list%active_np))

            CALL dist_ipart(my_particle_list%active_np, ipart_arr)

            CALL dist_part(my_particle_list%active_np, p_igrid_arr, x, y, z)

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*,*) ' '
                        WRITE(*, '("Initializing ", I0, " Particle(s):")') my_particle_list%active_np
                    CASE ("verbose")
                        WRITE(*,*) ' '
                        WRITE(*, '("Initializing ", I0, " Particle(s):")') my_particle_list%active_np
                END SELECT
            END IF

            my_particle_list%ifinal = my_particle_list%active_np

            DO i = 1, my_particle_list%active_np

                CALL set_particle( my_particle_list%particles(i), &
                    ipart = ipart_arr(i), x = x(i), y = y(i), z = z(i), &
                    igrid = p_igrid_arr(i))

                CALL print_particle_status( my_particle_list%particles(i) )

            END DO

        END IF

        IF (myid == 0) THEN
            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    WRITE(*, '("Initialization of Particle(s) finished successfully.")')
                    WRITE(*,*) ' '
                CASE ("verbose")
                    WRITE(*, '("Initialization of Particle(s) finished successfully.")')
                    WRITE(*,*) ' '
            END SELECT
        END IF

    END SUBROUTINE init_particle_list

    !-----------------------------------

    SUBROUTINE dist_ipart(npart, ipart_arr) ! this routine is supposed to hand out a list of unique particle ids (ipart) to every process ! ONLY NON MPI RUNS FOR NOW

        ! subroutine arguments
        INTEGER(intk), INTENT(IN) :: npart
        INTEGER(intk), INTENT(out) :: ipart_arr(npart)

        ! local variables
        INTEGER(intk) :: i, counter

        counter = 1
        DO i = myid * npart + 1, myid * npart + npart
            ipart_arr(counter) = i
            counter = counter + 1_intk
        END DO

        ! MPI barrier

    END SUBROUTINE dist_ipart

    !-----------------------------------

    SUBROUTINE dist_part(npart, p_igrid_arr, x, y, z)

        ! subroutine arguments...
        INTEGER(intk), INTENT(in) :: npart
        INTEGER(intk), INTENT(inout) :: p_igrid_arr(npart)
        REAL(realk), INTENT(inout) :: x(npart), y(npart), z(npart)

        ! local variables...
        INTEGER(intk) :: i, j, igrid
        REAL(realk) :: myvolume, volume_fractions(nmygrids), grid_rn
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        ! compute combined volume of all grids this process owns (for now assumning ther is only one level)...

        myvolume = 0.0_realk

        DO i = 1, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            myvolume = myvolume + (maxx - minx) * (maxy - miny) * (maxz - minz)

        END DO

        ! compute fraction of the combined volume of each grids volume (for now assumning ther is only one level)...

        igrid = mygrids(1)
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

        volume_fractions(1) = (maxx - minx) * (maxy - miny) * (maxz - minz) / myvolume

        volume_fractions(1) = MAX(0.0_realk, volume_fractions(1))
        volume_fractions(1) = MIN(1.0_realk, volume_fractions(1))

        DO i = 2, nmygrids

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            volume_fractions(i) = volume_fractions(i-1) + ((maxx - minx) * (maxy - miny) * (maxz - minz) / myvolume)

            volume_fractions(i) = MAX(0.0_realk, volume_fractions(i))
            volume_fractions(i) = MIN(1.0_realk, volume_fractions(i))

        END DO

        ! distribute particles along all grids this process owns (uniformely distibuted)

        DO j = 1, npart

            i = 1_intk

            CALL RANDOM_NUMBER(grid_rn)

            DO WHILE (grid_rn > volume_fractions(i))

                i = i + 1_intk

            END DO

            igrid = mygrids(i)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            p_igrid_arr(j) = igrid

            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(x(j))
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(y(j))
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(z(j))

            x(j) = minx + x(j) * (maxx - minx)
            y(j) = miny + y(j) * (maxy - miny)
            z(j) = minz + z(j) * (maxz - minz)

        END DO

    END SUBROUTINE dist_part

    !===================================

END MODULE
