MODULE particle_connect_mod

    USE, INTRINSIC :: ISO_C_BINDING
    USE MPI_f08
    USE core_mod
    USE comms_mod, ONLY: myid
    USE particle_list_mod
    USE particlecore_mod

    IMPLICIT NONE (type, external)

    PRIVATE

    ! Maximum number of connections on one single process, either
    ! outgoing or incomming, on any single grid level
    INTEGER(intk) :: maxConns

    ! Particle type (not a class, as otherwise polymorphism implied)
    TYPE(baseparticle_t), ALLOCATABLE :: sendBufParticle(:)
    TYPE(baseparticle_t), ALLOCATABLE :: recvBufParticle(:)

    ! Sizes of the buffers
    INTEGER(intk) :: sizeSendBuf
    INTEGER(intk) :: sizeRecvBuf

    ! MPI type for the particle
    TYPE(MPI_Datatype) :: particle_mpitype

    ! Lists that hold the Send and Recv connections per grid level
    ! This list must be pre-compiled before the first call to 'connect'
    ! is being made. The reason for this is because it is expensive
    ! to compute every single time a connect is being made.
    !
    ! Dimensions contain:
    !   Dim 1: Information about a specific connection
    !   Dim 2: The different connections
    !
    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of receiving process
    !   Field 2: Rank of sending process
    !   Field 3: ID of receiving grid
    !   Field 4: ID of sending grid
    !   Field 5: Which face (1..26) to receive
    !   Field 6: Which face (1..26) to send
    !   Field 7: Message tag (for MPI)
    !   Field 8: Geometry exchange flag
    INTEGER(intk), ALLOCATABLE :: sendConns(:,:), recvConns(:,:)

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendReqs(:), recvReqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nSend, nRecv, nRecvFaces
    INTEGER(int32), ALLOCATABLE :: sendList(:), recvList(:)
    INTEGER(intk), ALLOCATABLE :: recvIdxList(:,:)

    ! Number of send and receive connections
    INTEGER(intk) :: iSend = 0, iRecv = 0

    ! Counters for send- and receive operations (for locations in
    ! send and receive buffers)
    INTEGER(intk) :: sendCounter, recvCounter

    ! Variable to indicate if the connection information has
    ! been created.
    LOGICAL :: isInit = .FALSE.

    INTEGER(intk), PARAMETER :: facelist(4,26) = RESHAPE((/ &
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
        3, 2, 4, 6 /), SHAPE(facelist))

    INTEGER(intk), PARAMETER :: facenbr(26) = (/ &
        2, &
        1, &
        4, &
        3, &
        6, &
        5, &
        12, &
        11, &
        14, &
        13, &
        8, &
        7, &
        10, &
        9, &
        18, &
        17, &
        16, &
        15, &
        26, &
        25, &
        24, &
        23, &
        22, &
        21, &
        20, &
        19 /)

    ! These patterns come from a Python program, however, it was discovered,
    ! that the GC fcorr-stencils need in inner corners a special rescue
    ! neighbor. fcorr-stencils exist from 3..ii-1 and so on. In the inner
    ! corners in a three-grid configuration, in the front-left (8),
    ! front-top (10) and right-top (16), these need data to come from
    ! the grids that treat the face-normal velocity also face-normal. This is
    ! because the face-normal and face-tangential velocities are
    ! treated diffrently.
    !
    ! The original un-altered order of these faces are in the comment
    ! behind the adapted ones.
    INTEGER(intk), PARAMETER :: rescue_dir(7, 26) = RESHAPE((/ &
        1, 0, 0, 0, 0, 0, 0, &
        2, 0, 0, 0, 0, 0, 0, &
        3, 0, 0, 0, 0, 0, 0, &
        4, 0, 0, 0, 0, 0, 0, &
        5, 0, 0, 0, 0, 0, 0, &
        6, 0, 0, 0, 0, 0, 0, &
        7, 1, 3, 0, 0, 0, 0, &
        8, 4, 1, 0, 0, 0, 0, &  ! 8, 1, 4, 0, 0, 0, 0, &
        9, 1, 5, 0, 0, 0, 0, &
        10, 6, 1, 0, 0, 0, 0, &  ! 10, 1, 6, 0, 0, 0, 0, &
        11, 2, 3, 0, 0, 0, 0, &
        12, 2, 4, 0, 0, 0, 0, &
        13, 2, 5, 0, 0, 0, 0, &
        14, 2, 6, 0, 0, 0, 0, &
        15, 3, 5, 0, 0, 0, 0, &
        16, 6, 3, 0, 0, 0, 0, &  ! 16, 3, 6, 0, 0, 0, 0, &
        17, 4, 5, 0, 0, 0, 0, &
        18, 4, 6, 0, 0, 0, 0, &
        19, 7, 9, 15, 1, 3, 5, &
        20, 7, 10, 16, 1, 3, 6, &
        21, 8, 9, 17, 1, 4, 5, &
        22, 8, 10, 18, 1, 4, 6, &
        23, 11, 13, 15, 2, 3, 5, &
        24, 11, 14, 16, 2, 3, 6, &
        25, 12, 13, 17, 2, 4, 5, &
        26, 12, 14, 18, 2, 4, 6 /), SHAPE(rescue_dir))

    INTEGER(intk), PARAMETER :: rescue_nbr(7, 26) = RESHAPE((/ &
        2, 0, 0, 0, 0, 0, 0, &
        1, 0, 0, 0, 0, 0, 0, &
        4, 0, 0, 0, 0, 0, 0, &
        3, 0, 0, 0, 0, 0, 0, &
        6, 0, 0, 0, 0, 0, 0, &
        5, 0, 0, 0, 0, 0, 0, &
        12, 11, 8, 0, 0, 0, 0, &
        11, 7, 12, 0, 0, 0, 0, &  ! 11, 12, 7, 0, 0, 0, 0, &
        14, 13, 10, 0, 0, 0, 0, &
        13, 9, 14, 0, 0, 0, 0, &  ! 13, 14, 9, 0, 0, 0, 0, &
        8, 7, 12, 0, 0, 0, 0, &
        7, 8, 11, 0, 0, 0, 0, &
        10, 9, 14, 0, 0, 0, 0, &
        9, 10, 13, 0, 0, 0, 0, &
        18, 17, 16, 0, 0, 0, 0, &
        17, 15, 18, 0, 0, 0, 0, &  ! 17, 18, 15, 0, 0, 0, 0, &
        16, 15, 18, 0, 0, 0, 0, &
        15, 16, 17, 0, 0, 0, 0, &
        26, 25, 24, 22, 23, 21, 20, &
        25, 26, 23, 21, 24, 22, 19, &
        24, 23, 26, 20, 25, 19, 22, &
        23, 24, 25, 19, 26, 20, 21, &
        22, 21, 20, 26, 19, 25, 24, &
        21, 22, 19, 25, 20, 26, 23, &
        20, 19, 22, 24, 21, 23, 26, &
        19, 20, 21, 23, 22, 24, 25 /), SHAPE(rescue_nbr))

        ! Publicly callable functions of module
        PUBLIC :: init_particle_connect, particle_connect, finish_particle_connect, get_target_grid

CONTAINS


    SUBROUTINE particle_connect(particle_list)

        IMPLICIT NONE

        ! subroutine argument
        TYPE(particle_list_t), INTENT(inout) :: particle_list

        !local variables
        INTEGER(intk) :: i, j, iproc, pos, num
        INTEGER(intk) :: destgrid, destproc
        INTEGER(intk) :: iprocnbr, cSend, cRecv
        INTEGER(intk) :: active_np_old  ! for safety checks

        ! we use "intk" instead of "ifk" (limits numbers)
        INTEGER(intk), ALLOCATABLE :: npsend(:)
        INTEGER(intk), ALLOCATABLE :: nprecv(:)
        INTEGER(intk), ALLOCATABLE :: ndispsend(:)
        INTEGER(intk), ALLOCATABLE :: ndisprecv(:)
        INTEGER(intk), ALLOCATABLE :: sendind(:)

        ! for periodic boundaries
        REAL(realk) :: old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, &
         new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz

        IF ( .NOT. isInit ) THEN
            WRITE(*,*) 'Particle connect not initialized'
            CALL errr(__FILE__, __LINE__)
        END IF

        active_np_old = particle_list%active_np

        DO i = 1, particle_list%ifinal

            ! jumping inactive particles
            IF ( particle_list%particles(i)%is_active /= 1 ) THEN
                CYCLE
            END IF

            ! setting the destination of particle (quo vadis, particle?)
            CALL get_target_grid(particle_list%particles(i), destgrid, destproc)


            IF ( destproc > numprocs .OR. destproc < 0 ) THEN
                WRITE(*,*) 'Obviously ill-addressed particle to proc', destproc
                CALL errr(__FILE__, __LINE__)
            END IF

            ! triage of particles
            IF ( particle_list%particles(i)%igrid == destgrid ) THEN
                ! particle stays on grid
                CALL update_particle_cell( particle_list%particles(i) )
            ELSE
                ! get old and new grid boundaries for periodic boundary handling
                ! (at this point igrid is still NOT updated, meaning particle%igrid = old igrid)
                CALL get_bbox(old_minx, old_maxx, old_miny, old_maxy, old_minz, old_maxz, particle_list%particles(i)%igrid)
                CALL get_bbox(new_minx, new_maxx, new_miny, new_maxy, new_minz, new_maxz, destgrid)

                ! coordinate manipulation of particles passing periodic boundaries (JULIUS: should this be sourced out into its own routine?)
                IF (particle_list%particles(i)%x < new_minx) THEN
                    particle_list%particles(i)%x = new_maxx - ABS(particle_list%particles(i)%x - old_minx)
                END IF

                IF (new_maxx < particle_list%particles(i)%x) THEN
                    particle_list%particles(i)%x = new_minx + ABS(particle_list%particles(i)%x - old_maxx)
                END IF

                IF (particle_list%particles(i)%y < new_miny) THEN
                    particle_list%particles(i)%y = new_maxy - ABS(particle_list%particles(i)%y - old_miny)
                END IF

                IF (new_maxy < particle_list%particles(i)%y) THEN
                    particle_list%particles(i)%y = new_miny + ABS(particle_list%particles(i)%y - old_maxy)
                END IF

                IF (particle_list%particles(i)%z < new_minz) THEN
                    particle_list%particles(i)%z = new_maxz - ABS(particle_list%particles(i)%z - old_minz)
                END IF

                IF (new_maxz < particle_list%particles(i)%z) THEN
                    particle_list%particles(i)%z = new_minz + ABS(particle_list%particles(i)%z - old_maxz)
                END IF

                ! particle changes the grid
                IF ( destproc == myid ) THEN
                    ! particle remains on process
                    particle_list%particles(i)%igrid = destgrid
                    CALL set_particle_cell( particle_list%particles(i) )
                ELSE
                    ! particle is marked for MPI transfer
                    particle_list%particles(i)%iproc = destproc
                    particle_list%particles(i)%igrid = destgrid
                END IF

            END IF
        END DO

        ! --- step 1: The marking is done (grid and proc indicate destination). Done.

        ALLOCATE( npsend(iSend) )
        npsend = 0

        DO i = 1, particle_list%ifinal
            ! jumping inactive particles
            IF ( particle_list%particles(i)%is_active /= 1 ) THEN
                CYCLE
            END IF
            ! jumping local particles
            IF ( particle_list%particles(i)%iproc == myid ) THEN
                CYCLE
            END IF
            ! search for the process to send to (only checks few "neighbor processes")
            DO iproc = 1, iSend
                IF ( sendConns(1,iproc) == particle_list%particles(i)%iproc ) THEN
                    npsend(iproc) = npsend(iproc) + 1
                END IF
            END DO
        END DO

        ! --- step 2: The counting is done. Done.

        ALLOCATE( nprecv(iRecv) )
        nprecv = -1

        ! posting non-blocking (!) receives
        ! int MPI_Irecv(void *buf, int count,
        !     MPI_Datatype datatype, int source,
        !     int tag, MPI_Comm comm, MPI_Request *request)
        DO i = 1, iRecv
            iprocnbr = recvConns(2, i)
            CALL MPI_Irecv( nprecv(i), 1, mglet_mpi_int, &
            iprocnbr, 123, MPI_COMM_WORLD, recvreqs(i) )
        END DO

        ! posting non-blocking (!) sends
        ! int MPI_Isend(const void *buf, int count,
        !     MPI_Datatype datatype, int dest, int tag,
        !     MPI_Comm comm, MPI_Request *request)
        DO i = 1, iSend
            iprocnbr = sendConns(1, i)
            CALL MPI_Isend( npsend(i), 1, mglet_mpi_int, &
            iprocnbr, 123, MPI_COMM_WORLD, sendreqs(i) )
        END DO

        ! --- step 3: The communication has been launched (not finished!). Open.

        ! displacements for start of section for one destination
        ALLOCATE( ndispsend(iSend) )
        ndispsend = -1; ndispsend(1) = 1
        DO i = 2, iSend
            IF ( npsend(i-1) < 0 ) THEN
                WRITE(*,*) 'Negative npsend value'
                CALL errr(__FILE__, __LINE__)
            END IF
            ndispsend(i) = ndispsend(i-1) + npsend(i-1)
        END DO

        ! allocate send buffer and copy particles insections
        sizeSendBuf = SUM(npsend)
        ALLOCATE(sendind(sizeSendBuf))
        ALLOCATE(sendBufParticle(sizeSendBuf))

        j = 1
        DO i = 1, particle_list%ifinal
            ! jumping inactive particles
            IF ( particle_list%particles(i)%is_active /= 1 ) THEN
                CYCLE
            END IF
            ! jumping local particles
            IF ( particle_list%particles(i)%iproc == myid ) THEN
                CYCLE
            END IF

            ! collect indices of particles list entries that will be empty after MPI send
            sendind(j) = i
            j = j + 1

            ! buffer is filled
            DO iproc = 1, iSend
                IF ( sendConns(1, iproc) == particle_list%particles(i)%iproc ) THEN
                    pos = ndispsend(iproc)
                    IF ( pos > sizeSendBuf ) THEN
                        WRITE(*,*) 'Send buffer size exceeded'
                        CALL errr(__FILE__, __LINE__)
                    ELSE IF ( pos < 1 ) THEN
                        WRITE(*,*) 'Invalid buffer index', pos
                        CALL errr(__FILE__, __LINE__)
                    ELSE
                        ! copy particle into buffer
                        sendBufParticle(pos) = particle_list%particles(i)
                    END IF
                    ! increment the position wherer future particle for
                    ! this destination process will be stored in the buffer
                    ndispsend(iproc) = ndispsend(iproc) + 1

                    ! setting the local particle as inactive (active in buffer)
                    particle_list%particles(i)%is_active = 0
                    particle_list%active_np = particle_list%active_np - 1
                END IF
            END DO
        END DO

        ! resetting after incrementation
        ndispsend = -1; ndispsend(1) = 1
        DO i = 2, iSend
            ndispsend(i) = ndispsend(i-1) + npsend(i-1)
        END DO

        ! buffer must be full with valid particles without gaps
        DO i = 1, sizeSendBuf
            IF ( sendBufParticle(pos)%is_active < 1 ) THEN
                WRITE(*,*) 'Invalid send buffer entry at ', i
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        ! WRITE(*,*) sendBufParticle

        ! --- step 4: Displacements determined and send buffer filled. Done.

        ! checking if communication done (one call should suffice...)
        CALL MPI_Waitall(iSend, sendreqs, MPI_STATUSES_IGNORE)
        CALL MPI_Waitall(iRecv, recvreqs, MPI_STATUSES_IGNORE)

        ! WRITE(*,*) myid, 'npsend:', npsend, 'to', sendConns(1, 1:iSend), 'nprecv:', nprecv, 'from', recvConns(2, 1:iRecv)

        ! --- step 5: Finishing the communication of particle numbers. Done.


        ! displacements for start of section for one source
        ALLOCATE( ndisprecv(iRecv) )
        ndisprecv = -1; ndisprecv(1) = 1
        DO i = 2, iRecv
            IF ( nprecv(i-1) < 0 ) THEN
                WRITE(*,*) 'Invalid number of received particles'
                CALL errr(__FILE__, __LINE__)
            END IF
            ndisprecv(i) = ndisprecv(i-1) + nprecv(i-1)
        END DO

        sizeRecvBuf = SUM(nprecv)
        ALLOCATE( recvBufParticle(sizeRecvBuf) )

        ! TO DO: Ab jetzt kann gestestet werden, ob all ankommenden Partikel in die Liste passen!
        ! Entsprechend kann die List erweitert oder sogar gekürzt werden, während MPI für Partciel läuft.

        ! Check if list is long enough and add additional space if not
        IF ( particle_list%max_np - particle_list%active_np < sizeRecvBuf) THEN

            CALL enlarge_particle_list(particle_list, INT(1.0 * (sizeRecvBuf - (particle_list%max_np - particle_list%active_np))))

        END IF

        ! --- step 6: Allocating the recieve buffer and displacements. Done.


        ! posting non-blocking (!) receives
        ! int MPI_Irecv(void *buf, int count,
        !     MPI_Datatype datatype, int source,
        !     int tag, MPI_Comm comm, MPI_Request *request)
        cRecv = 0
        ! WRITE(*,*) myid, 'cons:', recvConns(2, 1:iRecv), ndisprecv, nprecv
        DO i = 1, iRecv
            iprocnbr = recvConns(2, i)
            pos = ndisprecv(i)
            num = nprecv(i)
            IF ( num > 0 ) THEN
                cRecv = cRecv + 1
                CALL MPI_Irecv( recvBufParticle(pos), num, particle_mpitype, &
                iprocnbr, 321, MPI_COMM_WORLD, recvreqs(cRecv) )
            END IF
        END DO

        ! posting non-blocking (!) sends
        ! int MPI_Isend(const void *buf, int count,
        !     MPI_Datatype datatype, int dest, int tag,
        !     MPI_Comm comm, MPI_Request *request)
        cSend = 0
        DO i = 1, iSend
            iprocnbr = sendConns(1, i)
            pos = ndispsend(i)
            num = npsend(i)
            IF ( num > 0 ) THEN
                cSend = cSend + 1
                CALL MPI_Isend( sendBufParticle(pos), num, particle_mpitype, &
                iprocnbr, 321, MPI_COMM_WORLD, sendreqs(cSend) )
            END IF
        END DO

        ! --- step 7: The communication has been launched (not finished!). Open.


        ! checking if communication done (one call should suffice...)
        CALL MPI_Waitall(cSend, sendreqs, MPI_STATUSES_IGNORE)
        CALL MPI_Waitall(cRecv, recvreqs, MPI_STATUSES_IGNORE)

        ! WRITE(*,*) recvBufParticle

        ! assigning the new cell indices
        IF ( sizeRecvBuf > 0 ) THEN
            DO i = 1, sizeRecvBuf
                ! check if correctly delivered
                IF ( recvBufParticle(i)%is_active /= 1 ) THEN
                    WRITE(*,*) "Inactive particle delivered"
                    CALL errr(__FILE__, __LINE__)
                END IF
                ! check if correctly delivered
                IF ( recvBufParticle(i)%iproc /= myid ) THEN
                    WRITE(*,*) "Particle delivered to wrong proc", i, recvBufParticle(i)%iproc, myid
                    CALL errr(__FILE__, __LINE__)
                END IF
                CALL set_particle_cell( recvBufParticle(i) )
            END DO
        END IF

        ! --- step 8: Finishing the communication of actual particles. Done.

        ! Copy recieved particles into the list
        ! CAUTION: at this point, particle_list%particles(particle_list%ifinal)%is_active might be /= 1 ("empty")

        CALL integrate_particles(particle_list, sendind)

        ! Some safety checks

        ! BARRIER ONLY FOR DEGUGGING -- TEMPORARY <----------------------------------------------- TODO : remove
        CALL MPI_Barrier(MPI_COMM_WORLD)

        CALL print_list_status(particle_list)

        IF (particle_list%active_np < active_np_old + sizeRecvBuf - sizeSendBuf) THEN
            SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, '("WARNING on proc ", I0, ": Particle list holds FEWER particles than expected!")') myid
                    CASE ("verbose")
                        WRITE(*, '("WARNING on proc ", I0, ": Particle list holds FEWER particles than expected!")') myid
            END SELECT
        END IF

        IF (particle_list%active_np > active_np_old + sizeRecvBuf - sizeSendBuf) THEN
            SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, '("WARNING on proc ", I0, ": Particle list holds MORE particles than expected!")') myid
                    CASE ("verbose")
                        WRITE(*, '("WARNING on proc ", I0, ": Particle list holds MORE particles than expected!")') myid
            END SELECT
        END IF

        IF (particle_list%ifinal /= particle_list%active_np) THEN
            SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, '("WARNING on proc ", I0, ": my_particle_list%active_np (", I0, ") does not coincide with my_particle_list%ifinal (", I0, ")" )') &
                        myid, particle_list%active_np, particle_list%ifinal
                    CASE ("verbose")
                        WRITE(*, '("WARNING on proc ", I0, ": my_particle_list%active_np (", I0, ") does not coincide with my_particle_list%ifinal (", I0, ")" )') &
                        myid, particle_list%active_np, particle_list%ifinal
            END SELECT
        END IF

        ! --- step 9: Received particles have been copied into list. Done.

        DEALLOCATE( npsend )
        DEALLOCATE( nprecv )
        DEALLOCATE( ndispsend )
        DEALLOCATE( sendBufParticle )
        DEALLOCATE( ndisprecv )
        DEALLOCATE( recvBufParticle )

        ! --- step 10: Clearing the buffers. Done.

    END SUBROUTINE particle_connect


    SUBROUTINE init_particle_connect()
        INTEGER(intk) :: i, iface, igrid
        INTEGER(intk) :: iface1, iface2, iface3
        INTEGER(intk) :: itypbc1, itypbc2, itypbc3
        INTEGER(intk) :: iprocnbr, itypbc, inbrface, inbrgrid

        INTEGER(int32), ALLOCATABLE :: maxTag(:)
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)

        INTEGER(intk) :: neighbours(26)
        INTEGER :: iexchange

        ! Maximum number of connections for "simple" cases is number
        ! of grids*26. However, due to the possible prescence of
        ! precursors etc, we add a few more.
        maxConns = INT((nMyGrids+1)*26.0*1.2, intk)
        ALLOCATE(sendConns(8, maxConns))
        ALLOCATE(recvConns(8, maxConns))
        sendConns = 0
        recvConns = 0

        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(recvIdxList(3, maxConns))
        ALLOCATE(sendList(numprocs))
        ALLOCATE(recvList(numprocs))
        ALLOCATE(sendReqs(numprocs))
        ALLOCATE(recvReqs(numprocs))
        recvIdxList = 0
        sendList = 0
        recvList = 0

        ALLOCATE(maxTag(0:numprocs-1))
        ALLOCATE(sendcounts(0:numprocs-1))
        ALLOCATE(sdispls(0:numprocs-1))
        ALLOCATE(recvcounts(0:numprocs-1))
        ALLOCATE(rdispls(0:numprocs-1))
        maxTag = 0
        sendcounts = 0
        sdispls = 0
        recvcounts = 0
        rdispls = 0
        nRecv = 0

        ! SIMON: Hier sammeln wir zunächst zu viel Information, die für die
        ! Partikel nicht benötigt wird. Style to be improved...

        DO i = 1, nMyGrids

            ! getting the grid parameters
            igrid = myGrids(i)
            CALL get_neighbours(neighbours, igrid)

            ! Check surfaces of grid
            DO iface = 1, 6
                itypbc = itypboconds(1, iface, igrid)

                IF (itypbc == 7 .OR. itypbc == 19) THEN

                    CALL get_nbrs(iface, neighbours, inbrgrid, inbrface)
                    IF (inbrgrid == 0) THEN
                        CYCLE
                    END IF
                    iprocnbr = idprocofgrd(inbrgrid)
                    IF (iprocnbr == myid) THEN
                        CYCLE
                    END IF

                    ! only if neighbor not already listed
                    IF ( sendcounts(iprocnbr) == 0 ) THEN
                        iexchange = 1
                        nRecv = nRecv + 1
                        maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                        recvConns(1, nRecv) = myid      ! Receiving process (this process)
                        recvConns(2, nRecv) = iprocnbr  ! Sending process (neighbour process)
                        recvConns(3, nRecv) = igrid     ! Receiving grid (on current process)
                        recvConns(4, nRecv) = inbrgrid  ! Sending grid (on neighbour process)
                        recvConns(5, nRecv) = iface     ! Which face receive (1..26)
                        recvConns(6, nRecv) = inbrface  ! Which face receive from (sending face) (1..26)
                        recvConns(7, nRecv) = maxTag(iprocnbr)  ! Message tag
                        recvConns(8, nRecv) = iexchange  ! Geometry exchange flag

                        sendcounts(iprocnbr) = SIZE(recvConns, 1)   ! not an increment
                    END IF
                END IF
            END DO

            ! Check lines fo grid
            DO iface = 7, 18
                iface1 = facelist(2, iface)
                iface2 = facelist(3, iface)

                itypbc1 = itypboconds(1, iface1, igrid)
                itypbc2 = itypboconds(1, iface2, igrid)

                IF (itypbc1 == 7 .OR. itypbc2 == 7 .OR. &
                    itypbc1 == 19 .OR. itypbc2 == 19) THEN

                    CALL get_nbrs(iface, neighbours, inbrgrid, inbrface)
                    IF (inbrgrid == 0) THEN
                        CYCLE
                    END IF
                    iprocnbr = idprocofgrd(inbrgrid)
                    IF (iprocnbr == myid) THEN
                        CYCLE
                    END IF

                    ! only if neighbor not already listed
                    IF ( sendcounts(iprocnbr) == 0 ) THEN
                        iexchange = 1
                        nRecv = nRecv + 1
                        maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                        recvConns(1, nRecv) = myid      ! Receiving process (this process)
                        recvConns(2, nRecv) = iprocnbr  ! Sending process (neighbour process)
                        recvConns(3, nRecv) = igrid     ! Receiving grid (on current process)
                        recvConns(4, nRecv) = inbrgrid  ! Sending grid (on neighbour process)
                        recvConns(5, nRecv) = iface     ! Which face receive (1..26)
                        recvConns(6, nRecv) = inbrface  ! Which face receive from (sending face) (1..26)
                        recvConns(7, nRecv) = maxTag(iprocnbr)  ! Message tag
                        recvConns(8, nRecv) = iexchange  ! Geometry exchange flag

                        sendcounts(iprocnbr) = SIZE(recvConns, 1)   ! not an increment
                    END IF
                END IF
            END DO

            ! Check corners of grid
            DO iface = 19, 26
                iface1 = facelist(2, iface)
                iface2 = facelist(3, iface)
                iface3 = facelist(4, iface)
                itypbc1 = itypboconds(1, iface1, igrid)
                itypbc2 = itypboconds(1, iface2, igrid)
                itypbc3 = itypboconds(1, iface3, igrid)

                IF (itypbc1 == 7 .OR. itypbc2 == 7 .OR. itypbc3 == 7 .OR. &
                    itypbc1 == 19 .OR. itypbc2 == 19 .OR. itypbc3 == 19) THEN

                    CALL get_nbrs(iface, neighbours, inbrgrid, inbrface)
                    IF (inbrgrid == 0) THEN
                        CYCLE
                    END IF
                    iprocnbr = idprocofgrd(inbrgrid)
                    IF (iprocnbr == myid) THEN
                        CYCLE
                    END IF

                    ! only if neighbor not already listed
                    IF ( sendcounts(iprocnbr) == 0 ) THEN
                        iexchange = 1
                        nRecv = nRecv + 1
                        maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                        recvConns(1, nRecv) = myid      ! Receiving process (this process)
                        recvConns(2, nRecv) = iprocnbr  ! Sending process (neighbour process)
                        recvConns(3, nRecv) = igrid     ! Receiving grid (on current process)
                        recvConns(4, nRecv) = inbrgrid  ! Sending grid (on neighbour process)
                        recvConns(5, nRecv) = iface     ! Which face receive (1..26)
                        recvConns(6, nRecv) = inbrface  ! Which face receive from (sending face) (1..26)
                        recvConns(7, nRecv) = maxTag(iprocnbr)  ! Message tag
                        recvConns(8, nRecv) = iexchange  ! Geometry exchange flag

                        sendcounts(iprocnbr) = SIZE(recvConns, 1)   ! not an increment
                    END IF
                END IF
            END DO
        END DO

        ! Sort recvConns by process ID
        CALL sort_conns_unique( recvConns(:,1:nRecv) )
        iRecv = nRecv

        ! Calculate sdispl offset (send)
        DO i=1,numprocs-1
            ! = value is either 0 (not a neighbor) or 8 (a neighbor)
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Calculate rdispl offset (receive)
        DO i=1,numprocs-1
            ! = value is either 0 (not a neighbor) or 8 (a neighbor)
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        ! Check that number of connections fit in array
        iSend = (rdispls(numprocs-1) + recvcounts(numprocs-1)) / SIZE(sendConns, 1)
        IF (iSend > maxConns) THEN
            WRITE(*,*) "Number of connections exceeded on process ", myid
            WRITE(*,*) "maxConns =", maxConns, "nMyGrids =", nMyGrids, &
                "iSend = ", iSend
            CALL errr(__FILE__, __LINE__)
        END IF

        ! int MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
        ! const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
        ! const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
        ! MPI_Comm comm)

        ! Exchange connection information
        CALL MPI_Alltoallv( &
            recvConns(1, 1), sendcounts, sdispls, MPI_INTEGER, &
            sendConns(1, 1), recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        !IF ( myid == 0 ) THEN
            WRITE(*,*) 'I am proc:', myid
            WRITE(*,*) 'I own grids: '
            DO i = 1, nmygrids
                WRITE(*,*) '    - grid ', mygrids(i)
            END DO
            WRITE(*,*) ' - I receive from the following ', iRecv, 'processes:'
            DO i = 1, iRecv
                WRITE(*,*) '    - proc ', recvConns(2, i)
            END DO
            WRITE(*,*) ' - I send to the following ', iSend, 'processes:'
            DO i = 1, iSend
                WRITE(*,*) '    - proc ', sendConns(1, i)
            END DO

        !END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

        nRecv = 0

        ! creating the MPI data type
        CALL create_particle_mpitype( particle_mpitype )
        isInit = .TRUE.

    END SUBROUTINE init_particle_connect


    SUBROUTINE get_nbrs(iface, neighbours, nbrgrid, nbrface)
        INTEGER(intk), INTENT(IN) :: iface
        INTEGER(intk), INTENT(IN) :: neighbours(26)
        INTEGER(intk), INTENT(OUT) :: nbrgrid
        INTEGER(intk), INTENT(OUT) :: nbrface

        INTEGER(intk) :: n_rescue, i, dir
        INTEGER(intk) :: iface1, iface2, iface3
        INTEGER(intk) :: itypbc1, itypbc2, itypbc3

        ! Should be 7...
        n_rescue = SIZE(rescue_nbr, 1)

        ! 0 means no connect
        nbrgrid = 0
        nbrface = 0
        DO i = 1, n_rescue
            dir = rescue_dir(i, iface)

            ! rescue_dir is ordered and when a 0 is encountered there is
            ! nothing more to do...
            IF (dir == 0) THEN
                EXIT
            END IF

            ! If there is a neighbour in this position, use this
            IF (neighbours(dir) > 0) THEN
                nbrgrid = neighbours(dir)
                nbrface = rescue_nbr(i, iface)

                ! Check if this is suited for a connect (symmetry req.)
                ! This require knowledge of the global grid structure -
                ! currently this is OK.
                IF (nbrface > 18) THEN
                    ! Get adjacent primary faces
                    iface1 = facelist(2, nbrface)
                    iface2 = facelist(3, nbrface)
                    iface3 = facelist(4, nbrface)

                    ! Get type of BC on these
                    itypbc1 = itypboconds(1, iface1, nbrgrid)
                    itypbc2 = itypboconds(1, iface2, nbrgrid)
                    itypbc3 = itypboconds(1, iface3, nbrgrid)

                    ! If none of the neighboring faces are CON or CO1, the connect
                    ! should not be carried out - check next neighbour
                    IF ((.NOT. (itypbc1 == 7 .OR. itypbc1 == 19)) .AND. &
                            (.NOT. (itypbc2 == 7 .OR. itypbc2 == 19)) .AND. &
                            (.NOT. (itypbc3 == 7 .OR. itypbc3 == 19))) THEN
                        ! Reset neighbour information and cycle loop
                        nbrgrid = 0
                        nbrface = 0
                        CYCLE
                    END IF
                END IF

                ! If sofar, all is good!
                EXIT
            END IF
        END DO
    END SUBROUTINE get_nbrs



    SUBROUTINE finish_particle_connect()
        isInit = .FALSE.

        DEALLOCATE(sendConns)
        DEALLOCATE(recvConns)

        DEALLOCATE(recvIdxList)
        DEALLOCATE(sendList)
        DEALLOCATE(recvList)
        DEALLOCATE(sendReqs)
        DEALLOCATE(recvReqs)

    END SUBROUTINE finish_particle_connect



    SUBROUTINE sort_conns_unique(list)
        ! Input array to be sorted
        INTEGER(int32), INTENT(inout) :: list(:,:)

        INTEGER(intk) :: i,j

        ! Temporary storage
        INTEGER(int32) :: temp(8)

        IF (SIZE(list, 1) /= SIZE(temp)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Sort by sending processor number (field 2)
        DO i = 2, SIZE(list, 2)
            j = i - 1
            temp(:) = list(:,i)
            DO WHILE (j >= 1)
                IF (list(2,j) > temp(2)) THEN
                    list(:,j+1) = list(:,j)
                    j = j - 1
                ELSE
                    EXIT
                END IF
            END DO
            list(:,j+1) = temp(:)
        END DO

        ! Check for redundant entries
        DO i = 2, SIZE(list, 2)
            IF ( list(2,i) == list(2,i-1) ) THEN
                WRITE(*,*) 'Redundant listing of neighbor process ', list(2,i)
                CALL errr(__FILE__, __LINE__)
            END IF
            IF ( list(2,i) == myid ) THEN
                WRITE(*,*) 'Self connection listed at ', list(2,i)
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

    END SUBROUTINE sort_conns_unique



    SUBROUTINE get_target_grid( particle, destgrid, destproc )

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(in) :: particle
        INTEGER(intk), INTENT(out) :: destgrid
        INTEGER(intk), INTENT(out) :: destproc

        ! local variables
        INTEGER(intk) :: iface, neighbours(26)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        ! getting the box of last grid the particla
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, particle%igrid)

        ! intialization (will be overwritten is particle left grid)
        iface = -1

        ! checking the geometrical relation
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

        IF ( iface == 0 ) THEN
            ! particle stays on the same grid
            destgrid = particle%igrid
            destproc = particle%iproc

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, '("Destination Proc: ", I0)') destproc
                    WRITE(*, '("Destination grid: ", I0)') destgrid
                    IF (myid == 0) THEN
                        WRITE(*, *) " "
                    END IF
            END SELECT

            IF ( destproc /= myid ) THEN
                WRITE(*,*) 'Inconsistent particle parameters'
                CALL errr(__FILE__, __LINE__)
            END IF

        ELSE IF ( iface > 0 ) THEN
            ! particle moves across grid boundary
            CALL get_neighbours(neighbours, particle%igrid)
            destgrid = neighbours(iface)
            destproc = idprocofgrd(destgrid)

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, '("Destination Proc: ", I0)') destproc
                    WRITE(*, '("Destination grid: ", I0)') destgrid
                    IF (myid == 0) THEN
                        WRITE(*, *) " "
                    END IF
            END SELECT

            IF ( destproc == myid ) THEN
                destproc = particle%iproc
            END IF

        ELSE
            WRITE(*,*) 'Undefined behaviour'
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE get_target_grid



    SUBROUTINE create_particle_mpitype(dtype)
        ! Subrouitine arguments
        TYPE(MPI_Datatype), INTENT(inout) :: dtype

        ! Local variables
        INTEGER(intk) :: i
        TYPE(baseparticle_t) :: foo
        INTEGER(MPI_ADDRESS_KIND) :: base, disp(particle_mpi_elems)
        INTEGER(int32) :: blocklen(particle_mpi_elems)
        TYPE(MPI_Datatype) :: types(particle_mpi_elems)
        TYPE(MPI_Datatype) :: triple_int_mpi_type

        CALL MPI_Type_contiguous(3, mglet_mpi_int, triple_int_mpi_type)

        CALL MPI_Get_address(foo%is_active, disp(1))
        CALL MPI_Get_address(foo%ipart, disp(2))
        CALL MPI_Get_address(foo%iproc, disp(3))
        CALL MPI_Get_address(foo%igrid, disp(4))
        CALL MPI_Get_address(foo%ijkcell, disp(5))
        CALL MPI_Get_address(foo%facepath, disp(6))
        CALL MPI_Get_address(foo%x, disp(7))
        CALL MPI_Get_address(foo%y, disp(8))
        CALL MPI_Get_address(foo%z, disp(9))

        types(1) = mglet_mpi_int    ! is_active
        types(2) = mglet_mpi_int    ! ipart
        types(3) = mglet_mpi_int    ! iproc
        types(4) = mglet_mpi_int    ! igrid
        types(5) = triple_int_mpi_type  ! ijkcell(3)
        types(6) = triple_int_mpi_type  ! facepath(3)
        types(7) = mglet_mpi_real    ! x
        types(8) = mglet_mpi_real    ! y
        types(9) = mglet_mpi_real    ! z

        ! computing the displacements in byte
        base = disp(1)
        DO i = 1, particle_mpi_elems
            disp(i) = disp(i) - base
        END DO

        ! creating and submitting type
        blocklen = 1
        CALL MPI_Type_create_struct(particle_mpi_elems, &
            blocklen, disp, types, dtype)
        CALL MPI_Type_commit(dtype)

        ! cleaning up the auxiliary type
        CALL MPI_Type_free(triple_int_mpi_type)
    END SUBROUTINE create_particle_mpitype


    SUBROUTINE is_true_neigbour()



    END SUBROUTINE is_true_neigbour

    ! copy particles from recieve Buffer into passed particle list
    SUBROUTINE integrate_particles(particle_list, sendind)

        ! subroutine argument
        TYPE(particle_list_t), INTENT(inout) :: particle_list
        INTEGER(intk), INTENT(in) :: sendind(sizeSendBuf)

        !local variables
        INTEGER(intk) :: i, j

        IF (sizeSendBuf == 0 .AND. sizeRecvBuf == 0) THEN

            RETURN

        END IF

        particle_list%active_np = particle_list%active_np + sizeRecvBuf

        IF (sizeSendBuf <= sizeRecvBuf) THEN

            DO i = 1, sizeSendBuf
                particle_list%particles(sendind(i)) = recvBufParticle(i)
            END DO

            DO i = i, sizeRecvBuf ! i = sizeSendBuf
                particle_list%ifinal = particle_list%ifinal + 1
                particle_list%particles(particle_list%ifinal) = recvBufParticle(i)
            END DO

        ELSE

            DO i = 1, sizeRecvBuf
                particle_list%particles(sendind(i)) = recvBufParticle(i)
            END DO

            ! i = MAX(1, sizeRecvBuf + 1)
            DO i = i, sizeSendBuf

                IF (particle_list%ifinal == sendind(i)) THEN
                    particle_list%ifinal = particle_list%ifinal - 1
                    EXIT
                END IF

                IF (particle_list%ifinal < sendind(i)) THEN
                    EXIT
                END IF

                DO j = 1, particle_list%ifinal - sendind(i)
                    IF (particle_list%particles(particle_list%ifinal)%is_active == 1) THEN

                        particle_list%particles(sendind(i)) = particle_list%particles(particle_list%ifinal)
                        particle_list%particles(particle_list%ifinal)%is_active = 0
                        particle_list%ifinal = particle_list%ifinal - 1
                        EXIT

                    ELSE

                        particle_list%ifinal = particle_list%ifinal - 1

                    END IF
                END DO

                IF (particle_list%particles(sendind(i))%is_active /= 1) THEN
                    particle_list%ifinal = particle_list%ifinal - 1
                END IF

            END DO

        END IF

    END SUBROUTINE integrate_particles

    ! alternative way of copying recieved particles into passed particle list
    ! (dont know if this works properly due to the do concurrent ...)
    SUBROUTINE integrate_particles_concurrent(particle_list, sendind)

        ! subroutine argument
        TYPE(particle_list_t), INTENT(inout) :: particle_list
        INTEGER(intk), INTENT(in) :: sendind(sizeSendBuf)

        !local variables
        INTEGER(intk) :: i, j

        IF (sizeSendBuf == 0 .AND. sizeRecvBuf == 0) THEN

            RETURN

        END IF

        particle_list%active_np = particle_list%active_np + sizeRecvBuf

        IF (sizeSendBuf <= sizeRecvBuf) THEN

            IF (1 <= MIN(sizeSendBuf, sizeRecvBuf - sizeSendBuf)) THEN
                DO CONCURRENT (i = 1:MIN(sizeSendBuf, sizeRecvBuf - sizeSendBuf))
                    particle_list%particles(sendind(i)) = recvBufParticle(i)
                    particle_list%particles(particle_list%ifinal + i) = recvBufParticle(sizeRecvBuf + 1 - i)
                END DO
                i = MIN(sizeSendBuf, sizeRecvBuf - sizeSendBuf) + 1
            ELSE
                i = 1
            END IF

            IF (i <= sizeSendBuf) THEN
                DO CONCURRENT (j = i:sizeSendBuf)
                    particle_list%particles(sendind(j)) = recvBufParticle(j)
                END DO
            END IF

            IF (i <= sizeRecvBuf - sizeSendBuf) THEN
                DO CONCURRENT (j = i:sizeRecvBuf - sizeSendBuf)
                    particle_list%particles(particle_list%ifinal + j) = recvBufParticle(sizeRecvBuf + 1 - j)
                END DO
            END IF

            particle_list%ifinal = particle_list%ifinal + sizeRecvBuf - sizeSendBuf

        ELSE

            IF (1 <= sizeRecvBuf) THEN
                DO CONCURRENT (i = 1:sizeRecvBuf)
                    particle_list%particles(sendind(i)) = recvBufParticle(i)
                END DO
            END IF

            DO i = i, sizeSendBuf

                IF (particle_list%ifinal == sendind(i)) THEN
                    particle_list%ifinal = particle_list%ifinal - 1
                    EXIT
                END IF

                IF (particle_list%ifinal < sendind(i)) THEN
                    EXIT
                END IF

                DO j = 1, particle_list%ifinal - sendind(i)
                    IF (particle_list%particles(particle_list%ifinal)%is_active == 1) THEN

                        particle_list%particles(sendind(i)) = particle_list%particles(particle_list%ifinal)
                        particle_list%particles(particle_list%ifinal)%is_active = 0
                        particle_list%ifinal = particle_list%ifinal - 1
                        EXIT

                    ELSE

                        particle_list%ifinal = particle_list%ifinal - 1

                    END IF
                END DO

                IF (particle_list%particles(sendind(i))%is_active /= 1) THEN
                    particle_list%ifinal = particle_list%ifinal - 1
                    EXIT
                END IF

            END DO

        END IF

    END SUBROUTINE integrate_particles_concurrent

END MODULE particle_connect_mod