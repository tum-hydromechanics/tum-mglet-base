MODULE particle_exchange_mod

    USE, INTRINSIC :: ISO_C_BINDING

    USE particle_list_mod
    USE particle_statistics_mod

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
    !   Field 3: Message tag (for MPI)
    !   Field 4: Geometry exchange flag
    INTEGER(intk), ALLOCATABLE :: sendConns(:,:), recvConns(:,:)

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendReqs(:), recvReqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nRecv
    INTEGER(int32), ALLOCATABLE :: sendList(:), recvList(:)
    INTEGER(intk), ALLOCATABLE :: recvIdxList(:,:)

    ! Number of send and receive connections
    INTEGER(intk) :: iSend = 0, iRecv = 0

    ! we use "intk" instead of "ifk" (limits numbers)
    INTEGER(intk), ALLOCATABLE :: nprecv(:)
    INTEGER(intk), ALLOCATABLE :: npsend(:)
    INTEGER(intk), ALLOCATABLE :: ndispsend(:)
    INTEGER(intk), ALLOCATABLE :: ndisprecv(:)

    ! Variable to indicate if the connection information has
    ! been created.
    LOGICAL :: isInit = .FALSE.

    ! The first column just contains the number 2D faces the face (id = line)
    ! is uniquily defined by. Column 2-4 hold those 2D faces.
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

        ! Publicly callable functions of module
        PUBLIC :: init_particle_exchange, exchange_particles, finish_particle_exchange, get_target_grid !, prepare_particle_exchange

CONTAINS

!    SUBROUTINE prepare_particle_exchange(particle)
!
!        ! subroutine arguments
!        TYPE(baseparticle_t), INTENT(in) :: particle
!
!        ! local variables
!        INTEGER(intk) :: iproc
!
!        DO iproc = 1, iSend
!            IF (sendConns(1, iproc) == particle%iproc) THEN
!                npsend(iproc) = npsend(iproc) + 1
!            END IF
!        END DO
!
!    END SUBROUTINE prepare_particle_exchange

    SUBROUTINE exchange_particles(particle_list, itstep)

        IMPLICIT NONE

        ! subroutine argument
        TYPE(particle_list_t), INTENT(inout) :: particle_list
        INTEGER(intk), INTENT(in) :: itstep

        !local variables
        INTEGER(intk) :: i, j, iproc, pos, num
        INTEGER(intk) :: destgrid, destproc, iface
        INTEGER(intk) :: iprocnbr, cSend, cRecv
        INTEGER(intk) :: active_np_old  ! for safety checks

        ! we use "intk" instead of "ifk" (limits numbers)
        ! INTEGER(intk), ALLOCATABLE :: npsend(:)
        INTEGER(intk), ALLOCATABLE :: sendind(:)

        CALL start_timer(900)
        CALL start_timer(930)

        IF (.NOT. isInit) THEN
            WRITE(*,*) 'Particle connect not initialized'
            CALL errr(__FILE__, __LINE__)
        END IF

        active_np_old = particle_list%active_np

        npsend = 0
        nprecv = -1

        DO i = 1, particle_list%ifinal

            ! jumping inactive particles
            IF (particle_list%particles(i)%state < 1) THEN
                CYCLE
            END IF

            ! for particle slice statistics (must be called before update_coordinates !!!)
            CALL stop_timer(930)
            CALL associate_new_slice(particle_list%particles(i), itstep)
            CALL start_timer(930)

            ! setting the destination of particle (quo vadis, particle?)
            CALL get_target_grid(particle_list%particles(i), destgrid, destproc, iface)

            IF (destproc > numprocs .OR. destproc < 0) THEN
                WRITE(*,*) 'Obviously ill-addressed particle to proc', destproc
                CALL errr(__FILE__, __LINE__)
            END IF

            ! coordinate manipulation of particles passing periodic boundaries
            ! (at this point igrid is still NOT updated, meaning particle%igrid is still the "old" grid)
            CALL update_coordinates(particle_list%particles(i), destgrid, iface)

            ! triage of particles
            IF (particle_list%particles(i)%igrid == destgrid) THEN

                ! particle stays on grid
                CALL update_particle_cell(particle_list%particles(i))

            ELSE

                ! for particle statistics
                CALL stop_timer(930)
                CALL deregister_particle(particle_list%particles(i), itstep)
                CALL start_timer(930)

                ! particle changes the grid
                IF (destproc == myid) THEN

                    ! particle remains on process
                    particle_list%particles(i)%igrid = destgrid
                    CALL set_particle_cell(particle_list%particles(i))

                    ! for particle statistics
                    CALL stop_timer(930)
                    CALL register_particle(particle_list%particles(i), itstep)
                    CALL start_timer(930)

                ELSE

                    ! particle is marked for MPI transfer
                    particle_list%particles(i)%iproc = destproc
                    particle_list%particles(i)%igrid = destgrid

                    ! search for the process to send to (only checks few "neighbor processes")
                    DO iproc = 1, iSend
                        IF (sendConns(1, iproc) == particle_list%particles(i)%iproc) THEN
                            npsend(iproc) = npsend(iproc) + 1
                        END IF
                    END DO

                END IF

            END IF
        END DO

        ! --- step 1: The marking is done (grid and proc indicate destination). Done.
        ! --- step 2: The counting is done. Done.

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

        ! JULIUS: would be nice to not iterate over the whole particle list twice. Maybe there is a way to allocate sendind before the first iteration?
        ! Maybe use sending from previous exchnage ?
        j = 1
        DO i = 1, particle_list%ifinal
            ! jumping inactive particles
            IF (particle_list%particles(i)%state < 1) THEN
                CYCLE
            END IF
            ! jumping local particles
            IF (particle_list%particles(i)%iproc == myid) THEN
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
                    particle_list%particles(i)%state = -1
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
            IF ( sendBufParticle(i)%state < 1 ) THEN
                WRITE(*,*) 'Proc', myid, ': Invalid send buffer entry at ', i
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
        ndisprecv = -1; ndisprecv(1) = 1
        DO i = 2, iRecv
            IF ( nprecv(i-1) < 0 ) THEN
                WRITE(*,*) 'Invalid number of received particles'
                CALL errr(__FILE__, __LINE__)
            END IF
            ndisprecv(i) = ndisprecv(i-1) + nprecv(i-1)
        END DO

        sizeRecvBuf = SUM(nprecv)
        ALLOCATE(recvBufParticle(sizeRecvBuf))

        ! TO DO: Ab jetzt kann gestestet werden, ob all ankommenden Partikel in die Liste passen!
        ! Entsprechend kann die List erweitert oder sogar gekürzt werden, während MPI für Partciel läuft.

        ! Check if list is long enough and add additional space if not
        IF ( particle_list%max_np - particle_list%active_np < sizeRecvBuf) THEN

            CALL reallocate_particle_list(particle_list, INT(1.0 * (sizeRecvBuf - (particle_list%max_np - particle_list%active_np))))

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
                IF ( recvBufParticle(i)%state < 1 ) THEN
                    WRITE(*,*) "Inactive particle delivered"
                    CALL errr(__FILE__, __LINE__)
                END IF
                ! check if correctly delivered
                IF ( recvBufParticle(i)%iproc /= myid ) THEN
                    WRITE(*,*) "Particle delivered to wrong proc", i, recvBufParticle(i)%iproc, myid
                    CALL errr(__FILE__, __LINE__)
                END IF

                CALL set_particle_cell(recvBufParticle(i))

                ! for gridstat
                CALL stop_timer(930)
                CALL register_particle(recvBufParticle(i), itstep)
                CALL start_timer(930)

            END DO
        END IF

        ! --- step 8: Finishing the communication of actual particles. Done.

        ! Copy recieved particles into the list
        ! CAUTION: at this point, particle_list%particles(particle_list%ifinal)%state might be < 1 ("empty")

        CALL integrate_particles(particle_list, sendind)

        ! Some safety checks

        ! BARRIER ONLY FOR DEGUGGING -- TEMPORARY <----------------------------------------------- TODO : remove
        CALL MPI_Barrier(MPI_COMM_WORLD)

        SELECT CASE (TRIM(particle_terminal))
            CASE ("none")
                CONTINUE
            CASE ("normal")
                CALL print_list_status(particle_list)
                WRITE(*, '()')
            CASE ("verbose")
                CALL print_list_status(particle_list)
                WRITE(*, '()')
        END SELECT

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

        DEALLOCATE(sendBufParticle)
        DEALLOCATE(recvBufParticle)
        DEALLOCATE(sendind)

        ! --- step 10: Clearing the buffers. Done.

        ! BARRIER ONLY FOR DEGUGGING -- TEMPORARY <----------------------------------------------- TODO : remove
        CALL MPI_Barrier(MPI_COMM_WORLD)

        CALL stop_timer(930)
        CALL stop_timer(900)

    END SUBROUTINE exchange_particles


    SUBROUTINE init_particle_exchange()

        ! local variables
        INTEGER(intk) :: i, iface, igrid
        INTEGER(intk) :: iface1, iface2, iface3
        INTEGER(intk) :: itypbc1, itypbc2, itypbc3
        INTEGER(intk) :: iprocnbr, itypbc, inbrgrid

        INTEGER(int32), ALLOCATABLE :: maxTag(:)
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)

        INTEGER(intk) :: neighbours(26)
        INTEGER :: iexchange

        CALL start_timer(900)
        CALL start_timer(910)

        ! Maximum number of connections for "simple" cases is number
        ! of grids*26. However, due to the possible prescence of
        ! precursors etc, we add a few more.
        maxConns = INT((nMyGrids+1)*26.0*1.2, intk)
        ALLOCATE(sendConns(4, maxConns))
        ALLOCATE(recvConns(4, maxConns))
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
        ! JULIUS: partly cleaned, but possibly more to improve

        ! <--------------------------

        DO i = 1, nMyGrids

            ! getting the grid parameters
            igrid = myGrids(i)
            CALL get_neighbours(neighbours, igrid)

            ! Check surfaces of grid
            DO iface = 1, 6

                itypbc = itypboconds(1, iface, igrid)
                inbrgrid = neighbours(iface)

                IF (itypbc == 7 .OR. itypbc == 19) THEN

                    IF (inbrgrid == 0) THEN
                        CYCLE
                    END IF

                    iprocnbr = idprocofgrd(neighbours(iface))

                    IF (iprocnbr == myid) THEN
                        CYCLE
                    END IF

                    ! only if neighbor not already listed
                    IF ( sendcounts(iprocnbr) == 0 ) THEN

                        iexchange = 1
                        nRecv = nRecv + 1
                        maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                        recvConns(1, nRecv) = myid              ! Receiving process (this process)
                        recvConns(2, nRecv) = iprocnbr          ! Sending process (neighbour process)
                        recvConns(3, nRecv) = maxTag(iprocnbr)  ! Message tag
                        recvConns(4, nRecv) = iexchange         ! Geometry exchange flag

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

                inbrgrid = neighbours(iface)

                IF (itypbc1 == 7 .OR. itypbc2 == 7 .OR. &
                    itypbc1 == 19 .OR. itypbc2 == 19) THEN

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

                        recvConns(1, nRecv) = myid              ! Receiving process (this process)
                        recvConns(2, nRecv) = iprocnbr          ! Sending process (neighbour process)
                        recvConns(3, nRecv) = maxTag(iprocnbr)  ! Message tag
                        recvConns(4, nRecv) = iexchange         ! Geometry exchange flag

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

                inbrgrid = neighbours(iface)

                IF (itypbc1 == 7 .OR. itypbc2 == 7 .OR. itypbc3 == 7 .OR. &
                    itypbc1 == 19 .OR. itypbc2 == 19 .OR. itypbc3 == 19) THEN

                    IF (inbrgrid == 0) THEN
                        CYCLE
                    END IF

                    iprocnbr = idprocofgrd(inbrgrid)

                    IF (iprocnbr == myid) THEN
                        CYCLE
                    END IF

                    ! only if neighbor not already listed
                    IF (sendcounts(iprocnbr) == 0) THEN

                        iexchange = 1
                        nRecv = nRecv + 1
                        maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                        recvConns(1, nRecv) = myid              ! Receiving process (this process)
                        recvConns(2, nRecv) = iprocnbr          ! Sending process (neighbour process)
                        recvConns(3, nRecv) = maxTag(iprocnbr)  ! Message tag
                        recvConns(4, nRecv) = iexchange         ! Geometry exchange flag

                        ! JULIUS: why not assign to recvcounts here? Naming seems confusing to me.
                        sendcounts(iprocnbr) = SIZE(recvConns, 1)   ! not an increment

                    END IF
                END IF
            END DO
        END DO

        ! <--------------------------

        ! Sort recvConns by process ID
        CALL sort_conns_unique(recvConns(:,1:nRecv))
        iRecv = nRecv

        ! JULIUS: whats the point the following (up to  CALL create_particle_mpitype)?
        ! Would sendConn(1,i) = recvCon(2,i) / sendConn(2,i) = recvCon(1,i) not suffice? And why is sendConn needed anyways if symmetric to recvConn?

        ! Calculate sdispl offset (send)
        DO i=1,numprocs-1
            ! = value is either 0 (not a neighbor) or 4 (a neighbor)
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Calculate rdispl offset (receive)
        DO i=1,numprocs-1
            ! = value is either 0 (not a neighbor) or 4 (a neighbor)
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

        ! only for debugging, doesnt even work reliable (?)
        DO i = 0, numprocs - 1

            IF (myid == i) THEN

                WRITE(*,*) 'I am proc:', myid
                WRITE(*,*) 'I own grids: '

                WRITE(*,*) mygrids(:)

                WRITE(*,*) ' - I receive from the following ', iRecv, 'processes (recvConns):'
                WRITE(*,*) recvConns(2, 1:iRecv)

                WRITE(*,*) ' - I send to the following ', iSend, 'processes (sendConns):'
                WRITE(*,*) sendConns(1, 1:iSend)

            END IF

            CALL MPI_Barrier(MPI_COMM_WORLD)

        END DO

        nRecv = 0

        ALLOCATE(npsend(iSend))
        npsend = 0

        ALLOCATE(ndispsend(iSend))
        ndispsend = 0

        ALLOCATE(nprecv(iRecv))
        nprecv = 0

        ALLOCATE(ndisprecv(iRecv))
        ndisprecv = 0

        ! creating the MPI data type
        CALL create_particle_mpitype(particle_mpitype)
        isInit = .TRUE.

        !DEALLOCATE(maxTag)
        !DEALLOCATE(sendcounts)
        !DEALLOCATE(sdispls)
        !DEALLOCATE(recvcounts)
        !DEALLOCATE(rdispls)

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_exchange

    SUBROUTINE finish_particle_exchange()

        isInit = .FALSE.

        DEALLOCATE(sendConns)
        DEALLOCATE(recvConns)
        DEALLOCATE(recvIdxList)
        DEALLOCATE(sendList)
        DEALLOCATE(recvList)
        DEALLOCATE(sendReqs)
        DEALLOCATE(recvReqs)
        DEALLOCATE(npsend)
        DEALLOCATE(ndispsend)
        DEALLOCATE(nprecv)
        DEALLOCATE(ndisprecv)

    END SUBROUTINE finish_particle_exchange

    SUBROUTINE sort_conns_unique(list)
        ! Input array to be sorted
        INTEGER(int32), INTENT(inout) :: list(:,:)

        INTEGER(intk) :: i, j

        ! Temporary storage
        INTEGER(int32) :: temp(4)

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

    SUBROUTINE get_target_grid(particle, destgrid, destproc, iface)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(in) :: particle
        INTEGER(intk), INTENT(out) :: destgrid
        INTEGER(intk), INTENT(out) :: destproc
        INTEGER(intk), INTENT(out) :: iface

        ! local variables
        INTEGER(intk) :: neighbours(26)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz, dist

        ! getting the box of last grid the particla
        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, particle%igrid)

        ! intialization (will be overwritten is particle left grid)
        iface = -1

        ! check if particle is still on grid first as this will be the case for most particles (assuming a reasonable grid size)
        ! to reduce operations
        CALL get_exit_face(particle, dist, iface)

        ! if the distance to the grid "dist" is 0, the particle is still on the grid
        ! however, get_exit_face might still return (iface > 0) if the particle is exactly on any grid boundary
        ! if so, set iface to 0
        IF (dist == 0) THEN
            iface = 0
        END IF

        IF (iface == 0) THEN
            ! particle stays on the same grid
            destgrid = particle%igrid
            destproc = particle%iproc

            SELECT CASE (TRIM(particle_terminal))
                CASE ("none")
                    CONTINUE
                CASE ("normal")
                    CONTINUE
                CASE ("verbose")
                    WRITE(*, '("Proc ", I0 ," Destination Proc: ", I0)') myid, destproc
                    WRITE(*, '("Proc ", I0 ," Destination Grid: ", I0)') myid, destgrid
                    IF (myid == 0) THEN
                        WRITE(*, *) " "
                    END IF
            END SELECT

            IF (destproc /= myid) THEN
                WRITE(*,*) 'Inconsistent particle parameters'
                CALL errr(__FILE__, __LINE__)
            END IF

        ELSE IF (iface > 0) THEN
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

            IF (destproc == myid) THEN
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

        CALL MPI_Get_address(foo%state, disp(1))
        ! JULIUS: isnt the following disp declaration unnessecary?
        CALL MPI_Get_address(foo%ipart, disp(2))
        CALL MPI_Get_address(foo%iproc, disp(3))
        CALL MPI_Get_address(foo%igrid, disp(4))
        CALL MPI_Get_address(foo%islice, disp(5))
        CALL MPI_Get_address(foo%gitstep, disp(6))
        CALL MPI_Get_address(foo%sitstep, disp(7))
        CALL MPI_Get_address(foo%ijkcell, disp(8))
        CALL MPI_Get_address(foo%x, disp(9))
        CALL MPI_Get_address(foo%y, disp(10))
        CALL MPI_Get_address(foo%z, disp(11))

        types(1) = mglet_mpi_int    ! state
        types(2) = mglet_mpi_int    ! ipart
        types(3) = mglet_mpi_int    ! iproc
        types(4) = mglet_mpi_int    ! igrid
        types(5) = mglet_mpi_int    ! islice
        types(6) = mglet_mpi_int    ! gitstep
        types(7) = mglet_mpi_int    ! sitstep
        types(8) = triple_int_mpi_type  ! ijkcell(3)
        types(9) = mglet_mpi_real     ! x
        types(10) = mglet_mpi_real    ! y
        types(11) = mglet_mpi_real    ! z

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

            DO i = i, sizeRecvBuf ! i = sizeSendBuf + 1
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
                    IF (particle_list%particles(particle_list%ifinal)%state >= 1) THEN

                        particle_list%particles(sendind(i)) = particle_list%particles(particle_list%ifinal)
                        particle_list%particles(particle_list%ifinal)%state = -1
                        particle_list%ifinal = particle_list%ifinal - 1
                        EXIT

                    ELSE

                        particle_list%ifinal = particle_list%ifinal - 1

                    END IF
                END DO

                IF (particle_list%particles(sendind(i))%state < 1) THEN
                    particle_list%ifinal = particle_list%ifinal - 1
                END IF

            END DO

        END IF

    END SUBROUTINE integrate_particles

END MODULE particle_exchange_mod