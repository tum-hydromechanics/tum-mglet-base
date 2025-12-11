MODULE particle_exchange_mod

    USE, INTRINSIC :: ISO_C_BINDING
    USE MPI_f08
    USE comms_mod

    USE particle_runtimestat_mod, ONLY: psim_n_sent
    USE particle_list_mod
    USE particle_statistics_mod
    USE particle_utils_mod
    USE particle_loadbalance_mod

    IMPLICIT NONE

    PRIVATE

    ! Maximum number of connections on one single process, either
    ! outgoing or incomming, on any single grid level
    INTEGER(intk) :: maxConns

    ! Lists that hold the symmetric (Send and Recv) connections on the particle level
    ! Dimensions contain:
    !   Dim 1: Information about a specific connection
    !   Dim 2: The different connections
    ! The information in the first dimension is sorted as follows:
    !   Field 1: Rank of this process
    !   Field 2: Rank of neighbour process
    !   Field 3: Message tag (for MPI)
    !   Field 4: Geometry exchange flag
    INTEGER(intk), ALLOCATABLE :: symConns(:,:)

    ! Lists that hold the send and receive request arrays
    TYPE(MPI_Request), ALLOCATABLE :: sendReqs(:), recvReqs(:)

    ! Lists that hold the messages that are ACTUALLY sendt and received
    INTEGER(intk) :: nConns

    INTEGER(intk), ALLOCATABLE :: nprecv(:)
    INTEGER(intk), ALLOCATABLE :: npsend(:)
    INTEGER(intk), ALLOCATABLE :: ndispsend(:)
    INTEGER(intk), ALLOCATABLE :: ndisprecv(:)

    ! MPI type for the particle
    TYPE(MPI_Datatype) :: particle_mpitype
    
    ! Particle type (not a class, as otherwise polymorphism implied)
    TYPE(baseparticle_t), ALLOCATABLE :: sendBufParticle(:)
    TYPE(baseparticle_t), ALLOCATABLE :: recvBufParticle(:)
    INTEGER(intk), ALLOCATABLE :: sendind(:)

    ! Sizes of the buffers
    INTEGER(intk) :: sizeSendBuf
    INTEGER(intk) :: sizeRecvBuf

    ! Variable to indicate if the connection information has been created
    LOGICAL :: isInit = .FALSE.

    LOGICAL :: high_mem_sorting = .FALSE.

    PUBLIC :: init_particle_exchange, exchange_particles, finish_particle_exchange, get_target_grid

CONTAINS

    SUBROUTINE exchange_particles(particle_list, ittot, itstep)

        IMPLICIT NONE

        ! subroutine argument
        TYPE(particle_list_t), INTENT(inout) :: particle_list
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: itstep

        !local variables
        INTEGER(intk) :: i, j, iproc, pos, num, dummy
        INTEGER(intk) :: destgrid, destproc, iface
        INTEGER(intk) :: iprocnbr, cSend, cRecv
        INTEGER(intk) :: active_np_old  ! for safety checks
        INTEGER(intk) :: err_local = 0, err_global = 0

        CALL start_timer(900)
        CALL start_timer(940)

        IF (.NOT. isInit) THEN
            WRITE(*,*) 'Particle connect not initialized'
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (MOD(ittot, loadbalance_step) == 0) CALL set_loadbalance_connections()

        active_np_old = particle_list%active_np

        npsend = 0
        nprecv = -1

        DO i = 1, particle_list%ifinal

            ! jumping inactive particles
            IF (particle_list%particles(i)%state < 1) THEN
                IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                        WRITE(*, '("WARNING on proc ", I0, ": Particle list entry ", I0, " unexpectately holds and inactive Partcle!")') myid, i
                END IF
                err_local = 1
                CYCLE
            END IF

            ! setting the destination of particle
            CALL get_target_grid(particle_list%particles(i), destgrid, destproc, iface)

            IF (destproc > numprocs .OR. destproc < 0) THEN
                WRITE(*,*) 'Obviously ill-addressed particle to proc', destproc
                CALL errr(__FILE__, __LINE__)
            END IF

            ! manipulation of relative coordinates at periodic boundaries
            ! at this point igrid is still NOT updated, meaning particle%igrid is still the "old" grid
            CALL update_coordinates(particle_list%particles(i), destgrid, iface)

            ! for particle slice statistics: must be called after update_coordinates !!!
            CALL stop_timer(940)
            CALL associate_new_slice(particle_list%particles(i), ittot, itstep)
            CALL start_timer(940)

            ! triage of particles
            IF (particle_list%particles(i)%igrid == destgrid) THEN

                ! particle stays on grid => only update cell
                CALL update_particle_cell(particle_list%particles(i))

            ELSE

                ! for particle statistics
                CALL stop_timer(940)
                CALL deregister_particle(particle_list%particles(i), ittot, itstep)
                CALL start_timer(940)

                ! particle changes the grid
                IF (destproc == myid) THEN

                    ! particle remains on process
                    particle_list%particles(i)%igrid = destgrid
                    CALL set_particle_cell(particle_list%particles(i))

                    ! for particle statistics
                    CALL stop_timer(940)
                    CALL register_particle(particle_list%particles(i), itstep)
                    CALL start_timer(940)

                ELSE

                    ! particle is marked for MPI transfer
                    particle_list%particles(i)%iproc = destproc
                    particle_list%particles(i)%igrid = destgrid

                    ! search for the process to send to (only checks "neighbor processes")
                    DO iproc = 1, nConns
                        IF (symConns(2, iproc) == particle_list%particles(i)%iproc) THEN
                            npsend(iproc) = npsend(iproc) + 1
                        END IF
                    END DO

                END IF

            END IF
        END DO

        IF (numprocs == 1) THEN
            CALL stop_timer(940)
            CALL stop_timer(900)
            RETURN
        END IF

        ! posting NON-blocking receives
        ! int MPI_Irecv(void *buf, int count,
        !     MPI_Datatype datatype, int source,
        !     int tag, MPI_Comm comm, MPI_Request *request)
        DO i = 1, nConns
            iprocnbr = symConns(2, i)
            CALL MPI_Irecv( nprecv(i), 1, mglet_mpi_int, &
            iprocnbr, 123, MPI_COMM_WORLD, recvreqs(i) )
        END DO

        ! posting non-blocking (!) sends
        ! int MPI_Isend(const void *buf, int count,
        !     MPI_Datatype datatype, int dest, int tag,
        !     MPI_Comm comm, MPI_Request *request)
        DO i = 1, nConns
            iprocnbr = symConns(2, i)
            CALL MPI_Isend( npsend(i), 1, mglet_mpi_int, &
            iprocnbr, 123, MPI_COMM_WORLD, sendreqs(i) )
        END DO

        ! displacements for start of section for one destination
        ndispsend = -1
        IF (SIZE(ndispsend) > 0) ndispsend(1) = 1
        DO i = 2, nConns
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

        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
            psim_n_sent = psim_n_sent + SUM(npsend)
        END IF

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

            ! buffer is filled
            DO iproc = 1, nConns
                IF ( symConns(2, iproc) == particle_list%particles(i)%iproc ) THEN

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
                    particle_list%particles(i)%ipart = -1
                    particle_list%particles(i)%state = -1
                    particle_list%active_np = particle_list%active_np - 1

                    ! collect indices of particles list entries that will be empty after MPI send
                    sendind(j) = i
                    j = j + 1

                END IF
            END DO

        END DO

        ! resetting after incrementation
        ndispsend = -1; ndispsend(1) = 1
        DO i = 2, nConns
            ndispsend(i) = ndispsend(i-1) + npsend(i-1)
        END DO

        ! buffer must be full with valid particles without gaps
        DO i = 1, sizeSendBuf
            IF ( sendBufParticle(i)%state < 1 ) THEN
                WRITE(*,*) 'Proc', myid, ': Invalid send buffer entry at ', i
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        ! checking if communication done (one call should suffice...)
        CALL MPI_Waitall(nConns, sendreqs, MPI_STATUSES_IGNORE)
        CALL MPI_Waitall(nConns, recvreqs, MPI_STATUSES_IGNORE)

        ! displacements for start of section for one source
        ndisprecv = -1;
        IF (SIZE(ndisprecv) > 0) ndisprecv(1) = 1
        DO i = 2, nConns
            IF ( nprecv(i-1) < 0 ) THEN
                WRITE(*,*) 'Invalid number of received particles'
                CALL errr(__FILE__, __LINE__)
            END IF
            ndisprecv(i) = ndisprecv(i-1) + nprecv(i-1)
        END DO

        sizeRecvBuf = SUM(nprecv)
        ALLOCATE(recvBufParticle(sizeRecvBuf))

        ! Check if list is long enough and add additional space if not
        IF (particle_list%max_np - particle_list%active_np < sizeRecvBuf) THEN
            CALL reallocate_particle_list(particle_list, INT(1.0 * (sizeRecvBuf - (particle_list%max_np - particle_list%active_np))))
        END IF

        ! posting NON-blocking receives
        ! int MPI_Irecv(void *buf, int count,
        !     MPI_Datatype datatype, int source,
        !     int tag, MPI_Comm comm, MPI_Request *request)
        cRecv = 0
        DO i = 1, nConns
            iprocnbr = symConns(2, i)
            pos = ndisprecv(i)
            num = nprecv(i)
            IF ( num > 0 ) THEN
                cRecv = cRecv + 1
                CALL MPI_Irecv( recvBufParticle(pos), num, particle_mpitype, &
                iprocnbr, 321, MPI_COMM_WORLD, recvreqs(cRecv) )
            END IF
        END DO

        ! posting NON-blocking sends
        ! int MPI_Isend(const void *buf, int count,
        !     MPI_Datatype datatype, int dest, int tag,
        !     MPI_Comm comm, MPI_Request *request)
        cSend = 0
        DO i = 1, nConns
            iprocnbr = symConns(2, i)
            pos = ndispsend(i)
            num = npsend(i)
            IF ( num > 0 ) THEN
                cSend = cSend + 1
                CALL MPI_Isend( sendBufParticle(pos), num, particle_mpitype, &
                iprocnbr, 321, MPI_COMM_WORLD, sendreqs(cSend) )
            END IF
        END DO

        ! checking if communication done (one call should suffice...)
        CALL MPI_Waitall(cSend, sendreqs, MPI_STATUSES_IGNORE)
        CALL MPI_Waitall(cRecv, recvreqs, MPI_STATUSES_IGNORE)

        ! some checks and assigning the new cell indices
        IF (sizeRecvBuf > 0) THEN
            DO i = 1, sizeRecvBuf
                ! check if correctly delivered
                IF (recvBufParticle(i)%state < 1) THEN
                    WRITE(*,*) "Inactive particle delivered"
                    CALL errr(__FILE__, __LINE__)
                END IF
                ! check if correctly delivered
                IF (recvBufParticle(i)%iproc /= myid) THEN
                    WRITE(*,*) "Particle delivered to wrong proc", i, recvBufParticle(i)%iproc, myid
                    CALL errr(__FILE__, __LINE__)
                END IF

                CALL set_particle_cell(recvBufParticle(i))

                ! for gridstat
                CALL stop_timer(940)
                CALL register_particle(recvBufParticle(i), itstep)
                CALL start_timer(940)

            END DO
        END IF

        ! Copy recieved particles into the list
        ! CAUTION: up to here, particle_list%particles(particle_list%ifinal)%state might be < 1 ("empty")
        IF (.NOT. dparticle_sorting) THEN
            CALL integrate_particles_unsorted(particle_list, sendind)
            CALL check_plist(particle_list, abort = .TRUE.)
        ELSE
            IF (.NOT. high_mem_sorting) THEN
                CALL integrate_particles_unsorted(particle_list, sendind)
                CALL check_plist(particle_list, abort = .TRUE.)

                CALL sort_by_grid(particle_list)
                CALL check_plist(particle_list, abort = .TRUE.)
            ELSE
                CALL integrate_particles_sorted(particle_list, sendind)
                CALL check_plist(particle_list, abort = .TRUE.)
            END IF
        END IF
    
        ! Some safety checks
        IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
            IF (myid /= 0) THEN
                CALL MPI_Recv(dummy, 1, mglet_mpi_int, myid - 1, 900, &
                MPI_COMM_WORLD, MPI_STATUS_IGNORE)
            END IF
            CALL print_list_status(particle_list)
            WRITE(*, '()')
            IF (myid /= numprocs - 1) THEN
                CALL MPI_Send(dummy, 1, mglet_mpi_int, myid + 1, 900, &
                MPI_COMM_WORLD)
            END IF
        END IF

        IF (particle_list%active_np < active_np_old + sizeRecvBuf - sizeSendBuf) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("WARNING on proc ", I0, ": Particle list holds FEWER particles than expected!")') myid
            END IF
            err_local = 1
        END IF

        IF (particle_list%active_np > active_np_old + sizeRecvBuf - sizeSendBuf) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("WARNING on proc ", I0, ": Particle list holds MORE particles than expected!")') myid
            END IF
            err_local = 1
        END IF

        IF (particle_list%ifinal /= particle_list%active_np) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("WARNING on proc ", I0, ": my_particle_list%active_np (", I0, ") does not coincide with my_particle_list%ifinal (", I0, ")" )') &
                 myid, particle_list%active_np, particle_list%ifinal
            END IF
            err_local = 1
        END IF

        ! TODO: make the following error gathering conditional for compilation as a debugging feature
        IF (TRIM(particle_terminal) == "verbose") THEN
            CALL MPI_Allreduce(err_local, err_global, 1, mglet_mpi_int, MPI_MAX, MPI_COMM_WORLD)
            IF (err_global == 0) THEN
                CALL write_particle_list_txt(ittot)
                CALL write_buffer(ittot, "Send")
                CALL write_buffer(ittot, "Recv")
            ELSE
                CALL write_particle_list_txt(ittot, "err")
                CALL write_buffer(ittot, "Send", "err")
                CALL write_buffer(ittot, "Recv", "err")
            END IF
            IF (err_global == 1) THEN
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

        IF (ALLOCATED(sendBufParticle)) DEALLOCATE(sendBufParticle)
        IF (ALLOCATED(recvBufParticle)) DEALLOCATE(recvBufParticle)
        IF (ALLOCATED(sendind)) DEALLOCATE(sendind)

        CALL stop_timer(940)
        CALL stop_timer(900)

    END SUBROUTINE exchange_particles

    SUBROUTINE init_particle_exchange()

        ! local variables
        INTEGER(intk) :: i, iface, igrid
        INTEGER(intk) :: iprocnbr, dummy

        INTEGER(int32), ALLOCATABLE :: maxTag(:)
        INTEGER(int32), ALLOCATABLE :: is_connected(:)

        INTEGER(intk) :: neighbour
        INTEGER :: iexchange

        CALL start_timer(900)
        CALL start_timer(910)

        ! Maximum number of connections for "simple" cases is number
        ! of grids*26. However, due to the possible prescence of
        ! precursors etc, we add a few more.
        maxConns = INT((nMyGrids+1)*26.0*1.2, intk)
        ALLOCATE(symConns(4, maxConns))
        symConns = 0

        ! The maximum number of concurrent communications are the number
        ! of processes
        ALLOCATE(sendReqs(numprocs))
        ALLOCATE(recvReqs(numprocs))

        ALLOCATE(maxTag(0:numprocs-1))
        ALLOCATE(is_connected(0:numprocs-1))
        maxTag = 0
        is_connected = 0
        nConns = 0

        ! ---------------------------

        DO i = 1, nmygridslvl(particle_level)

            ! getting the grid parameters
            igrid = mygridslvl(i, particle_level)

            ! Check surfaces of grid
            DO iface = 1, 26

                neighbour = particle_boundaries%face_neighbours(iface, igrid)

                IF (neighbour == igrid) THEN
                    CYCLE
                END IF

                iprocnbr = idprocofgrd(neighbour)

                IF (iprocnbr == myid) THEN
                    CYCLE
                END IF

                ! only if neighbor not already listed
                IF (is_connected(iprocnbr) == 0) THEN

                    iexchange = 1
                    nConns = nConns + 1
                    maxTag(iprocnbr) = maxTag(iprocnbr) + 1

                    symConns(1, nConns) = myid              ! Receiving process (this process)
                    symConns(2, nConns) = iprocnbr          ! Sending process (neighbour process)
                    symConns(3, nConns) = maxTag(iprocnbr)  ! Message tag
                    symConns(4, nConns) = iexchange         ! Geometry exchange flag

                    is_connected(iprocnbr) = 1

                END IF
            END DO
        END DO

        ! Sort symConns by process ID
        CALL sort_conns_unique(symConns(:,1:nConns), 2, .TRUE., myid)

        IF (TRIM(particle_terminal) == "verbose") THEN
            IF (myid /= 0) THEN
                CALL MPI_Recv(dummy, 1, mglet_mpi_int, myid - 1, 900, &
                MPI_COMM_WORLD, MPI_STATUS_IGNORE)
            END IF
            WRITE(*,*) 'I am proc:', myid
            WRITE(*,*) 'I own grids (on particle_level): '
            WRITE(*,*) mygridslvl(:, particle_level)
            WRITE(*,*) ' - I connect to the following ', nConns, 'processes (symConns):'
            WRITE(*,*) symConns(2, 1:nConns)
            WRITE(*, '()')
            IF (myid /= numprocs - 1) THEN
                CALL MPI_Send(dummy, 1, mglet_mpi_int, myid + 1, 900, &
                MPI_COMM_WORLD)
            END IF
        END IF

        ALLOCATE(npsend(nConns))
        npsend = 0

        ALLOCATE(ndispsend(nConns))
        ndispsend = 0

        ALLOCATE(nprecv(nConns))
        nprecv = 0

        ALLOCATE(ndisprecv(nConns))
        ndisprecv = 0

        ! creating the MPI data type
        CALL create_particle_mpitype(particle_mpitype)
        isInit = .TRUE.

        DEALLOCATE(maxTag)
        DEALLOCATE(is_connected)

        CALL stop_timer(910)
        CALL stop_timer(900)

    END SUBROUTINE init_particle_exchange


    SUBROUTINE finish_particle_exchange()

        CALL start_timer(900)
        CALL start_timer(990)
        isInit = .FALSE.

        IF (ALLOCATED(symConns)) DEALLOCATE(symConns)
        IF (ALLOCATED(sendReqs)) DEALLOCATE(sendReqs)
        IF (ALLOCATED(recvReqs)) DEALLOCATE(recvReqs)
        IF (ALLOCATED(npsend)) DEALLOCATE(npsend)
        IF (ALLOCATED(ndispsend)) DEALLOCATE(ndispsend)
        IF (ALLOCATED(nprecv)) DEALLOCATE(nprecv)
        IF (ALLOCATED(ndisprecv)) DEALLOCATE(ndisprecv)

        CALL stop_timer(990)
        CALL stop_timer(900)

    END SUBROUTINE finish_particle_exchange


    SUBROUTINE get_target_grid(particle, destgrid, destproc, iface)

        ! subroutine arguments
        TYPE(baseparticle_t), INTENT(in) :: particle
        INTEGER(intk), INTENT(out) :: destgrid
        INTEGER(intk), INTENT(out) :: destproc
        INTEGER(intk), INTENT(out) :: iface

        ! local variables
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
        IF (dist <= 0.0_realk) THEN
            iface = 0
        END IF

        IF (iface == 0) THEN
            ! particle stays on the same grid
            destgrid = particle%igrid
            destproc = particle%iproc

            IF (TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("Proc ", I0 ," Destination Proc: ", I0)') myid, destproc
                WRITE(*, '("Proc ", I0 ," Destination Grid: ", I0)') myid, destgrid
                IF (myid == 0) THEN
                    WRITE(*, *) " "
                END IF
            END IF

            IF (destproc /= myid) THEN
                WRITE(*,*) 'Inconsistent particle parameters'
                CALL errr(__FILE__, __LINE__)
            END IF

        ELSE IF (iface > 0) THEN
            ! particle moves across grid boundary
            destgrid = particle_boundaries%face_neighbours(iface, particle%igrid)
            destproc = idprocofgrd(destgrid)

            IF (TRIM(particle_terminal) == "verbose") THEN
                WRITE(*, '("Destination Proc: ", I0)') destproc
                WRITE(*, '("Destination grid: ", I0)') destgrid
                IF (myid == 0) THEN
                    WRITE(*, *) " "
                END IF
            END IF

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
        TYPE(MPI_Datatype) :: triple_real_mpi_type

        CALL MPI_Type_contiguous(3, mglet_mpi_int, triple_int_mpi_type)
        CALL MPI_Type_contiguous(3, mglet_mpi_real, triple_real_mpi_type)

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
        CALL MPI_Get_address(foo%xyz_abs, disp(12))
        CALL MPI_Get_address(foo%xyz_sentry, disp(13))

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
        types(12) = triple_real_mpi_type ! xyz_abs
        types(13) = triple_real_mpi_type ! xyt_sentry

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
    ! ifinal input not adapted yet
    SUBROUTINE integrate_particles_unsorted(particle_list, sendind)

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

            ! i = sizeRecvBuf + 1
            DO i = i, sizeSendBuf

                IF (particle_list%ifinal < sendind(i)) THEN
                    EXIT
                END IF

                IF (particle_list%ifinal == sendind(i)) THEN
                    particle_list%ifinal = particle_list%ifinal - 1
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

    END SUBROUTINE integrate_particles_unsorted

    SUBROUTINE integrate_particles_sorted(particle_list, sendind)

        ! subroutine argument
        TYPE(particle_list_t), INTENT(inout) :: particle_list
        INTEGER(intk), INTENT(in) :: sendind(sizeSendBuf)

        !local variables
        INTEGER(intk) :: i, j

        CALL errr(__FILE__, __LINE__)

    END SUBROUTINE integrate_particles_sorted

    ! for debugging
    SUBROUTINE write_buffer(ittot, btyp, suffix)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        CHARACTER(len = 4), INTENT(in) :: btyp ! "Send" or "Recv"
        CHARACTER(len = 3), INTENT(in), OPTIONAL :: suffix

        ! local varibales
        CHARACTER(len = mglet_filename_max) :: filename
        INTEGER(intk) :: unit, i
        LOGICAL :: exists

        IF (PRESENT(suffix)) THEN
            WRITE(filename,'(A, "Buffer-", I0, "-", A, ".txt")') btyp, myid, suffix
        ELSE
            WRITE(filename,'(A, "Buffer-", I0, ".txt")') btyp, myid
        END IF

        INQUIRE(file = TRIM(filename), exist = exists)

        IF (exists) THEN
            OPEN(newunit = unit, file = TRIM(filename), status = 'OLD', action = 'WRITE')
        ELSE
            OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')
        END IF

        WRITE(unit, '(A, "Buffer ", I0, " - Timestep ", I0)') btyp, myid, ittot
        WRITE(unit, '(" ")')
        WRITE(unit, '("PARTICLES")')

        IF (btyp == "Send") THEN
            DO i = 1, SIZE(sendBufParticle)
                    WRITE(unit, '("sendind = ", I0)') sendind(i)
                    WRITE(unit, '("ipart = ", I9, ", iproc", I3, ", igrid = ", I3, ", state = ", I3)') sendBufParticle(i)%ipart, &
                    sendBufParticle(i)%iproc, sendBufParticle(i)%igrid, sendBufParticle(i)%state
                    WRITE(unit, '("i/j/k cell :", 3I9)') sendBufParticle(i)%ijkcell(1), &
                    sendBufParticle(i)%ijkcell(2), sendBufParticle(i)%ijkcell(3)
                    WRITE(unit, '("x/y/z      :", 3F9.6)') sendBufParticle(i)%x, &
                    sendBufParticle(i)%y, sendBufParticle(i)%z
                    WRITE(unit, '(" ")')
            END DO
        ELSEIF (btyp == "Recv") THEN
            DO i = 1, SIZE(recvBufParticle)
                    WRITE(unit, '("ipart = ", I9, ", iproc", I3, ", igrid = ", I3, ", state = ", I3)') recvBufParticle(i)%ipart, &
                    recvBufParticle(i)%iproc, recvBufParticle(i)%igrid, recvBufParticle(i)%state
                    WRITE(unit, '("i/j/k cell :", 3I9)') recvBufParticle(i)%ijkcell(1), &
                    recvBufParticle(i)%ijkcell(2), recvBufParticle(i)%ijkcell(3)
                    WRITE(unit, '("x/y/z      :", 3F9.6)') recvBufParticle(i)%x, &
                    recvBufParticle(i)%y, recvBufParticle(i)%z
                    WRITE(unit, '(" ")')
            END DO
        ELSE
            WRITE(*, '("Unknown Particle Buffer Type")')
            CALL errr(__FILE__, __LINE__)
        END IF

        CLOSE(unit)

    END SUBROUTINE write_buffer

END MODULE particle_exchange_mod