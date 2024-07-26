MODULE particle_connect_mod

USE core_mod      ! provides:
USE particlecore_mod

IMPLICIT NONE

CONTAINS

    SUBROUTINE get_target_grid(particle, destgrid, destproc)

        ! subroutine arguments
        CLASS(baseparticle_t), INTENT(in) :: particle
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

        IF ( iface > 0 ) THEN
            ! particle moves across boundary
            CALL get_neighbours(neighbours, particle%igrid)
            destgrid = neighbours(iface)
            destproc = idprocofgrd(destgrid)
        ELSE
            ! particle stays on the same grid
            destgrid = particle%igrid
            destproc = particle%iproc
            IF ( destproc /= myid ) THEN
                WRITE(*,*) 'Inconsistent particle parameters'
                CALL errr(__FILE__, __LINE__)
            END IF
        END IF

    END SUBROUTINE get_target_grid


    SUBROUTINE init_particle_connect()
        INTEGER(intk) :: i, iface, igrid
        INTEGER(intk) :: iface1, iface2, iface3
        INTEGER(intk) :: itypbc1, itypbc2, itypbc3
        INTEGER(intk) :: iprocnbr, itypbc, inbrface, inbrgrid

        INTEGER(int32), ALLOCATABLE :: maxTag(:)
        INTEGER(int32), ALLOCATABLE :: sendcounts(:), sdispls(:)
        INTEGER(int32), ALLOCATABLE :: recvcounts(:), rdispls(:)

        INTEGER(intk) :: nFaceTot
        INTEGER(intk) :: nFaceGeom
        INTEGER(intk) :: nLineTot
        INTEGER(intk) :: nLineGeom
        INTEGER(intk) :: nCornerTot
        INTEGER(intk) :: nCornerGeom

        INTEGER(intk) :: neighbours(26)

        LOGICAL :: exchange
        INTEGER :: iexchange

        nFaceTot = 0
        nFaceGeom = 0
        nLineTot = 0
        nLineGeom = 0
        nCornerTot = 0
        nCornerGeom = 0

        CALL set_timer(150, "CONNECT2")

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

        ! It is really important that nplane = 2 also in preconnect
        nplane = 2
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
                    nFaceTot = nFaceTot + 1
                    iexchange = 1
                    nFaceGeom = nFaceGeom + 1
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

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) + SIZE(recvConns, 1)
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
                    nLineTot = nLineTot + 1
                    iexchange = 1
                    nLineGeom = nLineGeom + 1
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

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) + SIZE(recvConns, 1)
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
                    nCornerTot = nCornerTot + 1
                    iexchange = 1
                    nCornerGeom = nCornerGeom + 1
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

                    sendcounts(iprocnbr) = sendcounts(iprocnbr) + SIZE(recvConns, 1)
                END IF
            END DO
        END DO

        iRecv = nRecv

        ! Sort recvConns by process ID
        CALL sort_conns(recvConns(:,1:nRecv))

        ! Calculate sdispl offset
        DO i=1,numprocs-1
            sdispls(i) = sdispls(i-1) + sendcounts(i-1)
        END DO

        ! First exchange NUMBER OF ELEMENTS TO RECEIVE, to be able to
        ! calculate rdispls array
        CALL MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, &
            MPI_INTEGER, MPI_COMM_WORLD)

        ! Calculate rdispl offset
        DO i=1,numprocs-1
            rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        END DO

        ! Check that number of connections fit in array
        iSend = (rdispls(numprocs-1) + recvcounts(numprocs-1))/SIZE(sendConns, 1)
        IF (iSend > maxConns) THEN
            write(*,*) "Number of connections exceeded on process ", myid
            write(*,*) "maxConns =", maxConns, "nMyGrids =", nMyGrids, &
                "iSend = ", iSend
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Exchange connection information
        CALL MPI_Alltoallv(recvConns(1, 1), sendcounts, sdispls, MPI_INTEGER, &
            sendConns(1, 1), recvcounts, rdispls, MPI_INTEGER, &
            MPI_COMM_WORLD)

        isInit = .TRUE.

        ! Agglomerate statistics
        CALL MPI_Allreduce(MPI_IN_PLACE, nFaceTot, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nFaceGeom, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nLineTot, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nLineGeom, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nCornerTot, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, nCornerGeom, 1, mglet_mpi_int, &
            MPI_SUM, MPI_COMM_WORLD)

        IF (myid == 0) THEN
            WRITE(*, '("PARTCILE CONNECT STATISTICS:")')
            WRITE(*, '(4X, "Faces:         ", I7, 4X, "with geometry: ", I7)') &
                nFaceTot, nFaceGeom
            WRITE(*, '(4X, "Lines:         ", I7, 4X, "with geometry: ", I7)') &
                nLineTot, nLineGeom
            WRITE(*, '(4X, "Corners:       ", I7, 4X, "with geometry: ", I7)') &
                nCornerTot, nCornerGeom
            WRITE(*, '()')
        END IF

        nRecv = 0

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


    SUBROUTINE finish_connect2()
        isInit = .FALSE.

        DEALLOCATE(sendConns)
        DEALLOCATE(recvConns)

        DEALLOCATE(recvIdxList)
        DEALLOCATE(sendList)
        DEALLOCATE(recvList)
        DEALLOCATE(sendReqs)
        DEALLOCATE(recvReqs)
    END SUBROUTINE finish_connect2


    SUBROUTINE sort_conns(list)
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

    END SUBROUTINE sort_conns




























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

END MODULE particle_connect_mod