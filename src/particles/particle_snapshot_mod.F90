MODULE particle_snapshot_mod

    ! This module is responsible for:
    ! Writing particle positions at predifined timesteps as VTK files (particle snapshots).

    USE utils_mod

    USE particle_list_mod

    IMPLICIT NONE

    TYPE psnapshot_info_t

        INTEGER(intk) :: iproc
        INTEGER(intk) :: nprocs

        INTEGER(intk) :: itstep
        !should depend on the domain lengths and be determined in init_psnapshots
        CHARACTER(7) :: coordinate_format = '(F12.6)'

        INTEGER(intk) :: nsnapshots
        ! stores the number of particles for each snapshot; will be allocated to length = nsnapshots
        INTEGER(intk), ALLOCATABLE :: nparticles(:)
        ! stores the integer of all timesteps for which a particle snapshot will be produced
        INTEGER(intk), ALLOCATABLE :: timesteps(:)
        ! stores the time for each snapshot; will be allocated to length = nsnapshots
        REAL(realk), ALLOCATABLE :: times(:)

        ! stores the id (starting from 1 and rising contiguously) of the current snapshot
        INTEGER(intk) :: counter = 0_intk

    END TYPE psnapshot_info_t

    TYPE(psnapshot_info_t) :: psnapshot_info

CONTAINS

    !------------------------------

    SUBROUTINE init_psnapshots(mtstep, dt)

        ! subroutine arguments ...
        INTEGER(intk) :: mtstep
        REAL(realk) :: dt

        ! local variables ...
        INTEGER(intk) :: i
        LOGICAL :: snapshots_exist

        CALL start_timer(920)

        INQUIRE(file = 'Particle_Snapshots/snapshot0.pvtp', exist = snapshots_exist)

        IF (snapshots_exist) THEN

            IF (myid == 0) THEN
                SELECT CASE (TRIM(particle_terminal))
                    CASE ("none")
                        CONTINUE
                    CASE ("normal")
                        WRITE(*, *) ' '
                        WRITE(*, *) "ERROR: Directory Particle_Snaphots already exists. Terminating Process!"
                    CASE ("verbose")
                        WRITE(*, *) ' '
                        WRITE(*, *) "ERROR: Directory Particle_Snaphots already exists. Terminating Process!"
                END SELECT
            END IF

            CALL errr(__FILE__, __LINE__)

        END IF

        IF (myid == 0) THEN
            CALL create_directory("Particle_Snapshots") ! ! ! realtive to working directory ! ! !
        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

        psnapshot_info%iproc = myid
        psnapshot_info%nprocs = numprocs
        psnapshot_info%itstep = psnapshot_step

        IF (psnapshot_info%itstep == 1_intk) THEN

            psnapshot_info%nsnapshots = mtstep + 1

        ELSE

            psnapshot_info%nsnapshots = CEILING(REAL(mtstep / psnapshot_info%itstep), intk) + 1_intk

        END IF

        ALLOCATE(psnapshot_info%nparticles(psnapshot_info%nsnapshots))
        ALLOCATE(psnapshot_info%timesteps(psnapshot_info%nsnapshots))
        ALLOCATE(psnapshot_info%times(psnapshot_info%nsnapshots))

        psnapshot_info%timesteps(1) = 0_intk
        psnapshot_info%timesteps(psnapshot_info%nsnapshots) = mtstep


        DO i = 2, psnapshot_info%nsnapshots - 1

            ! assumes that dt is constant for the whole run
            psnapshot_info%timesteps(i) = (i - 1_intk) * psnapshot_info%itstep

        END DO

        IF (myid == 0) THEN
            WRITE(*,*) 'Writing Particle Snapshots for timesteps: '
            WRITE(*,*) ' '

            DO i = 2, psnapshot_info%nsnapshots

                IF (MOD(i - 1_intk, 10) == 0) THEN

                    WRITE(*, '(I0)', advance="yes") psnapshot_info%timesteps(i)

                ELSE

                    WRITE(*, '(I0)', advance="no") psnapshot_info%timesteps(i)
                    WRITE(*, '(A)', advance="no") ' '

                END IF

            END DO
            WRITE(*,*) ' '
        END IF

        CALL stop_timer(920)

    END SUBROUTINE init_psnapshots

    !------------------------------

    SUBROUTINE write_psnapshot(itstep, timeph)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: timeph

        CALL start_timer(920)

        IF (psnapshot_info%timesteps(psnapshot_info%counter + 1) == itstep) THEN

            psnapshot_info%counter = psnapshot_info%counter + 1_intk

            ! Racing condition is prevented within the create_psnapshot_subfolder routine
            CALL create_psnapshot_subfolder()

            CALL write_psnapshot_piece()

            CALL write_psnapshot_master(timeph)

        END IF

        CALL stop_timer(920)

    END SUBROUTINE write_psnapshot

    !------------------------------

    SUBROUTINE create_psnapshot_subfolder()

        ! local variables
        INTEGER(intk) :: i, dummy
        CHARACTER(len = mglet_filename_max) :: subfolder
        TYPE(MPI_Request) :: request

        IF (myid == 0) THEN

            WRITE(subfolder, '("Particle_Snapshots/snapshot", I0)') psnapshot_info%counter - 1_intk
            CALL create_directory(TRIM(subfolder)) ! ! ! realtive to working directory ! ! !

            DO i = 1, numprocs - 1
                CALL MPI_ISend(dummy, 1, mglet_mpi_int, &
                i, 0, MPI_COMM_WORLD, request)
            END DO

        ELSE

            CALL MPI_Recv(dummy, 1, mglet_mpi_int, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE)

        END IF

        ! alternatively
        ! CALL MPI_Barrier(MPI_COMM_WORLD)

    END SUBROUTINE create_psnapshot_subfolder

    !------------------------------

    ! writes vtp (xml) file containing all particles of the respective process
    SUBROUTINE write_psnapshot_piece()

        ! local variables
        INTEGER(intk) :: i, unit = 162
        CHARACTER(len = mglet_filename_max) :: subfolder, filename, active_np_char
        !CHARACTER(:), ALLOCATABLE :: active_np_char

        WRITE(subfolder, '("Particle_Snapshots/snapshot", I0)') psnapshot_info%counter - 1_intk

        WRITE(filename, '(A, "/piece", I0, ".vtp")') TRIM(subfolder), myid

        !ALLOCATE(CHARACTER(CEILING(LOG10(REAL(my_particle_list%max_np))) :: active_np_char)
        WRITE(active_np_char, '(I0)') my_particle_list%active_np

        OPEN(unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

        WRITE(unit, '(A)') '<?xml version="1.0"?>'
        WRITE(unit, '(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
        WRITE(unit, '(A)') '  <PolyData>'
        WRITE(unit, '(A)') '    <Piece NumberOfPoints="' // TRIM(active_np_char) // &
                              '" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
        WRITE(unit, '(A)') '      <PointData Name="particle_id">'
        WRITE(unit, '(A)') '        <DataArray type="Int32" format="ascii" NumberOfComponents="1" Name="particle_id">'

        DO i = 1, my_particle_list%ifinal

            IF ( my_particle_list%particles(i)%state < 1 ) THEN
                CYCLE
            END IF

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, '(I0)') my_particle_list%particles(i)%ipart

        END DO

        WRITE(unit, '(A)') '        </DataArray>'
        WRITE(unit, '(A)') '      </PointData>'
        WRITE(unit, '(A)') '      <Points>'
        WRITE(unit, '(A)') '        <DataArray type="Float32" NumberOfComponents="3">'

        DO i = 1, my_particle_list%ifinal

            IF ( my_particle_list%particles(i)%state < 1 ) THEN
                CYCLE
            END IF

            WRITE(unit, '("        ")', advance="no")
            WRITE(unit, psnapshot_info%coordinate_format, advance="no") my_particle_list%particles(i)%x
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, psnapshot_info%coordinate_format, advance="no") my_particle_list%particles(i)%y
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, psnapshot_info%coordinate_format, advance="yes") my_particle_list%particles(i)%z

        END DO

        WRITE(unit, '(A)') '        </DataArray>'
        WRITE(unit, '(A)') '      </Points>'
        WRITE(unit, '(A)') '    </Piece>'
        WRITE(unit, '(A)') '  </PolyData>'
        WRITE(unit, '(A)') '</VTKFile>'

    END SUBROUTINE write_psnapshot_piece

    !------------------------------

    ! writes pvtp (xml) file (master file for all pieces of the same time)
    SUBROUTINE write_psnapshot_master(timeph)

        !subroutine arguments
        REAL(realk), INTENT(in) :: timeph

        !local variables
        INTEGER(intk) :: proc, unit = 163
        CHARACTER(len = mglet_filename_max) :: filename, piece

        psnapshot_info%times(psnapshot_info%counter) = timeph

        IF (myid == 0) THEN

            WRITE(filename,'("Particle_Snapshots/snapshot", I0, ".pvtp")') (psnapshot_info%counter - 1_intk)

            OPEN(unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

            WRITE(unit, '(A)') '<?xml version="1.0"?>'
            WRITE(unit, '(A)') '<VTKFile type="PPolyData" version="0.1" byte_order="LittleEndian">'
            WRITE(unit, '(A)') '  <PPolyData>'
            WRITE(unit, '(A)') '    <PPointData Name="particle_id">'
            WRITE(unit, '(A)') '       <PDataArray type="Int32" format="ascii" NumberOfComponents="1" Name="particle_id"/>'
            WRITE(unit, '(A)') '    </PPointData>'
            WRITE(unit, '(A)') '    <PPoints>'
            WRITE(unit, '(A)') '      <PDataArray type="Float32" format="ascii" NumberOfComponents="3"/>'
            WRITE(unit, '(A)') '    </PPoints>'

            DO proc = 0, numprocs - 1

                WRITE(piece, '("snapshot", I0, "/piece", I0, ".vtp")') (psnapshot_info%counter - 1_intk), proc
                WRITE(unit, '(A)') '    <Piece Source="' // TRIM(piece) // '"/>'

            END DO

            WRITE(unit, '(A)') '  </PPolyData>'
            WRITE(unit, '(A)') '</VTKFile>'

        END IF

    END SUBROUTINE write_psnapshot_master

    !------------------------------

    SUBROUTINE write_psnapshot_timeinfo()

        INTEGER(intk) :: i, unit = 164
        CHARACTER(len = mglet_filename_max) :: filename

        CALL start_timer(920)

        IF (myid == 0) THEN

            WRITE(filename,'("Particle_Snapshots/timeinfo.txt")')

            OPEN(unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

            WRITE(unit, '("Number of timesteps: ", I0)') psnapshot_info%nsnapshots

            WRITE(unit, '("-------- snapshot times: ---------")')

            DO i = 1, psnapshot_info%nsnapshots

                WRITE(unit, '("timestep ", I0, ": time = ")', advance = "no") i - 1
                WRITE(unit, psnapshot_info%coordinate_format) psnapshot_info%times(i)

            END DO

        END IF

        CALL stop_timer(920)

    END SUBROUTINE write_psnapshot_timeinfo

    !------------------------------

    ! TODO: particle trajectories ?

END MODULE particle_snapshot_mod
