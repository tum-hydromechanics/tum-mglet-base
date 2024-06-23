! a module for particle input/output operations, if initial condition (ic)
! or boundary conditions (bc) should be read from a file or information should
! be written to a file

MODULE particle_io_mod

    ! TODO: psnapshot_info_t to store number of serial files per time, times and accuracy/format of point coordinates
    ! (sort of as a log file)

    USE particle_list_mod
    USE MPI_f08


    IMPLICIT NONE

    TYPE psnapshot_info_t

        INTEGER(intk) :: iproc
        INTEGER(intk) :: nprocs

        INTEGER(intk) :: itstep = 8
        CHARACTER(6) :: coordinate_format = '(F6.2)' !should depend on the domain lengths and be determined in init_psnapshots

        INTEGER(intk) :: nsnapshots
        INTEGER(intk), ALLOCATABLE :: nparticles(:) ! stores the number of particles for each snapshot; will be allocated to length = nsnapshots
        INTEGER(intk), ALLOCATABLE :: timesteps(:) ! stores the integer of all timesteps for which a particle snapshot will be produced
        REAL(realk), ALLOCATABLE :: times(:) ! stores the time for each snapshot; will be allocated to length = nsnapshots

        INTEGER(intk) :: counter = 0_intk ! stores the id (starting from 1 and rising contiguously) of the current snapshot

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

        IF (myid == 0) THEN
            CALL create_directory("PARTICLE_SNAPSHOTS") ! ! ! realtive to working directory ! ! !
        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

        psnapshot_info%iproc = myid
        psnapshot_info%nprocs = numprocs

        psnapshot_info%nsnapshots = mtstep / psnapshot_info%itstep + 1_intk

        ALLOCATE(psnapshot_info%nparticles(psnapshot_info%nsnapshots))
        ALLOCATE(psnapshot_info%timesteps(psnapshot_info%nsnapshots))
        ALLOCATE(psnapshot_info%times(psnapshot_info%nsnapshots))

        psnapshot_info%timesteps(1) = 1_intk
        psnapshot_info%timesteps(psnapshot_info%nsnapshots) = mtstep


        DO i = 2, psnapshot_info%nsnapshots - 1

            psnapshot_info%timesteps(i) = (i - 1_intk) * psnapshot_info%itstep ! assumes that dt is constant for the whole run

        END DO

    END SUBROUTINE init_psnapshots

    !------------------------------

    SUBROUTINE write_psnapshots(itstep, timeph)

        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: timeph

        IF (psnapshot_info%timesteps(psnapshot_info%counter + 1) == itstep) THEN

            psnapshot_info%counter = psnapshot_info%counter + 1_intk ! should be broadcasted via MPI so that no missmatches occur ?

            CALL init_psnapshot_subfolder()

            CALL write_psnapshot_piece()

            CALL write_psnapshot_master(timeph)

        END IF

    END SUBROUTINE write_psnapshots

    !------------------------------

    SUBROUTINE init_psnapshot_subfolder()

        CHARACTER(len = mglet_filename_max) :: subfolder

        IF (myid == 0) THEN

            WRITE(subfolder, '("PARTICLE_SNAPSHOTS/snapshot", I0)') psnapshot_info%counter
            CALL create_directory(TRIM(subfolder)) ! ! ! realtive to working directory ! ! !

        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

    END SUBROUTINE init_psnapshot_subfolder

    !------------------------------

    SUBROUTINE write_psnapshot_piece()
    ! writes vtp (xml) file containing all particles of the respective process

        ! local variables
        INTEGER(intk) :: i, unit
        CHARACTER(len = mglet_filename_max) :: subfolder, filename, active_np_char
        !CHARACTER(:), ALLOCATABLE :: active_np_char

        WRITE(subfolder, '("PARTICLE_SNAPSHOTS/snapshot", I0)') psnapshot_info%counter

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

            IF (.NOT. my_particle_list%particle_stored(i)) THEN
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

            IF (.NOT. my_particle_list%particle_stored(i)) THEN
                CYCLE
            END IF

            WRITE(unit, '("          ")', advance="no")
            WRITE(unit, psnapshot_info%coordinate_format, advance="no") my_particle_list%particles(i)%x
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, psnapshot_info%coordinate_format, advance="no") my_particle_list%particles(i)%y
            WRITE(unit, '(A)', advance="no") ' '
            WRITE(unit, psnapshot_info%coordinate_format, advance="yes") my_particle_list%particles(i)%z

        END DO

        ! also write vertices for visualization?

        WRITE(unit, '(A)') '        </DataArray>'
        WRITE(unit, '(A)') '      </Points>'
        WRITE(unit, '(A)') '    </Piece>'
        WRITE(unit, '(A)') '  </PolyData>'
        WRITE(unit, '(A)') '</VTKFile>'

    END SUBROUTINE write_psnapshot_piece

    !------------------------------

    SUBROUTINE write_psnapshot_master(timeph)
    ! writes pvtp (xml) file (master file for all pieces of the same time)

        !subroutine arguments...
        REAL(realk), INTENT(in) :: timeph

        !local variables...
        INTEGER(intk) :: proc, unit
        CHARACTER(len = mglet_filename_max) :: filename, piece

        psnapshot_info%times(psnapshot_info%counter) = timeph

        IF (myid == 0) THEN

            WRITE(filename,'("PARTICLE_SNAPSHOTS/snapshot", I0, ".pvtp")') psnapshot_info%counter

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

                WRITE(piece, '("snapshot", I0, "/piece", I0, ".vtp")') psnapshot_info%counter, proc
                WRITE(unit, '(A)') '    <Piece Source="' // TRIM(piece) // '"/>'

            END DO

            WRITE(unit, '(A)') '  </PPolyData>'
            WRITE(unit, '(A)') '</VTKFile>'

        END IF

    END SUBROUTINE write_psnapshot_master

    !------------------------------

    SUBROUTINE write_psnapshot_timeinfo()

        INTEGER(intk) :: i, unit
        CHARACTER(len = mglet_filename_max) :: filename

        IF (myid == 0) THEN

            WRITE(filename,'("PARTICLE_SNAPSHOTS/timeinfo.txt")')

            OPEN(unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

            WRITE(unit, '("Number of times: ", I0)') psnapshot_info%nsnapshots

            WRITE(unit, '("-------- snapshot times: ---------")')

            DO i = 1, psnapshot_info%nsnapshots

                WRITE(unit, '("timestep ", I0, ": time = ")', advance = "no") i
                WRITE(unit, psnapshot_info%coordinate_format) psnapshot_info%times(i)

            END DO

        END IF

    END SUBROUTINE write_psnapshot_timeinfo

    !------------------------------

    ! TODO: particle trajectories

END MODULE
