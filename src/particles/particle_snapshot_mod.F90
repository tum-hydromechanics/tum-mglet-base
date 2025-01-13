MODULE particle_snapshot_mod

    ! This module is responsible for:
    ! Writing particle positions at predifined timesteps as VTK files (particle snapshots).

    USE utils_mod

    USE particle_list_mod

    IMPLICIT NONE

    TYPE psnapshot_info_t

        INTEGER(intk) :: iproc
        INTEGER(intk) :: nprocs

        ! should depend on the domain lengths and be determined in init_psnapshots
        CHARACTER(len = format_char_len) :: coordinate_format

        ! particles to write in snapshots
        INTEGER(intk) :: global_np
        INTEGER(intk) :: pstep
        INTEGER(intk), ALLOCATABLE :: particle_ids(:)

        INTEGER(intk) :: nsnapshots
        ! stores the number of particles for each snapshot; will be allocated to length = nsnapshots
        INTEGER(intk), ALLOCATABLE :: nparticles(:)
        ! stores the integer of all timesteps for which a particle snapshot will be produced
        INTEGER(intk), ALLOCATABLE :: timesteps(:)
        ! stores the time for each snapshot; will be allocated to length = nsnapshots
        REAL(realk), ALLOCATABLE :: phtimes(:)

        ! stores the id (starting from 1 and rising contiguously) of the current snapshot
        INTEGER(intk) :: counter = 0

    END TYPE psnapshot_info_t

    TYPE(psnapshot_info_t) :: psnapshot_info

CONTAINS

    !------------------------------

    SUBROUTINE init_psnapshots(ittot, mtstep, dt)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        INTEGER(intk), INTENT(in) :: mtstep
        REAL(realk), INTENT(in) :: dt

        ! local variables
        INTEGER(intk) :: i

        CALL start_timer(900)
        CALL start_timer(960)

       !INQUIRE(directory = './Particle_Snapshots', exist = snapshots_exist)

       !IF (snapshots_exist) THEN
       !    WRITE(*, *) "ERROR: Directory Particle_Snaphots already exists. Terminating Process!"
       !    CALL errr(__FILE__, __LINE__)
       !END IF

        IF (myid == 0) THEN
            CALL create_directory("Particle_Snapshots") ! ! ! realtive to working directory ! ! !
        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

        ! metainfo
        psnapshot_info%iproc = myid
        psnapshot_info%nprocs = numprocs
        psnapshot_info%coordinate_format = vtk_float_format

        ! particle info
        ! the following assumes that all ids within [1, 2, ..., global_np - 1, global_np] are assigned to particles unambiguously
        IF (psnapshot_npart >= global_np .OR. psnapshot_npart == psnapshot_write_all_particles_tag) THEN
            psnapshot_info%global_np = global_np
            ALLOCATE(psnapshot_info%particle_ids(0))
        ELSEIF (psnapshot_npart > 1) THEN
            psnapshot_info%global_np = psnapshot_npart
            psnapshot_info%pstep = NINT(REAL(global_np / psnapshot_npart))
            ALLOCATE(psnapshot_info%particle_ids(psnapshot_info%global_np))
            psnapshot_info%particle_ids(1) = 1
            DO i = 2, SIZE(psnapshot_info%particle_ids)
                psnapshot_info%particle_ids(i) = MIN(1 + (i - 1) * psnapshot_info%pstep, global_np)
            END DO
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

        ! time info
        IF (psnapshot_step == 1) THEN
            psnapshot_info%nsnapshots = mtstep + 1_intk
        ELSE
            psnapshot_info%nsnapshots = CEILING(REAL(mtstep / psnapshot_step) - 100 * EPSILON(1.0_realk), intk) + 1_intk
        END IF

        ALLOCATE(psnapshot_info%nparticles(psnapshot_info%nsnapshots))
        ALLOCATE(psnapshot_info%timesteps(psnapshot_info%nsnapshots))
        ALLOCATE(psnapshot_info%phtimes(psnapshot_info%nsnapshots))

        psnapshot_info%timesteps(1) = ittot
        psnapshot_info%timesteps(psnapshot_info%nsnapshots) = ittot + mtstep

        DO i = 2, psnapshot_info%nsnapshots - 1
            psnapshot_info%timesteps(i) = psnapshot_info%timesteps(i - 1) + psnapshot_step
        END DO

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*,*) "Writing Particle Snapshots with ", psnapshot_info%global_np, " particles at timesteps: "
                WRITE(*,*) ' '
                WRITE(*, '(I0, " + ")', advance="yes") ittot
                DO i = 2, psnapshot_info%nsnapshots
                    IF (MOD(i - 1_intk, 10) == 0) THEN
                        WRITE(*, '(I0, ",")', advance="yes") psnapshot_info%timesteps(i) - ittot
                    ELSE
                        WRITE(*, '(I0, ",")', advance="no") psnapshot_info%timesteps(i) - ittot
                        WRITE(*, '(A)', advance="no") ' '
                    END IF

                END DO
                WRITE(*, '()')
                WRITE(*, '()')
            END IF
        END IF

        CALL stop_timer(960)
        CALL stop_timer(900)

    END SUBROUTINE init_psnapshots

    !------------------------------

    SUBROUTINE write_psnapshot(ittot, timeph)

        ! subroutine arguments
        INTEGER(intk), INTENT(in) :: ittot
        REAL(realk), INTENT(in) :: timeph

        CALL start_timer(900)
        CALL start_timer(960)

        IF (psnapshot_info%timesteps(psnapshot_info%counter + 1) == ittot) THEN

            psnapshot_info%counter = psnapshot_info%counter + 1_intk

            ! Racing condition is prevented within the create_psnapshot_subfolder routine
            CALL create_psnapshot_subfolder()

            CALL write_psnapshot_piece()

            CALL write_psnapshot_master(timeph)

        END IF

        CALL stop_timer(960)
        CALL stop_timer(900)

    END SUBROUTINE write_psnapshot

    !------------------------------

    SUBROUTINE create_psnapshot_subfolder()

        ! local variables
        INTEGER(intk) :: i, dummy
        CHARACTER(len = mglet_filename_max) :: subfolder
        TYPE(MPI_Request) :: request

        IF (myid == 0) THEN

            WRITE(subfolder, '("Particle_Snapshots/snapshot", I0)') psnapshot_info%timesteps(psnapshot_info%counter)
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
        INTEGER(intk) :: i, j, jstart, counter, unit
        CHARACTER(len = mglet_filename_max) :: subfolder, filename, np_char
        LOGICAL, ALLOCATABLE :: dwrite_particle(:)

        WRITE(subfolder, '("Particle_Snapshots/snapshot", I0)') psnapshot_info%timesteps(psnapshot_info%counter)

        WRITE(filename, '(A, "/piece", I0, ".vtp")') TRIM(subfolder), myid

        ALLOCATE(dwrite_particle(my_particle_list%ifinal))
        IF (psnapshot_info%global_np /= global_np) THEN
            dwrite_particle = .FALSE.
            counter = 0
            DO i = 1, my_particle_list%ifinal
                jstart = FLOOR(REAL((my_particle_list%particles(i)%ipart - 1)) / psnapshot_info%pstep) + 1
                DO j = jstart, SIZE(psnapshot_info%particle_ids)
                    IF (psnapshot_info%particle_ids(j) == my_particle_list%particles(i)%ipart) THEN
                        dwrite_particle(i) = .TRUE.
                        counter = counter + 1
                    ELSEIF (psnapshot_info%particle_ids(j) > my_particle_list%particles(i)%ipart) THEN
                        EXIT
                    END IF
                END DO
            END DO
        ELSE
            dwrite_particle = .TRUE.
            counter = my_particle_list%active_np
        END IF

        WRITE(np_char, '(I0)') counter !my_particle_list%active_np

        OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

        WRITE(unit, '(A)') '<?xml version="1.0"?>'
        WRITE(unit, '(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
        WRITE(unit, '(A)') '  <PolyData>'
        WRITE(unit, '(A)') '    <Piece NumberOfPoints="' // TRIM(np_char) // &
                              '" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
        WRITE(unit, '(A)') '      <PointData Name="particle_id">'
        WRITE(unit, '(A)') '        <DataArray type="Int32" format="ascii" NumberOfComponents="1" Name="particle_id">'

        DO i = 1, my_particle_list%ifinal

            IF ( my_particle_list%particles(i)%state < 1 ) THEN
                CYCLE
            END IF

            IF (dwrite_particle(i)) THEN
                WRITE(unit, '("          ")', advance="no")
                WRITE(unit, '(I0)') my_particle_list%particles(i)%ipart
            END IF

        END DO

        WRITE(unit, '(A)') '        </DataArray>'
        WRITE(unit, '(A)') '      </PointData>'
        WRITE(unit, '(A)') '      <Points>'
        WRITE(unit, '(A)') '        <DataArray type="Float32" NumberOfComponents="3">'

        DO i = 1, my_particle_list%ifinal

            IF ( my_particle_list%particles(i)%state < 1 ) THEN
                CYCLE
            END IF

            IF (dwrite_particle(i)) THEN
                WRITE(unit, '("        ")', advance="no")
                WRITE(unit, psnapshot_info%coordinate_format, advance="no") my_particle_list%particles(i)%x
                WRITE(unit, '(A)', advance="no") ' '
                WRITE(unit, psnapshot_info%coordinate_format, advance="no") my_particle_list%particles(i)%y
                WRITE(unit, '(A)', advance="no") ' '
                WRITE(unit, psnapshot_info%coordinate_format, advance="yes") my_particle_list%particles(i)%z
            END IF

        END DO

        WRITE(unit, '(A)') '        </DataArray>'
        WRITE(unit, '(A)') '      </Points>'
        WRITE(unit, '(A)') '    </Piece>'
        WRITE(unit, '(A)') '  </PolyData>'
        WRITE(unit, '(A)') '</VTKFile>'

        CLOSE(unit)

    END SUBROUTINE write_psnapshot_piece

    !------------------------------

    ! writes pvtp (xml) file (master file for all pieces of the same time)
    SUBROUTINE write_psnapshot_master(timeph)

        !subroutine arguments
        REAL(realk), INTENT(in) :: timeph

        !local variables
        INTEGER(intk) :: proc, unit
        CHARACTER(len = mglet_filename_max) :: filename, piece

        psnapshot_info%phtimes(psnapshot_info%counter) = timeph

        IF (myid == 0) THEN

            WRITE(filename,'("Particle_Snapshots/snapshot", I0, ".pvtp")') psnapshot_info%timesteps(psnapshot_info%counter)

            OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

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

                WRITE(piece, '("snapshot", I0, "/piece", I0, ".vtp")') psnapshot_info%timesteps(psnapshot_info%counter), proc
                WRITE(unit, '(A)') '    <Piece Source="' // TRIM(piece) // '"/>'

            END DO

            WRITE(unit, '(A)') '  </PPolyData>'
            WRITE(unit, '(A)') '</VTKFile>'

            CLOSE(unit)

        END IF

    END SUBROUTINE write_psnapshot_master

    !------------------------------

    SUBROUTINE write_psnapshot_timeinfo()

        INTEGER(intk) :: i, unit
        CHARACTER(len = mglet_filename_max) :: filename

        CALL start_timer(900)
        CALL start_timer(960)

        IF (myid == 0) THEN

            WRITE(filename,'("Particle_Snapshots/timeinfo.txt")')

            OPEN(newunit = unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

            WRITE(unit, '("Number of timesteps: ", I0)') psnapshot_info%nsnapshots

            WRITE(unit, '("-------- snapshot times: ---------")')

            DO i = 1, psnapshot_info%nsnapshots

                WRITE(unit, '("timestep ", I0, ": time = ")', advance = "no") psnapshot_info%timesteps(i)
                WRITE(unit, psnapshot_info%coordinate_format) psnapshot_info%phtimes(i)

            END DO

            CLOSE(unit)

        END IF

        CALL stop_timer(960)
        CALL stop_timer(900)

    END SUBROUTINE write_psnapshot_timeinfo

    SUBROUTINE finish_particle_snapshots()

        IF (.NOT. dwrite_psnapshots) RETURN

        CALL start_timer(900)
        CALL start_timer(960)

        IF (ALLOCATED(psnapshot_info%particle_ids)) DEALLOCATE(psnapshot_info%particle_ids)
        IF (ALLOCATED(psnapshot_info%nparticles)) DEALLOCATE(psnapshot_info%nparticles)
        IF (ALLOCATED(psnapshot_info%timesteps)) DEALLOCATE(psnapshot_info%timesteps)
        IF (ALLOCATED(psnapshot_info%phtimes)) DEALLOCATE(psnapshot_info%phtimes)

        CALL stop_timer(960)
        CALL stop_timer(900)

    END SUBROUTINE finish_particle_snapshots

END MODULE particle_snapshot_mod
