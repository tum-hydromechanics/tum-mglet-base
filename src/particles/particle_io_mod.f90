! a module for particle input/output operations, if initial condition (ic) 
! or boundary conditions (bc) should be read from a file or information should 
! be written to a file 

MODULE particle_io_mod

! TODO: psnapshot_info_t to store number of serial files per time, times and accuracy/format of point coordinates
! (sort of as a log file) 

USE pvtp_mod
USE utils_mod
USE core_mod
USE timekeeper_mod
USE baseparticle_mod
USE particle_list_mod

!===================================

IMPLICIT NONE

TYPE psnapshot_info_t

	INTEGER(intk) :: iproc = myid
	INTEGER(intk) :: nprocs = numprocs

	REAL(realk) :: dt_psnapshot
	CHARACTER(6) :: coordinate_format = '(F6.2)' !should depend on the domain lengths and be determined in init_psnapshots
	INTEGER(intk) :: ntimesteps
	REAL(realk), ALLOCATABLE :: timesteps(:) ! stores the time for each snapshot; will be allocated to length = ntimesteps
	INTEGER, ALLOCATABLE :: nparticles(:) ! stores the number of particles for each snapshot; will be allocated to length = ntimesteps

	INTEGER(intk) :: counter = 0 ! stores the id (starting from 1 and rising contiguously) of the current snapshot

END TYPE psnapshot_info_t

!===================================

TYPE(psnapshot_info_t) :: psnapshot_info

!===================================

CONTAINS

	SUBROUTINE read_particle_ic()

	END SUBROUTINE read_particle_ic

	!------------------------------

	SUBROUTINE read_particle_bc()

	END SUBROUTINE read_particle_bc

	!------------------------------

	SUBROUTINE init_psnapshots()

	    IF (myid == 0) THEN
            CALL create_directory("PARTICLE_SNAPSHOTS") ! ! ! realtive to working directory ! ! !
        END IF
        
        CALL MPI_Barrier(MPI_COMM_WORLD)

        psnapshot_info%ntimesteps = CEILING((tend - timeph) / psnapshot_info%dt_psnapshot) + 1

        ALLOCATE(psnapshot_info%timesteps(psnapshot_info%ntimesteps))
        ALLOCATE(psnapshot_info%nparticles(psnapshot_info%ntimesteps))

	END SUBROUTINE init_psnapshots

	!------------------------------

	SUBROUTINE init_psnapshot_subfolder()

		CHARACTER(len = mglet_filename_max) :: subfolder

		psnapshot_info%counter = psnapshot_info%counter + 1 ! should be broadcasted via MPI so that no missmatches occur ?

		IF (myid == 0) THEN

	    	WRITE(subfolder, '("PARTICLE_SNAPSHOTS/snapshot", I0)') psnapshot_info%counter
            CALL create_directory(TRIM(subfolder)) ! ! ! realtive to working directory ! ! !

        END IF

        CALL MPI_Barrier(MPI_COMM_WORLD)

	END SUBROUTINE init_psnapshot_subfolder

	!------------------------------

	SUBROUTINE write_psnapshot_piece(my_particle_list, time)
	! writes vtp (xml) file containing all particles of the respective process   

		! subroutine arguments 
		TYPE(particle_list_t) :: my_particle_list
		REAL(realk) :: time

		! local variables 
		INTEGER :: i, unit
		CHARACTER(len = mglet_filename_max) :: subfolder, filename
		CHARACTER(:), ALLOCATABLE :: active_np_char

	    WRITE(subfolder, '("PARTICLE_SNAPSHOTS/snapshot", I0)') psnapshot_info%counter

		WRITE(filename, '(A, "/piece", I0, ".vtp")') TRIM(subfolder), myid

		ALLOCATE(CHARACTER(CEILING(LOG10(REAL(my_particle_list%max_np))) :: active_np_char)
		WRITE(active_np_char, '(I0)') my_particle_list%active_np

		OPEN(unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')
	
		WRITE(unit, '(A)') '<?xml version="1.0"?>'
		WRITE(unit, '(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
		WRITE(unit, '(A)') '  <PolyData>'
		WRITE(unit, '(A)') '    <Piece NumberOfPoints="' // TRIM(active_np_char) // & 
						  	'" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
		WRITE(unit, '(A)') '      <Points>'
		WRITE(unit, '(A)') '        <DataArray type="Float32" NumberOfComponents="3">'

		DO i = 1, my_particle_list%ifinal

			IF (.NOT. my_particle_list%particle_stored(i)) THEN
				CYCLE
			END IF

			WRITE(unit, psnapshot_info%coordinate_format, advance="no") my_particle_list%particles(i)%x
			WRITE(unit, '(A)', advance="no") ' '
			WRITE(unit, psnapshot_info%coordinate_format, advance="no") my_particle_list%particles(i)%y
			WRITE(unit, '(A)', advance="no") ' '
			WRITE(unit, psnapshot_info%coordinate_format, advance="yes") my_particle_list%particles(i)%z

		END DO 

		! also WRITE point data (e.g. velocity) and vertices for visualization? 

		WRITE(unit, '(A)') '        </DataArray>'
		WRITE(unit, '(A)') '      </Points>'
		WRITE(unit, '(A)') '    </Piece>'
		WRITE(unit, '(A)') '  </PolyData>'
		WRITE(unit, '(A)') '</VTKFile>'

	END SUBROUTINE write_psnapshot_piece

	!------------------------------

	SUBROUTINE write_psnapshot_master()
	! writes pvtp (xml) file (master file for all pieces of the same time)
	
		INTEGER :: proc, unit
		CHARACTER(len = mglet_filename_max) :: filename, piece

		IF (myid == 0) THEN

			WRITE(filename,'("PARTICLE_SNAPSHOTS/snapshot", I0, ".pvtp")') , psnapshot_info%counter

			OPEN(unit, file = TRIM(filename), status = 'NEW', action = 'WRITE')

			WRITE(unit, '(A)') '<?xml version="1.0"?>'
			WRITE(unit, '(A)') '<VTKFile type="PPolyData" version="0.1" byte_order="LittleEndian">'
			WRITE(unit, '(A)') '  <PPolyData>'
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

	SUBROUTINE write_psnapshot_info()

	END SUBROUTINE write_psnapshot_info

	!------------------------------

	! TODO: particle trajectories

END MODULE 
