! a module for particle input/output operations, if initial condition (ic) 
! or boundary conditions (bc) should be read from a file or information should 
! be written to a file 

MODULE particle_io_mod

! TODO: psnapshot_info_t to store number of serial files per time, times and accuracy/format of point coordinates
! (sort of as a log file) 

USE

!===================================

! VARIABLES

!===================================

CONTAINS

	SUBROUTINE read_particle_ic()

	END SUBROUTINE read_particle_ic

	!------------------------------

	SUBROUTINE read_particle_bc()

	END SUBROUTINE read_particle_bc

	!------------------------------

	SUBROUTINE write_psnapshot_serial(psnapshot_info, my_particle_list, time)
	! writes vtp (xml) file containing all particles of the respective process   

		! subroutine arguments 
		TYPE(particle_list_t) :: my_particle_list
		TYPE(psnapshot_info_t) :: psnapshot_info

		! local arguments 
		INTEGER :: i, myid_char_len, unit
		CARACTER(128) :: filepath 
		CARACTER(:), ALLOCATABLE :: myid_char, time_char, format, ifinal_char

		myid_char_len = FLOOR(log10(nprocs)) + 1
		ALLOCATE(myid_char(myid_char_len))

		CALL itoa(myid, myid_char)
		CALL itoa(time, time_char)
		CALL itoa(my_particle_list%ifinal, ifinal_char)

		filepath = './particle-vtk/' // 'snapshot_' // TRIM(time_char) // ' ' // TRIM(myid_char) // '.vtp'

		open(unit, file = filepath, status = 'new', action = 'write')
	
		write(unit,'(a)') '<?xml version="1.0"?>'
		write(unit,'(a)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
		write(unit,'(a)') '  <PolyData>'
		write(unit,'(a)') '    <Piece NumberOfPoints="' // TRIM(ifinal_char) // '" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
		write(unit,'(a)') '      <Points>'
		write(unit,'(a)') '        <DataArray type="Float32" NumberOfComponents="3">'

		DO i = 1, my_particle_list%ifinal

			IF (.NOT. my_particle_list%particle_stored(i)) THEN
				CYCLE
			END IF

			write(unit, format, advance="no") my_particle_list%particles(i)%x
			write(unit, '(a)', advance="no") ' '
			write(unit, format, advance="no") my_particle_list%particles(i)%y
			write(unit, '(a)', advance="no") ' '
			write(unit, format, advance="yes") my_particle_list%particles(i)%z

		END DO 

		! also write vertices? 

		write(unit,'(a)') '        </DataArray>'
		write(unit,'(a)') '      </Points>'
		write(unit,'(a)') '    </Piece>'
		write(unit,'(a)') '  </PolyData>'
		write(unit,'(a)') '</VTKFile>'

	END SUBROUTINE write_psnapshot_serial

	!------------------------------

	SUBROUTINE write_psnapshot_parallel(psnapshot_info)
	! writes pvtp (xml) file (master file for all snapshots)
	


	END SUBROUTINE write_psnapshot_parallel

	!------------------------------

	SUBROUTINE init_particle_trajectories()

	END SUBROUTINE init_particle_trajectories

	!------------------------------

	SUBROUTINE add_particle_trajectories()

	END SUBROUTINE add_particle_trajectories



END MODULE 
