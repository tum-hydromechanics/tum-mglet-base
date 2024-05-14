! this file is for particle lists

MODULE particle_list_mod

!===================================

USE precision_mod, ONLY: intk, realk
USE baseparticle_mod 

!===================================

IMPLICIT NONE

TYPE :: particle_list_t

	INTEGER(intk) :: iproc = myid	

	INTEGER(intk) :: max_np    			! max number of particles of this process/list
	INTEGER(intk) :: active_np 			! number of active particles of this process/list
	INTEGER(intk) :: ifinal 			! index of last entry of the list which holds an active particle

	TYPE(base_particle_t), ALLOCATABLE :: particles(:)
	LOGICAL, ALLOCATABLE :: particle_stored(:) 	! each logical value reflects whether a particle is stored in the list 
												! at the respective index. Is this a feasable and good way to keep track 
												! of particle storage (especially as is_init in particle_t carries the same information?
												

END TYPE particle_list_t

! local variables

INTEGER(intk), PARAMETER :: default_max_np = 1000 
INTEGER(intk), PARAMETER :: default_initial_np = 100 !ONLY DUMMY VALUE FOR NOW, SHOULD BE SCALES WITH THE SIZE OF THE SPATIAL DOMAIN THAT THE PROCESS HANDLES
TYPE(particle_list_t) :: my_particle_list

INTEGER(intk), ALLOCATABLE :: my_ipart_arr(:)

CONTAINS

	SUBROUTINE init_particles()

		! local variables
 		INTEGER(intk) :: i
 		REAL(realk) :: time = 0 
 		REAL(realk) :: x, y, z

		my_particle_list%max_np = default_max_np
 		my_particle_list%active_np = default_initial_np
		my_particle_list%ifinal = default_initial_np

		ALLOCATE(my_particle_list%particles(my_particle_list%max_np))
		ALLOCATE(my_particle_list%particle_stored(my_particle_list%max_np))

		my_particle_list%particle_stored = .FALSE.

 		DO i = 1, my_particle_list%active_np

 			CALL random_ic(x, y, z) ! ONLY DUMMY FOR NOW, DENPENDS ON THE PROCESS SPATIAL DOMAIN

 			CALL my_particle_list%particles(i)%init(ipart = my_ipart_arr(i), x = x, y = y, z = z, time = time)

 			my_particle_list%particle_stored(i) = .TRUE.

 		END DO
	

	END SUBROUTINE init_particles

	SUBROUTINE dist_ipart(ipart_arr) ! this routine is supposed to hand out a list of unique particle ids (ipart) to every process ! ONLY DUMMY FOR NOW

	! subroutine arguments
	INTEGER, ALLOCATABLE, INTENT(out) :: ipart_arr(:)
	
	! local variables
	INTEGER :: i

 	ALLOCATE(ipart_arr(default_initial_np)) 
	
	ipart_arr = (/(i, i = myid * default_initial_np + 1, myid * default_initial_np + default_initial_np)/)

	END SUBROUTINE dist_ipart

END MODULE 
