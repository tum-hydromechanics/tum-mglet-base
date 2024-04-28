! this file is for particle lists

MODULE particles_mod

!===================================

USE precision_mod, ONLY: intk, realk
USE baseparticle_mod 

!===================================

IMPLICIT NONE

TYPE :: particle_list_t

	INTEGER(intk) :: list_id = myid 	! not sure if useful!?
	INTEGER(intk) :: n_active_particles
	INTEGER(intk) :: list_length
	INTEGER(intk) :: n_final 			! index of last entry of the list which holds an active particle


	TYPE(base_particle_t), ALLOCATABLE :: particles(:)
	LOGICAL, ALLOCATABLE :: particle_stored(:) 	! each logical value reflects whether a particle is stored in the list 
												! at the respective index. Is this a feasable and good way to keep track 
												! of particle storage (especially as is_init in particle_t carries the same information?
												

END TYPE particle_list_t

! local variables

INTEGER(intk), PARAMETER :: default_list_length = 1000
INTEGER(intk), PARAMETER :: initial_n_particles_per_list = 100 ! RENAME ?
TYPE(particle_list_t) :: my_particle_list

CONTAINS

	SUBROUTINE init_particles()

		! For the beginning, every process creates his own actual list of particles.
		! Would it be better to have one global list of particles and give each process 
		! a list of pointers to elemnts of the global list, or does this lead to difficulties 
		! regarding memory access?


		! local variables
 		INTEGER(intk) :: i
 		REAL(realk) :: time = 0 
 		REAL(realk) :: x, y, z


 		my_particle_list%n_active_particles = initial_n_particles_per_list
		my_particle_list%list_length = default_list_length
		my_particle_list%n_final = initial_n_particles_per_list

		ALLOCATE(my_particle_list%particles(my_particle_list%list_length))
		ALLOCATE(my_particle_list%particle_stored(my_particle_list%list_length))

		my_particle_list%particle_stored = .FALSE.

 		DO i = 1, my_particle_list%n_active_particles

 			CALL random_ic(x, y, z)

 			CALL my_particle_list%particles(i)%init(ipart = xxx, iproc = my_particle_list%list_id, igrid = xxx, &
 			nearest_cell = xxx, x, y, z, time)

 			my_particle_list%particle_stored(i) = .TRUE.

 		END DO
	

	END SUBROUTINE init_particles

	SUBROUTINE dist_ipart()

	! this routine is supposed to hand out a list of unique particle ids (ipart) to every process

	END SUBROUTINE dist_ipart

END MODULE 
