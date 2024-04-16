! this file is for particle lists

MODULE particles_mod

!===================================

USE precision_mod, ONLY: intk, realk
USE baseparticle_mod 

!===================================

IMPLICIT NONE

INTEGER(intk), PARAMETER :: default_list_length = 1000
INTEGER(intk), PARAMETER :: default_n_particle_per_list = 100

TYPE :: particle_list_t

	INTEGER(intk) :: list_id = myid 	! not sure if useful!?
	INTEGER(intk) :: n_active_particles
	INTEGER(intk) :: list_length


	TYPE(base_particle_t), ALLOCATABLE :: particles(:)
	LOGICAL, ALLOCATABLE :: particle_stored(:) 	! each logical value reflects whether a particle is stored in the list 
												! at the respective index. Is this a feasable and good way to keep track 
												! of particle storage (especially as is_init in particle_t carries the same information?

	! CONTAINS ? 

END TYPE particle_list_t

TYPE(particle_list_t) :: my_particle_list

CONTAINS

	SUBROUTINE init_particles()

		! For the beginning, every process creates his own actual list of particles.
		! Would it be better to have one global list of particles and give each process 
		! a list of pointers to elemnts of the global list, or does this lead to difficulties 
		! regarding memory access?


		! local variables
 		INTEGER(intk) :: local_id
 		REAL(realk) :: particle_time = 0 
 		REAL(realk) :: x, y, z
		REAL(realk) :: u, v, w


 		my_particle_list%n_active_particles = default_n_particle_per_list
		my_particle_list%list_length = default_list_length

		ALLOCATE(my_particle_list%particles(my_particle_list%list_length))
		ALLOCATE(my_particle_list%particle_stored(my_particle_list%list_length))

		my_particle_list%particle_stored = .FALSE.

 		DO local_id = 1, my_particle_list%n_active_particles

 			CALL random_ic(x, y, z, u, v, w)

 			CALL my_particle_list%particles(local_id)%init(local_id, my_particle_list%list_id, x, y, z, u, v, w, particle_time)

 		END DO
	

	END SUBROUTINE init_particles


END MODULE 
