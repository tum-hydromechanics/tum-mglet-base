! this file is for TYPE definitions and general handling of particles 

MODULE particles_mod

!===================================

USE precision_mod, ONLY: intk, realk

!===================================

IMPLICIT NONE 

TYPE :: base_particle_t ! could be extended by a particle type that includes the mass
	
	LOGICAL :: is_init = .FALSE. ! if a list of particles is handled, 
									! this can be used to identify iif an entry is an actual particly
	INTEGER(intk) :: local_id
	INTEGER(intk) :: list_id 			! not sure if useful!?
	REAL(realk) :: x, y, z
	REAL(realk) :: u, v, w
	REAL(realk) :: particle_time

	!  :: local ?

	CONTAINS 

	PROCEDURE :: create
	PROCEDURE :: print_status
	! PROCEDURE :: get_ic  (could be another method) 

END TYPE base_particle_t

!===================================

CONTAINS

	SUBROUTINE init(this, local_id, list_id, x, y, z, u, v, w, particle_time)

			! Subroutine arguments
	    CLASS(base_particle_t), INTENT(out) :: this
	    INTEGER(intk), INTENT(in) :: local_id
	    INTEGER(intk), INTENT(in) :: list_id
	    REAL(realk), INTENT(in) :: x, y, z
		REAL(realk), INTENT(in) :: u, v, w
		REAL(realk) :: particle_time

		! Local variables
        ! none...

		this%is_created = .TRUE.

		this%local_id = local_id
		this%list_id = list_id

		this%x = x
		this%y = y
		this%z = z

		this%u = u
		this%v = v
		this%w = w

		this%particle_time = particle_time

	end SUBROUTINE

	SUBROUTINE random_ic(x, y, z, u, v, w)

	! should this be a method of baseparticle_t?

	! velocity of the massless particle should be equal to the veocity of the field at (x,y,z, t= particle_time)

	! x, y, and z should be within the spacial domain covvered by the working process

	! for now only random values between 0 and 1:

		REAL(realk), INTENT(out) :: x, y, z
		REAL(realk), INTENT(out) :: u, v, w

		CALL RANDOM_NUMBER(x)
		CALL RANDOM_NUMBER(y)
		CALL RANDOM_NUMBER(z)

		CALL RANDOM_NUMBER(u)
		CALL RANDOM_NUMBER(v)
		CALL RANDOM_NUMBER(w)

	END SUBROUTINE random_ic() 

	!------------------------------

	SUBROUTINE print_status()

	END SUBROUTINE print_status


END MODULE particles_mod
