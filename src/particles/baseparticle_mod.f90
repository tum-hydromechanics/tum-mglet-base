! this file is for TYPE definitions and general handling of particles 

MODULE baseparticle_mod

!===================================

USE precision_mod, ONLY: intk, realk
USE core_mod ! should be specified in more detail later
USE flowcore_mod ! should be specified in more detail later

!===================================

IMPLICIT NONE 

TYPE :: base_particle_t ! could be extended by a particle type that includes the mass
	
	LOGICAL :: is_init = .FALSE. ! if a list of particles is handled, this can be used to identify if an entry is an actual particly

	INTEGER(intk) :: ipart  ! particle ID in global scope (of all processes)
	INTEGER(intk) :: iproc	! id of process that currently handles the particle
	INTEGER(intk) :: igrid 	! stores the index of grid, which the particle is currently on 

	INTEGER(intk) :: ijk_cell(3) ! stores indices to the nearest grid cell 

	REAL(realk) :: x, y, z
	REAL(realk) :: time

	REAL(realk) :: diffusion_len ! or other diffusion specific parameter 

	CONTAINS 

	PROCEDURE :: init
	PROCEDURE :: print_status
	PROCEDURE :: get_p_igrid
	PROCEDURE :: get_p_cell

	!GENERIC, PUBLIC :: ijk_to_ptr => ijk_to_ptr1, ijk_to_ptr3 ! NOT NEEDED
	!PROCEDURE, PRIVATE :: ijk_to_ptr1, ijk_to_ptr3	! NOT NEEDED

	! PROCEDURE :: get_ic  (could be another method) 

END TYPE base_particle_t

!===================================

CONTAINS

	SUBROUTINE init(this, ipart, iproc, igrid, ijk_cell, x, y, z, time)

		! Subroutine arguments
	    CLASS(base_particle_t), INTENT(out) :: this
	    INTEGER(intk), INTENT(in) :: ipart
	    INTEGER(intk), INTENT(in) :: iproc
	    INTEGER(intk), INTENT(in) :: igrid
	    INTEGER(intk), INTENT(in) :: ijk_cell(3)
	    REAL(realk), INTENT(in) :: x, y, z

		REAL(realk), INTENT(in) :: time

		! Local variables
        	! none...

		this%is_init = .TRUE.

		this%ipart = ipart
		this%iproc = iproc
		this%igrid = !TODO

		this%ijk_cell = !TODO

		this%x = x
		this%y = y
		this%z = z

		this%time = time

	END SUBROUTINE init
	
	!------------------------------
	
	SUBROUTINE get_p_igrid(this)

		! Subroutine arguments
		CLASS(base_particle_t), INTENT(inout) :: this

		!local variables:
		LOGICAL :: found = .FALSE.
		INTEGER(intk) :: i, igrid
		REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        DO i = 1, nmygrids
            igrid = mygrids(i)

			minx = gridinfo(igrid)%bbox(1)
			maxx = gridinfo(igrid)%bbox(2)
			miny = gridinfo(igrid)%bbox(3)
			maxy = gridinfo(igrid)%bbox(4)
			minz = gridinfo(igrid)%bbox(5)
			maxz = gridinfo(igrid)%bbox(6)

			IF (this%x < minx) THEN
				CYCLE
			END IF

			IF (this%x > maxx) THEN
				CYCLE
			END IF

			IF (this%y < miny) THEN
				CYCLE
			END IF

			IF (this%y > maxy) THEN
				CYCLE
			END IF

			IF (this%z < minz) THEN
				CYCLE
			END IF

			IF (this%z > maxz) THEN
				CYCLE
			END IF

			found = .TRUE.
			this%igrid = igrid

		END DO 

		this%is_init = found ! if particle is not found, deactivate particle ???
	
	END SUBROUTINE get_p_igrid

	!------------------------------
	
	SUBROUTINE get_p_cell(this)

		! subroutine arguments
		CLASS(base_particle_t), INTENT(inout) :: this

		! local variables
		TYPE(field_t) :: x_f, y_f, z_f
        REAL(realk), POINTER, CONTIGUOUS :: x, y, z
        INTEGER(intk) :: k, j, i, kk, jj, ii
        REAL(realk) :: diff_old, diff_new

		CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL x_f%get_ptr(x, this%igrid)
        CALL y_f%get_ptr(y, this%igrid)
        CALL z_f%get_ptr(w, this%igrid)

        CALL get_mgdims(kk, jj, ii, this%igrid)

        ! the following assumes that the x/y/z values are sorted ! ! ! CHECK IF THIS IS TRUE

        diff_old = x(1) - this%x

        DO i = 2, ii
        	diff_new = x(i) - this%x

        	IF (diff_new > diff_old) THEN
        		this%ijk_cell(1) = i 
        		EXIT
        	ELSEIF (i == ii) THEN
        		this%ijk_cell(1) = i 
        	ELSE 
        		diff_old = diff_new
        		CYCLE
        	END IF
        END DO

        diff_old = y(1) - this%y

        DO j = 2, jj
        	diff_new = y(i) - this%y

        	IF (diff_new > diff_old) THEN
        		this%ijk_cell(2) = j
        		EXIT
        	ELSEIF (i == ii) THEN
        		this%ijk_cell(2) = j
        	ELSE 
        		diff_old = diff_new
        		CYCLE
        	END IF
        END DO

        diff_old = z(1) - this%z

        DO k = 2, kk
        	diff_new = z(i) - this%z

        	IF (diff_new > diff_old) THEN
        		this%ijk_cell(3) = k 
        		EXIT
        	ELSEIF (i == ii) THEN
        		this%ijk_cell(3) = k
        	ELSE 
        		diff_old = diff_new
        		CYCLE
        	END IF
        END DO
	
	END SUBROUTINE get_p_cell

	!------------------------------

	SUBROUTINE use_ijk_cell

		CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL x_f%get_ptr(x, this%igrid)
        CALL y_f%get_ptr(y, this%igrid)
        CALL z_f%get_ptr(w, this%igrid)

        nearest_cell_x_ptr = x(ijk_cell(1))
        nearest_cell_y_ptr = x(ijk_cell(2))
        nearest_cell_z_ptr = x(ijk_cell(3))

        !dxi/ddxi analogously ...

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        CALL x_f%get_ptr(u, this%igrid)
        CALL y_f%get_ptr(v, this%igrid)
        CALL z_f%get_ptr(w, this%igrid)

        nearest_cell_u_ptr = u(ijk_cell(1))
        nearest_cell_v_ptr = v(ijk_cell(2))
        nearest_cell_w_ptr = w(ijk_cell(3))

	END SUBROUTINE use_ijk_cell

	!------------------------------

	SUBROUTINE print_status(this)
	
		CLASS(base_particle_t), INTENT(in) :: this
	
		IF (this%is_init) THEN
		print *, 'Particle ID: ', this%ipart
		END IF
		
	END SUBROUTINE print_status
	
	!-----------------------------
	
	SUBROUTINE random_ic(x, y, z)

	! should this be a method of baseparticle_t?

	! velocity of the massless particle should be equal to the veocity of the field at (x,y,z, t= time)

	! x, y, and z should be within the spacial domain covvered by the working process

	! for now only random values between 0 and 1:

		REAL(realk), INTENT(out) :: x, y, z

		CALL RANDOM_NUMBER(x)
		CALL RANDOM_NUMBER(y)
		CALL RANDOM_NUMBER(z)

	END SUBROUTINE random_ic



END MODULE 
