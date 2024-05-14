! file for particle motion and timeintegration

MODULE particle_timeintegration_mod

USE precision_mod, ONLY: intk, realk
USE baseparticle_mod 
USE particles_mod
USE core_mod  
USE flowcore_mod


CONTAINS

	SUBROUTINE timeintegrate_particles()

		! local variables

		TYPE(field_t), POINTER :: x_f, y_f, z_f
		TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
		TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
		TYPE(field_t), POINTER :: u_f, v_f, w_f
		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: x, y, z
		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: dx, dy, dz
		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ddx, ddy, ddz
		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
		
		
		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: particles ! use contiguous here?
		REAL(realk), DIMENSION(6) :: p_bbox(6)
		REAL(realk) :: p_u, p_v, p_w, 
		INTEGER(intk) :: p, p_istag, p_jstag, p_kstag

		!is calling the geometric fields obsolete when using core_mode -> corefields_mod ?
		CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(xstag_f, "XSTAG")
        CALL get_field(ystag_f, "YSTAG")
        CALL get_field(zstag_f, "ZSTAG")

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

        CALL x_f%get_ptr(x, this%igrid)
        CALL y_f%get_ptr(y, this%igrid)
        CALL z_f%get_ptr(z, this%igrid)

        CALL xstag_f%get_ptr(xstag, this%igrid)
        CALL ystag_f%get_ptr(ystag, this%igrid)
        CALL zstag_f%get_ptr(zstag, this%igrid)

        CALL dx_f%get_ptr(dx, this%igrid)
        CALL dy_f%get_ptr(dy, this%igrid)
        CALL dz_f%get_ptr(dz, this%igrid)

        CALL ddx_f%get_ptr(ddx, this%igrid)
        CALL ddy_f%get_ptr(ddy, this%igrid)
        CALL ddz_f%get_ptr(ddz, this%igrid)

        CALL u_f%get_ptr(u, this%igrid)
        CALL v_f%get_ptr(v, this%igrid)
        CALL w_f%get_ptr(w, this%igrid)

        particles => my_particle_list%particles

    	DO p = 1, my_particle_list%ifinal

    		IF (.NOT. particles(p)%is_init) THEN

    			CYCLE

    		END IF 

    		p_istag = CEILING((particles(p)%x - x(ijkcell(1))) / ddx(ijkcell(1)))
    		p_jstag = CEILING((particles(p)%y - x(ijkcell(2))) / ddy(ijkcell(2)))
    		p_kstag = CEILING((particles(p)%z - x(ijkcell(3))) / ddz(ijkcell(3)))

    		p_bbox = /( x(ijkcell(1) + p_istag - 1), x(ijkcell(1) + p_istag), &
    					x(ijkcell(2) + p_jstag - 1), x(ijkcell(2) + p_jstag), &
    					x(ijkcell(3) + p_kstag - 1), x(ijkcell(3) + p_kstag) /) 

    		interp_trilinear( particles(p)%x, particles(p)%y, particles(p)%z, p_bbox, p_u,&
    						u(ijkcell(3) + p_kstag - 1, ijkcell(2) + p_jstag -1, ijkcell(1) - 1),&
    						u(ijkcell(3) + p_kstag    , ijkcell(2) + p_jstag -1, ijkcell(1) - 1),&
    						u(ijkcell(3) + p_kstag - 1, ijkcell(2) + p_jstag   , ijkcell(1) - 1),&
    						u(ijkcell(3) + p_kstag    , ijkcell(2) + p_jstag   , ijkcell(1) - 1),&
    						u(ijkcell(3) + p_kstag - 1, ijkcell(2) + p_jstag -1, ijkcell(1)    ),&
    						u(ijkcell(3) + p_kstag    , ijkcell(2) + p_jstag -1, ijkcell(1)    ),&
    						u(ijkcell(3) + p_kstag - 1, ijkcell(2) + p_jstag   , ijkcell(1)    ),&
    						u(ijkcell(3) + p_kstag    , ijkcell(2) + p_jstag   , ijkcell(1)    ) )

     		interp_trilinear( particles(p)%x, particles(p)%y, particles(p)%z, p_bbox, p_v,&
    						v(ijkcell(3) + p_kstag - 1, ijkcell(2) -1, ijkcell(1) + p_istag - 1),&
    						v(ijkcell(3) + p_kstag    , ijkcell(2) -1, ijkcell(1) + p_istag - 1),&
    						v(ijkcell(3) + p_kstag - 1, ijkcell(2)   , ijkcell(1) + p_istag - 1),&
    						v(ijkcell(3) + p_kstag    , ijkcell(2)   , ijkcell(1) + p_istag - 1),&
    						v(ijkcell(3) + p_kstag - 1, ijkcell(2) -1, ijkcell(1) + p_istag    ),&
    						v(ijkcell(3) + p_kstag    , ijkcell(2) -1, ijkcell(1) + p_istag    ),&
    						v(ijkcell(3) + p_kstag - 1, ijkcell(2)   , ijkcell(1) + p_istag    ),&
    						v(ijkcell(3) + p_kstag    , ijkcell(2)   , ijkcell(1) + p_istag    ) )

    		interp_trilinear( particles(p)%x, particles(p)%y, particles(p)%z, p_bbox, p_w,&
    						w(ijkcell(3) - 1, ijkcell(2) + p_jstag -1, ijkcell(1) + p_istag - 1),&
    						w(ijkcell(3)    , ijkcell(2) + p_jstag -1, ijkcell(1) + p_istag - 1),&
    						w(ijkcell(3) - 1, ijkcell(2) + p_jstag   , ijkcell(1) + p_istag - 1),&
    						w(ijkcell(3)    , ijkcell(2) + p_jstag   , ijkcell(1) + p_istag - 1),&
    						w(ijkcell(3) - 1, ijkcell(2) + p_jstag -1, ijkcell(1) + p_istag    ),&
    						w(ijkcell(3)    , ijkcell(2) + p_jstag -1, ijkcell(1) + p_istag    ),&
    						w(ijkcell(3) - 1, ijkcell(2) + p_jstag   , ijkcell(1) + p_istag    ),&
    						w(ijkcell(3)    , ijkcell(2) + p_jstag   , ijkcell(1) + p_istag    ) )

    		! TODO: proper timeintegration
    		! dummy method for now (explicit euler)
    		particles(p)%x = particles(p)%x + p_u * dt
    		particles(p)%y = particles(p)%y + p_v * dt
    		particles(p)%z = particles(p)%z + p_w * dt
    		
    		! TODO: use timekeeper for particle time ! ! !

    		! TODO: new ijkcell

        END DO

	END SUBROUTINE timeintegrate_particles 

	SUBROUTINE interp_trilinear(x, y, z, bbox &
		f, f19, f20, f21, f22, f23, f24, f25, f26) ! Trilinear interpolation for a rectilinear grid !

		! subroutine arguments 
		REAL(realk), INTENT(in) :: x, y, z
		REAL(realk), DIMENSION(6), INTENT(in) :: bbox
		REAL(realk), INTENT(in) :: f19, f20, f21, f22, f23, f24, f25, f26 ! number indicate the corners (see convention on the numbering of faces/edges/corners)
		REAL(realk), INTENT(out) :: f

		! local variables
		REAL(realk) :: xd, yd, zd

		xd = (x - bbox(1)) / (bbox(2) - bbox(1))
		yd = (y - bbox(3)) / (bbox(4) - bbox(3))
		zd = (z - bbox(5)) / (bbox(6) - bbox(5))

		f = (f19 * (1 - xd) + f23 * xd)  * (1 - yd) + (f21 * (1 - xd) + f25 * xd) * yd * (1 - zd) + &
			(f20 * (1 - xd) + f24 * xd) * (1 - yd) + (f22 * (1 - xd) + f26 * xd) * yd * zd

	END SUBROUTINE interpolate_trilinear

END MODULE 
