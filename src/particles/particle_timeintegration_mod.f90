! file for particle motion and timeintegration

MODULE particle_timeintegration_mod

!===================================

USE precision_mod, ONLY: intk, realk
USE core_mod  
USE flowcore_mod
USE baseparticle_mod 
USE particle_list_mod

!===================================

REAL(realk), PARAMETER :: D = 0.1 ! Diffusion constant in m²/s (homogeneous and isotropic diffusion for now)

!===================================

CONTAINS

	SUBROUTINE timeintegrate_particles(dt)

		! subroutine arguments 
		REAL(realk) :: dt

		! local variables		
		INTEGER(intk) :: igrid, p
		REAL(realk) :: p_u, p_v, p_w

		TYPE(field_t), POINTER :: x_f, y_f, z_f
		TYPE(field_t), POINTER :: xstag_f, ystag_f, zstag_f
		TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
		TYPE(field_t), POINTER :: u_f, v_f, w_f

		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: x, y, z
		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: xstag, ystag, zstag
		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ddx, ddy, ddz
		REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w

		REAL(realk) :: rand_dx, rand_dy, rand_dz, diffusion_dx, diffusion_dy, diffusion_dz

		!is calling the geometric fields obsolete when using core_mode -> corefields_mod ?
		CALL get_field(x_f, "X")
        CALL get_field(y_f, "Y")
        CALL get_field(z_f, "Z")

        CALL get_field(xstag_f, "XSTAG")
        CALL get_field(ystag_f, "YSTAG")
        CALL get_field(zstag_f, "ZSTAG")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")

    	DO p = 1, my_particle_list%ifinal

    		IF (.NOT. my_particle_list%particles(p)%is_init) THEN

    			CYCLE

    		END IF 

    		igrid = my_particle_list%particles(p)%igrid !just for clarity of the coming expressions

        	CALL x_f%get_ptr(x, igrid)
       		CALL y_f%get_ptr(y, igrid)
        	CALL z_f%get_ptr(z, igrid)

        	CALL ddx_f%get_ptr(ddx, igrid)
        	CALL ddy_f%get_ptr(ddy, igrid)
        	CALL ddz_f%get_ptr(ddz, igrid)

        	CALL xstag_f%get_ptr(xstag, igrid)
        	CALL ystag_f%get_ptr(ystag, igrid)
        	CALL zstag_f%get_ptr(zstag, igrid)       	

        	CALL u_f%get_ptr(u, igrid)
        	CALL v_f%get_ptr(v, igrid)
        	CALL w_f%get_ptr(w, igrid)

    		CALL interp_particle_uvw(my_particle_list%particles(p),&
    			p_u, p_v, p_w, xstag, ystag, zstag, ddx, ddy, ddz, u, v, w)

    		! TODO: proper timeintegration
    		! dummy method for now (explicit euler)

    		! how should 2 dimensional cases be handeled ? 
    		CALL RANDOM_NUMBER(rand_dx)
    		CALL RANDOM_NUMBER(rand_dy)
    		CALL RANDOM_NUMBER(rand_dz)

    		diffusion_dx = SQRT(2 * D * dt) * rand_dx / SQRT(rand_dx**(2) + rand_dy**(2) + rand_dz**(2))
    		diffusion_dy = SQRT(2 * D * dt) * rand_dx / SQRT(rand_dx**(2) + rand_dy**(2) + rand_dz**(2))
    		diffusion_dz = SQRT(2 * D * dt) * rand_dx / SQRT(rand_dx**(2) + rand_dy**(2) + rand_dz**(2))

    		particles(p)%x = particles(p)%x + p_u * dt + diffusion_dx
    		particles(p)%y = particles(p)%y + p_v * dt + diffusion_dy
    		particles(p)%z = particles(p)%z + p_w * dt + diffusion_dz

    		! TODO: exit domain handling / deactivate particles that leave the domain

        END DO

	END SUBROUTINE timeintegrate_particles 

	!------------------------------

	SUBROUTINE interp_particle_uvw(particle, p_u, p_v, p_w, xstag, ystag, zstag, ddx, ddy, ddz, u, v, w)

			! subroutine arguments
			TYPE(base_particle_t), INTENT(in) :: particle
			REAL(realk), INTENT(out), p_u, p_v, p_w
			REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :), INTENT(in) :: xstag, ystag, zstag
			REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :), INTENT(in) :: ddx, ddy, ddz
			REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :), INTENT(in) :: u, v, w

			! local variables
			INTEGER(intk) :: p_i, p_j, p_k, p_istag, p_jstag, p_kstag
			REAL(realk) :: p_x, p_y, p_z
			REAL(realk), DIMENSION(6) :: p_bbox_u, p_bbox_v, p_bbox_w


	    	p_i = particle%ijkcell(1) !just for clarity of the coming expressions
    		p_j = particle%ijkcell(2) !just for clarity of the coming expressions
    		p_k = particle%ijkcell(3) !just for clarity of the coming expressions

    		p_x = particle%x !just for clarity of the coming expressions
    		p_y = particle%y !just for clarity of the coming expressions
    		p_z = particle%z !just for clarity of the coming expressions

    		p_istag = CEILING((p_x - x(p_i)) / ddx(p_i))
    		p_jstag = CEILING((p_y - y(p_j)) / ddy(p_j))
    		p_kstag = CEILING((p_z - z(p_k)) / ddz(p_k))

		    p_bbox_u(1) = xstag(p_i - 1)
    		p_bbox_u(2) = xstag(p_i)
    		p_bbox_u(3) = ystag(p_j + p_jstag - 1)
    		p_bbox_u(4) = ystag(p_j + p_jstag)
    		p_bbox_u(5) = zstag(p_k + p_kstag - 1)
    		p_bbox_u(6) = zstag(p_k + p_kstag) 

    		p_bbox_v(1) = xstag(p_i + p_istag - 1)
    		p_bbox_v(2) = xstag(p_i + p_istag)
    		p_bbox_v(3) = ystag(p_j - 1)
    		p_bbox_v(4) = ystag(p_j)
    		p_bbox_v(5) = zstag(p_k + p_kstag - 1)
    		p_bbox_v(6) = zstag(p_k + p_kstag)

    		p_bbox_w(1) = xstag(p_i + p_istag - 1)
    		p_bbox_w(2) = xstag(p_i + p_istag)
    		p_bbox_w(3) = ystag(p_j + p_jstag - 1)
    		p_bbox_w(4) = ystag(p_j + p_jstag)
    		p_bbox_w(5) = zstag(p_k - 1)
    		p_bbox_w(6) = zstag(p_k)

    		interp_trilinear(p_x, p_y, p_z, p_bbox_u, p_u,&
    						u(p_k + p_kstag - 1, p_j + p_jstag -1, p_i - 1),&
    						u(p_k + p_kstag    , p_j + p_jstag -1, p_i - 1),&
    						u(p_k + p_kstag - 1, p_j + p_jstag   , p_i - 1),&
    						u(p_k + p_kstag    , p_j + p_jstag   , p_i - 1),&
    						u(p_k + p_kstag - 1, p_j + p_jstag -1, p_i    ),&
    						u(p_k + p_kstag    , p_j + p_jstag -1, p_i    ),&
    						u(p_k + p_kstag - 1, p_j + p_jstag   , p_i    ),&
    						u(p_k + p_kstag    , p_j + p_jstag   , p_i    ))

     		interp_trilinear(p_x, p_y, p_z, p_bbox_v, p_v,&
    						v(p_k + p_kstag - 1, p_j -1, p_i + p_istag - 1),&
    						v(p_k + p_kstag    , p_j -1, p_i + p_istag - 1),&
    						v(p_k + p_kstag - 1, p_j   , p_i + p_istag - 1),&
    						v(p_k + p_kstag    , p_j   , p_i + p_istag - 1),&
    						v(p_k + p_kstag - 1, p_j -1, p_i + p_istag    ),&
    						v(p_k + p_kstag    , p_j -1, p_i + p_istag    ),&
    						v(p_k + p_kstag - 1, p_j   , p_i + p_istag    ),&
    						v(p_k + p_kstag    , p_j   , p_i + p_istag    ))

    		interp_trilinear(p_x, p_y, p_z, p_bbox_w, p_w,&
    						w(p_k - 1, p_j + p_jstag -1, p_i + p_istag - 1),&
    						w(p_k    , p_j + p_jstag -1, p_i + p_istag - 1),&
    						w(p_k - 1, p_j + p_jstag   , p_i + p_istag - 1),&
    						w(p_k    , p_j + p_jstag   , p_i + p_istag - 1),&
    						w(p_k - 1, p_j + p_jstag -1, p_i + p_istag    ),&
    						w(p_k    , p_j + p_jstag -1, p_i + p_istag    ),&
    						w(p_k - 1, p_j + p_jstag   , p_i + p_istag    ),&
    						w(p_k    , p_j + p_jstag   , p_i + p_istag    ))

	SUBROUTINE interp_particle_uvw

	!------------------------------

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
