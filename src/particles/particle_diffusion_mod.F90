MODULE particle_diffusion_mod

    USE grids_mod
    USE field_mod
    USE fields_mod
    USE connect2_mod
    USE fort7_mod

    USE particle_config_mod

    IMPLICIT NONE

    INTERFACE get_particle_diffusion
        MODULE PROCEDURE :: get_particle_diffusion_automatic
        MODULE PROCEDURE :: get_particle_diffusion_homogeneous
        MODULE PROCEDURE :: get_particle_diffusion_nearest
        MODULE PROCEDURE :: get_particle_diffusion_interpolated
    END INTERFACE get_particle_diffusion

CONTAINS

    SUBROUTINE init_particle_diffusion()

        ! local_variables
        INTEGER(intk), PARAMETER :: units_diff(7) = [0, 2, -1, 0, 0, 0, 0]

        IF (.NOT. turb_diff) THEN
            RETURN
        END IF

        ! fields that hold the total (trubulent + diffusive) particle diffusion
        ! if trub_diff == .TRUE. fields i/o will try to read the diffusion fields,
        ! if they cannot be read, then generate_turbulent_diffusion will try to generate the diffsuion fields
        ! if diffsuion fields cannot be read nor generated, the simulation is terminated
        CALL set_field("P_DIFF_X", istag = 1, units = units_diff, &
         dread = dcont, required = dcont, dwrite = .TRUE., buffers = .TRUE.)
        CALL set_field("P_DIFF_Y", jstag = 1, units = units_diff, &
         dread = dcont, required = dcont, dwrite = .TRUE., buffers = .TRUE.)
        CALL set_field("P_DIFF_Z", kstag = 1, units = units_diff, &
         dread = dcont, required = dcont, dwrite = .TRUE., buffers = .TRUE.)

        ! TODO: how to find out if a field was actually read if required = .FALSE.?
        IF (.NOT. dcont) THEN
            CALL generate_diffusion_field
        END IF

    END SUBROUTINE init_particle_diffusion

    SUBROUTINE generate_diffusion_field()

        ! local_variables
        INTEGER(intk) :: igrid, ii, jj, kk, g, i, j, k

        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz
        TYPE(field_t), POINTER :: t1_avg_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: t1_avg
        TYPE(field_t), POINTER :: u_avg_f, v_avg_f, w_avg_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u_avg, v_avg, w_avg
        TYPE(field_t), POINTER :: ut1_avg_f, vt1_avg, wt1_avg
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ut1_avg, vt1_avg, wt1_avg
        TYPE(field_t), POINTER :: diffx_f, diffy_f, diffz_d
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: diffx, diffy, diffz

        ! if one of the following fields does not exist, get_field will call an error and terminate the process
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(t1_avg_f, "T1_AVG")

        CALL get_field(u_avg_f, "U_AVG")
        CALL get_field(v_avg_f, "V_AVG")
        CALL get_field(w_avg_f, "W_AVG")

        CALL get_field(ut1_avg_f, "UT1_AVG")
        CALL get_field(vt1_avg_f, "VT1_AVG")
        CALL get_field(wt1_avg_f, "WT1_AVG")

        CALL get_field(diffx_f, "P_DIFF_X")
        CALL get_field(diffy_f, "P_DIFF_Y")
        CALL get_field(diffz_f, "P_DIFF_Z")

        DO g = 1, nmygrids

            igrid = mygrids(g)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL t1_avg_f%get_ptr(t1_avg, igrid)

            CALL u_avg_f%get_ptr(u_avg, igrid)
            CALL v_avg_f%get_ptr(v_avg, igrid)
            CALL w_avg_f%get_ptr(w_avg, igrid)

            CALL ut1_avg_f%get_ptr(ut1_avg, igrid)
            CALL vt1_avg_f%get_ptr(vt1_avg, igrid)
            CALL wt1_avg_f%get_ptr(wt1_avg, igrid)

            CALL diffx_f%get_ptr(diffx, igrid)
            CALL diffy_f%get_ptr(diffy, igrid)
            CALL diffz_f%get_ptr(diffz, igrid)

            DO i = 3, ii - 2
                DO j = 3, jj - 2
                    DO k = 3, kk - 2
                        diffx(k, j, i) = (ut1_avg(k, j, i) - u_avg(k, j, i) * 0.5 (t1_avg(k, j, i + 1) + t1_avg(k, j, i))) &
                         * dx(k, j, i) / (t1_avg(k, j, i + 1) - t1_avg(k, j, i)) + D(1)
                    END DO
                END DO
            END DO

            DO i = 3, ii - 2
                DO j = 3, jj - 2
                    DO k = 3, kk - 2
                        diffy(k, j, i) = (vt1_avg(k, j, i) - v_avg(k, j, i) * 0.5 (t1_avg(k, j + 1, i) + t1_avg(k, j, i))) &
                         * dy(k, j, i) / (t1_avg(k, j + 1, i) - t1_avg(k, j, i)) + D(2)
                    END DO
                END DO
            END DO

            DO i = 3, ii - 2
                DO j = 3, jj - 2
                    DO k = 3, kk - 2
                        diffz(k, j, i) = (wt1_avg(k, j, i) - w_avg(k, j, i) * 0.5 (t1_avg(k + 1, j, i) + t1_avg(k, j, i))) &
                         * dz(k, j, i) / (t1_avg(k + 1, j, i) - t1_avg(k, j, i)) + D(3)
                    END DO
                END DO
            END DO

            ! TODO: check this approach for buffers; esp. levels ...
            DO ilevel = minlevel, maxlevel
                CALL connect(ilevel, 1, v1 = diffx, v2 = diffy, v3 = diffz, corners = .TRUE.)
            END DO
        END DO

    END SUBROUTINE generate_diffusion_field

    SUBROUTINE get_particle_diffusion_automatic(particle, dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff, dinterp_diffsuion)

            ! subroutine arguments
            TYPE(baseparticle_t), INTENT(in) :: particle
            REAL(realk), INTENT(in) :: dt
            LOGICAL, INTENT(in) :: dinterp_diffsuion
            REAL(realk), INTENT(out) :: pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff

            ! local variables
            INTEGER(intk), INTENT(in) :: kk, jj, ii
            TYPE(field_t), POINTER :: x_f, y_f, z_f
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: x, y, z
            TYPE(field_t), POINTER :: dx_f, dy_f, dz_f
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx, dy, dz
            TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx, ddy, ddz
            TYPE(field_t), POINTER :: diffx_f, diffy_f, diffz_d
            REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: diffx, diffy, diffz

            ! local variables
            REAL(realk) :: D_x, D_y, D_z, ranx, rany, ranz
            LOGICAL :: dusefield, dinterp

            IF (.NOT. ddiffusion) THEN
                RETURN
            END IF

            CALL start_timer(922)

            IF (dturb_diff)

                CALL get_field(x_f, "X")
                CALL get_field(y_f, "Y")
                CALL get_field(z_f, "Z")

                CALL x_f%get_ptr(x, particle%igrid)
                CALL y_f%get_ptr(y, particle%igrid)
                CALL z_f%get_ptr(z, particle%igrid)

                CALL get_field(diffx_f, "P_DIFF_X")
                CALL get_field(diffy_f, "P_DIFF_Y")
                CALL get_field(diffz_f, "P_DIFF_Z")

                CALL diffx_f%get_ptr(diffx, particle%igrid)
                CALL diffy_f%get_ptr(diffy, particle%igrid)
                CALL diffz_f%get_ptr(diffz, particle%igrid)

                IF (dinterp_diffsuion) THEN

                    CALL get_field(dx_f, "DX")
                    CALL get_field(dy_f, "DY")
                    CALL get_field(dz_f, "DZ")

                    CALL dx_f%get_ptr(dx, particle%igrid)
                    CALL dy_f%get_ptr(dy, particle%igrid)
                    CALL dz_f%get_ptr(dz, particle%igrid)

                    CALL get_field(ddx_f, "DDX")
                    CALL get_field(ddy_f, "DDY")
                    CALL get_field(ddz_f, "DDZ")

                    CALL ddx_f%get_ptr(ddx, particle%igrid)
                    CALL ddy_f%get_ptr(ddy, particle%igrid)
                    CALL ddz_f%get_ptr(ddz, particle%igrid)

                    CALL interpolate_lincon(kk, jj, ii, particle, &
                    D_x, D_y, D_z, diffx, diffy, diffz, x, y, z, dx, dy, dz, ddx, ddy, ddz)

                ELSE

                    CALL get_nearest_val(kk, jj, ii, particle, &
                     D_x, D_y, D_z, diffx, diffy, diffz, x, y, z)

                END IF

            ELSE

                D_x = D(1)
                D_y = D(2)
                D_z = D(3)

            END IF

            ranx = 0.0
            rany = 0.0
            ranz = 0.0

            CALL RANDOM_SEED()
            IF (D_x > 0) THEN
                CALL RANDOM_NUMBER(ranx)
                ranx = ranx - 0.5_realk
            END IF
            IF (D_y > 0) THEN
                CALL RANDOM_NUMBER(rany)
                rany = rany - 0.5_realk
            END IF
            IF (D_z > 0) THEN
                CALL RANDOM_NUMBER(ranz)
                ranz = ranz - 0.5_realk
            END IF

            ! diffusion velocity
            pu_diff = SQRT(2 * D_x / dt) * ranx / SQRT(ranx**2 + rany**2 + ranz**2)
            pv_diff = SQRT(2 * D_y / dt) * rany / SQRT(ranx**2 + rany**2 + ranz**2)
            pw_diff = SQRT(2 * D_z / dt) * ranz / SQRT(ranx**2 + rany**2 + ranz**2)

            ! diffusion length
            pdx_diff = SQRT(2 * D_x * dt) * ranx / SQRT(ranx**2 + rany**2 + ranz**2)
            pdy_diff = SQRT(2 * D_y * dt) * rany / SQRT(ranx**2 + rany**2 + ranz**2)
            pdz_diff = SQRT(2 * D_z * dt) * ranz / SQRT(ranx**2 + rany**2 + ranz**2)

            CALL stop_timer(922)

    END SUBROUTINE get_particle_diffusion_automatic

    SUBROUTINE get_particle_diffusion_homogeneous(particle, dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff)

            ! subroutine arguments
            TYPE(baseparticle_t), INTENT(in) :: particle
            REAL(realk), INTENT(in) :: dt
            REAL(realk), INTENT(out) :: pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff

            ! local variables
            REAL(realk) :: D_x, D_y, D_z, ranx, rany, ranz
            LOGICAL :: dusefield, dinterp

            IF (.NOT. ddiffusion) THEN
                RETURN
            END IF

            CALL start_timer(922)

            ranx = 0.0
            rany = 0.0
            ranz = 0.0

            D_x = D(1)
            D_y = D(2)
            D_z = D(3)

            CALL RANDOM_SEED()
            IF (D_x > 0) THEN
                CALL RANDOM_NUMBER(ranx)
                ranx = ranx - 0.5_realk
            END IF
            IF (D_y > 0) THEN
                CALL RANDOM_NUMBER(rany)
                rany = rany - 0.5_realk
            END IF
            IF (D_z > 0) THEN
                CALL RANDOM_NUMBER(ranz)
                ranz = ranz - 0.5_realk
            END IF

            ! diffusion velocity
            pu_diff = SQRT(2 * D_x / dt) * ranx / SQRT(ranx**2 + rany**2 + ranz**2)
            pv_diff = SQRT(2 * D_y / dt) * rany / SQRT(ranx**2 + rany**2 + ranz**2)
            pw_diff = SQRT(2 * D_z / dt) * ranz / SQRT(ranx**2 + rany**2 + ranz**2)

            ! diffusion length
            pdx_diff = SQRT(2 * D_x * dt) * ranx / SQRT(ranx**2 + rany**2 + ranz**2)
            pdy_diff = SQRT(2 * D_y * dt) * rany / SQRT(ranx**2 + rany**2 + ranz**2)
            pdz_diff = SQRT(2 * D_z * dt) * ranz / SQRT(ranx**2 + rany**2 + ranz**2)

            CALL stop_timer(922)

    END SUBROUTINE get_particle_diffusion_homogeneous

    SUBROUTINE get_particle_diffusion_nearest(particle, dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff, &
                kk, jj, ii, diffx, diffy, diffz, x, y, z)

            ! subroutine arguments
            TYPE(baseparticle_t), INTENT(in) :: particle
            REAL(realk), INTENT(in) :: dt
            REAL(realk), INTENT(out) :: pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff
            INTEGER(intk), INTENT(in) :: kk, jj, ii
            REAL(realk), INTENT(in) :: diffx(kk, jj, ii), diffy(kk, jj, ii), diffz(kk, jj, ii)
            REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk)

            ! local variables
            REAL(realk) :: D_x, D_y, D_z, ranx, rany, ranz
            LOGICAL :: dusefield, dinterp

            IF (.NOT. turb_diff) THEN
                CALL get_particle_diffusion_homogeneous(particle, dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff)
                RETURN
            END IF

            CALL start_timer(922)

            ranx = 0.0
            rany = 0.0
            ranz = 0.0

            CALL get_nearest_val(kk, jj, ii, particle, &
             D_x, D_y, D_z, diffx, diffy, diffz, x, y, z)

            CALL RANDOM_SEED()
            IF (D_x > 0) THEN
                CALL RANDOM_NUMBER(ranx)
                ranx = ranx - 0.5_realk
            END IF
            IF (D_y > 0) THEN
                CALL RANDOM_NUMBER(rany)
                rany = rany - 0.5_realk
            END IF
            IF (D_z > 0) THEN
                CALL RANDOM_NUMBER(ranz)
                ranz = ranz - 0.5_realk
            END IF

            ! diffusion velocity
            pu_diff = SQRT(2 * D_x / dt) * ranx / SQRT(ranx**2 + rany**2 + ranz**2)
            pv_diff = SQRT(2 * D_y / dt) * rany / SQRT(ranx**2 + rany**2 + ranz**2)
            pw_diff = SQRT(2 * D_z / dt) * ranz / SQRT(ranx**2 + rany**2 + ranz**2)

            ! diffusion length
            pdx_diff = SQRT(2 * D_x * dt) * ranx / SQRT(ranx**2 + rany**2 + ranz**2)
            pdy_diff = SQRT(2 * D_y * dt) * rany / SQRT(ranx**2 + rany**2 + ranz**2)
            pdz_diff = SQRT(2 * D_z * dt) * ranz / SQRT(ranx**2 + rany**2 + ranz**2)

            CALL stop_timer(922)

    END SUBROUTINE get_particle_diffusion_nearest

    SUBROUTINE get_particle_diffusion_interpolated(particle, dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff, &
                kk, jj, ii, diffx, diffy, diffz, x, y, z, dx, dy, dz, ddx, ddz, ddy)

            ! subroutine arguments
            TYPE(baseparticle_t), INTENT(in) :: particle
            REAL(realk), INTENT(in) :: dt
            REAL(realk), INTENT(out) :: pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff
            INTEGER(intk), INTENT(in) :: kk, jj, ii
            REAL(realk), INTENT(in) :: diffx(kk, jj, ii), diffy(kk, jj, ii), diffz(kk, jj, ii)
            REAL(realk), INTENT(in) :: x(ii), y(jj), z(kk), dx(ii), dy(jj), dz(kk), ddx(ii), ddy(jj), ddz(kk)

            ! local variables
            REAL(realk) :: D_x, D_y, D_z, ranx, rany, ranz

            IF (.NOT. turb_diff) THEN
                CALL get_particle_diffusion_homogeneous(particle, dt, pu_diff, pv_diff, pw_diff, pdx_diff, pdy_diff, pdz_diff)
                RETURN
            END IF

            CALL start_timer(922)

            ranx = 0.0
            rany = 0.0
            ranz = 0.0

            CALL interpolate_lincon(kk, jj, ii, particle, &
             D_x, D_y, D_z, diffx, diffy, diffz, x, y, z, dx, dy, dz, ddx, ddy, ddz)

            CALL RANDOM_SEED()
            IF (D_x > 0) THEN
                CALL RANDOM_NUMBER(ranx)
                ranx = ranx - 0.5_realk
            END IF
            IF (D_y > 0) THEN
                CALL RANDOM_NUMBER(rany)
                rany = rany - 0.5_realk
            END IF
            IF (D_z > 0) THEN
                CALL RANDOM_NUMBER(ranz)
                ranz = ranz - 0.5_realk
            END IF

            ! diffusion velocity
            pu_diff = SQRT(2 * D_x / dt) * ranx / SQRT(ranx**2 + rany**2 + ranz**2)
            pv_diff = SQRT(2 * D_y / dt) * rany / SQRT(ranx**2 + rany**2 + ranz**2)
            pw_diff = SQRT(2 * D_z / dt) * ranz / SQRT(ranx**2 + rany**2 + ranz**2)

            ! diffusion length
            pdx_diff = SQRT(2 * D_x * dt) * ranx / SQRT(ranx**2 + rany**2 + ranz**2)
            pdy_diff = SQRT(2 * D_y * dt) * rany / SQRT(ranx**2 + rany**2 + ranz**2)
            pdz_diff = SQRT(2 * D_z * dt) * ranz / SQRT(ranx**2 + rany**2 + ranz**2)

            CALL stop_timer(922)

    END SUBROUTINE get_particle_diffusion_interpolated

END MODULE