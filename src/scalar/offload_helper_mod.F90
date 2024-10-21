MODULE offload_helper_mod
    USE precision_mod, ONLY: intk, realk

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: N_DIMS = 3
    INTEGER(intk), PARAMETER :: N_BASB = 6
    INTEGER(intk), PARAMETER :: N_FACES = 6
    INTEGER(intk), PARAMETER :: N_BC_RANGE = 2

    INTEGER(intk), POINTER, CONTIGUOUS :: ip3d_offload(:), ip1d_offload(:)
    INTEGER(intk), POINTER, CONTIGUOUS :: mgdims_offload(:), mgbasb_offload(:)
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:, :) :: nboconds_offload

    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: gp_bc_ctyp_offload
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:) :: encoded_bc_ctyp_offload
 
    REAL(realk), POINTER, CONTIGUOUS :: rddx_offload(:), rddy_offload(:), rddz_offload(:)
    REAL(realk), POINTER, CONTIGUOUS :: ddx_offload(:), ddy_offload(:), ddz_offload(:)
    REAL(realk), POINTER, CONTIGUOUS :: rdx_offload(:), rdy_offload(:), rdz_offload(:)
    REAL(realk), POINTER, CONTIGUOUS :: bt_offload(:)
    REAL(realk), POINTER, CONTIGUOUS :: dx_offload(:), dy_offload(:), dz_offload(:)

    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: u_offload, v_offload, w_offload, t_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: qtu_offload, qtv_offload, qtw_offload

    ! ┌────────────────────────────────────────────────────────────────────────────┐
    ! | Keeps a copy of the data that is required on the target device             |
    ! | This is done so that omp directives are not scattered across the code      |
    ! | and interferes with the core flow implementation                           |
    ! └────────────────────────────────────────────────────────────────────────────┘
    !$omp declare target(ip3d_offload, ip1d_offload, mgdims_offload, mgbasb_offload, &
    !$omp& rddx_offload, rddy_offload, rddz_offload, rdx_offload, rdy_offload, rdz_offload, &
    !$omp& ddx_offload, ddy_offload, ddz_offload, bt_offload, u_offload, v_offload, w_offload, t_offload, &
    !$omp& dx_offload, dy_offload, dz_offload, qtu_offload, qtv_offload, qtw_offload, gp_bc_ctyp_offload, encoded_bc_ctyp_offload, nboconds_offload)
    
    PUBLIC :: ptr_to_grid_x, ptr_to_grid_y, ptr_to_grid_z, ptr_to_grid3, &
        offload_constants, finish_offload_constants, get_mgdims_target, get_mgbasb_target, &
        ddx_offload, ddy_offload, ddz_offload, rdx_offload, rdy_offload, rdz_offload, rddx_offload, rddy_offload, rddz_offload, &
        bt_offload, u_offload, v_offload, w_offload, t_offload, dx_offload, dy_offload, dz_offload, qtu_offload, qtv_offload, qtw_offload, gp_bc_ctyp_offload, encoded_bc_ctyp_offload, nboconds_offload, get_bc_ctyp_offload
CONTAINS
    SUBROUTINE offload_constants()
        USE pointers_mod, ONLY: ip3d, ip1d
        USE grids_mod, ONLY: nmygrids, get_mgdims, get_mgbasb, nboconds, get_bc_ctyp !front, back, right, left, bottom, top
        USE fields_mod
        USE realfield_mod
        USE scacore_mod

        ! Local variablees
        INTEGER(intk) :: igrid, iface, i, mgdims_arr_size, kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop, bc_counter, n_bo_conds, ibocd, bctypid, num_bcs
        TYPE(field_t), POINTER :: rddx_f, rddy_f, rddz_f, bt_f, ddx_f, ddy_f, ddz_f, rdx_f, rdy_f, rdz_f, dx_f, dy_f, dz_f, sca_f
        TYPE(field_t), POINTER :: u_f, v_f, w_f
        CHARACTER(len=8) :: ctyp

        ! Create grids_mod copy to offload
        mgdims_arr_size = N_DIMS * nmygrids
        ALLOCATE(mgdims_offload(mgdims_arr_size))
        DO igrid = 1, nmygrids
            i = (igrid - 1) * N_DIMS + 1
            CALL get_mgdims(kk, jj, ii, igrid)
            mgdims_offload(i) = ii
            mgdims_offload(i+1) = jj
            mgdims_offload(i+2) = kk
        END DO

        ! Create pointers_mod copy to offload
        ALLOCATE(ip3d_offload, source=ip3d)
        ALLOCATE(ip1d_offload, source=ip1d)

        ! Create copy for grid constants
        CALL get_field(bt_f, "BT")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")
        CALL get_field(rdx_f, "RDX")
        CALL get_field(rdy_f, "RDY")
        CALL get_field(rdz_f, "RDZ")
        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        ALLOCATE(bt_offload, source=bt_f%arr)
        ALLOCATE(ddx_offload, source=ddx_f%arr)
        ALLOCATE(ddy_offload, source=ddy_f%arr)
        ALLOCATE(ddz_offload, source=ddz_f%arr)
        ALLOCATE(rdx_offload, source=rdx_f%arr)
        ALLOCATE(rdy_offload, source=rdy_f%arr)
        ALLOCATE(rdz_offload, source=rdz_f%arr)
        ALLOCATE(rddx_offload, source=rddx_f%arr)
        ALLOCATE(rddy_offload, source=rddy_f%arr)
        ALLOCATE(rddz_offload, source=rddz_f%arr)
        ALLOCATE(dx_offload, source=dx_f%arr)
        ALLOCATE(dy_offload, source=dy_f%arr)
        ALLOCATE(dz_offload, source=dz_f%arr)

        ! Precompute boundary conditions
        ALLOCATE(mgbasb_offload(N_BASB * nmygrids))
        DO igrid = 1, nmygrids
            i = (igrid - 1) * N_BASB + 1
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            mgbasb_offload(i) = nfro
            mgbasb_offload(i+1) = nbac
            mgbasb_offload(i+2) = nrgt
            mgbasb_offload(i+3) = nlft
            mgbasb_offload(i+4) = nbot
            mgbasb_offload(i+5) = ntop
        END DO

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(sca_f, scalar(1)%name)
        ALLOCATE(u_offload, source=u_f%arr)
        ALLOCATE(v_offload, source=v_f%arr)
        ALLOCATE(w_offload, source=w_f%arr)
        ALLOCATE(qtu_offload, source=u_f%arr)
        ALLOCATE(qtv_offload, source=v_f%arr)
        ALLOCATE(qtw_offload, source=w_f%arr)
        ALLOCATE(t_offload, source=sca_f%arr)
        ALLOCATE(nboconds_offload, source=nboconds)

        CALL count_num_bc(num_bcs)
        ALLOCATE(gp_bc_ctyp_offload(N_BC_RANGE, N_FACES, nmygrids))
        ALLOCATE(encoded_bc_ctyp_offload(num_bcs))

        gp_bc_ctyp_offload = 0
        bc_counter = 1
        DO igrid = 1, nmygrids
            DO iface = 1, N_FACES
                n_bo_conds = nboconds(iface, igrid)
                gp_bc_ctyp_offload(1, iface, igrid) = bc_counter

                DO ibocd = 1, n_bo_conds
                    CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                    CALL map_bc_type(bctypid, ctyp)

                    encoded_bc_ctyp_offload(bc_counter) = bctypid
                    bc_counter = bc_counter + 1
                END DO

                gp_bc_ctyp_offload(2, iface, igrid) = bc_counter
            END DO
        END DO

        print *, "-----------------------------------"

        !$omp target enter data map(to: ip3d_offload, ip1d_offload, mgdims_offload, mgbasb_offload, &
        !$omp& rddx_offload, rddy_offload, rddz_offload, ddx_offload, ddy_offload, ddz_offload, &
        !$omp& rdx_offload, rdy_offload, rdz_offload, bt_offload, u_offload, v_offload, w_offload, t_offload, &
        !$omp& dx_offload, dy_offload, dz_offload, qtu_offload, qtv_offload, qtw_offload, gp_bc_ctyp_offload, encoded_bc_ctyp_offload, nboconds_offload)
    END SUBROUTINE offload_constants

    SUBROUTINE get_bc_ctyp_offload(ctyp_encoded, ibocd, iface, igrid)
        INTEGER(intk), INTENT(out) :: ctyp_encoded
        INTEGER(intk), INTENT(in) :: ibocd
        INTEGER(intk), INTENT(in) :: iface
        INTEGER(intk), INTENT(in) :: igrid

        INTEGER(intk) :: istart

        istart = gp_bc_ctyp_offload(1, iface, igrid)

        ctyp_encoded = encoded_bc_ctyp_offload(istart + ibocd - 1)
    END SUBROUTINE

    SUBROUTINE map_bc_type(bctypid, ctyp)
        INTEGER(intk), INTENT(OUT) :: bctypid
        CHARACTER(len=8), INTENT(IN) :: ctyp

        SELECT CASE(ctyp)
        CASE ("FIX")
            bctypid = 1
        CASE ("SIO")
            bctypid = 2
        CASE ("CON")
            bctypid = 3
        CASE ("SLI")
            bctypid = 4
        CASE ("SWA")
            bctypid = 5
        CASE ("NOS")
            bctypid = 6
        CASE ("OP1")
            bctypid = 7
        CASE DEFAULT
            print *, "Could not map BC: ", ctyp
        END SELECT
    END SUBROUTINE map_bc_type

    SUBROUTINE count_num_bc(num_bc)
        USE grids_mod, ONLY: nmygrids, nboconds

        INTEGER(intk), INTENT(OUT) :: num_bc

        INTEGER(intk) :: igrid, iface, n_bo_conds

        num_bc = 0

        DO igrid = 1, nmygrids
            DO iface = 1, N_FACES
                n_bo_conds = nboconds(iface, igrid)
                
                num_bc = num_bc + n_bo_conds
            END DO
        END DO
    END SUBROUTINE count_num_bc

    SUBROUTINE finish_offload_constants()
        !$omp target exit data map(delete: ip3d_offload, ip1d_offload, mgdims_offload, mgbasb_offload, &
        !$omp& rddx_offload, rddy_offload, rddz_offload, ddx_offload, ddy_offload, ddz_offload, &
        !$omp& rdx_offload, rdy_offload, rdz_offload, bt_offload, u_offload, v_offload, w_offload, t_offload, &
        !$omp& dx_offload, dy_offload, dz_offload, qtu_offload, qtv_offload, qtw_offload, gp_bc_ctyp_offload, encoded_bc_ctyp_offload, nboconds_offload)

        DEALLOCATE(ip3d_offload)
        DEALLOCATE(ip1d_offload)
        DEALLOCATE(mgdims_offload)
        DEALLOCATE(ddx_offload)
        DEALLOCATE(ddy_offload)
        DEALLOCATE(ddz_offload)
        DEALLOCATE(rdx_offload)
        DEALLOCATE(rdy_offload)
        DEALLOCATE(rdz_offload)
        DEALLOCATE(rddx_offload)
        DEALLOCATE(rddy_offload)
        DEALLOCATE(rddz_offload)
        DEALLOCATE(dx_offload)
        DEALLOCATE(dy_offload)
        DEALLOCATE(dz_offload)
        DEALLOCATE(bt_offload)
        DEALLOCATE(qtu_offload)
        DEALLOCATE(qtv_offload)
        DEALLOCATE(qtw_offload)
        DEALLOCATE(t_offload)
        DEALLOCATE(nboconds_offload)

        DEALLOCATE(mgbasb_offload)
    END SUBROUTINE finish_offload_constants

    SUBROUTINE get_mgbasb_target(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * N_BASB + 1

        nfro = mgbasb_offload(i)
        nbac = mgbasb_offload(i+1)
        nrgt = mgbasb_offload(i+2)
        nlft = mgbasb_offload(i+3)
        nbot = mgbasb_offload(i+4)
        ntop = mgbasb_offload(i+5)
    END SUBROUTINE get_mgbasb_target

    SUBROUTINE get_mgdims_target(kk, jj, ii, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: kk, jj, ii
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * N_DIMS + 1

        ii = mgdims_offload(i)
        jj = mgdims_offload(i+1)
        kk = mgdims_offload(i+2)
    END SUBROUTINE get_mgdims_target

    SUBROUTINE get_len_ii_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i)
    END SUBROUTINE get_len_ii_target

    SUBROUTINE get_len_jj_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i+1)
    END SUBROUTINE get_len_jj_target

    SUBROUTINE get_len_kk_target(len, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(OUT) :: len
        INTEGER(intk), INTENT(IN) :: igrid

        ! Local variablees
        INTEGER(intk) :: i
        i = (igrid - 1) * 3 + 1

        len = mgdims_offload(i+2)
    END SUBROUTINE get_len_kk_target

    SUBROUTINE get_ip1_target(ip1, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip1

        ip1 = ip1d_offload(igrid)
    END SUBROUTINE get_ip1_target

    SUBROUTINE get_ip3_target(ip3, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: ip3

        ip3 = ip3d_offload(igrid)
    END SUBROUTINE get_ip3_target

    SUBROUTINE ptr_to_grid_x(arr_ptr, n_grid, grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(inout) :: grid_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid

        ! Local variables
        INTEGER(intk) :: ip, len
        
        CALL get_ip1_target(ip, n_grid)
        CALL get_len_ii_target(len, n_grid)
        grid_ptr(1:len) => arr_ptr(ip:ip+len-1)
    END SUBROUTINE ptr_to_grid_x
    
    SUBROUTINE ptr_to_grid_y(arr_ptr, n_grid, grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(inout) :: grid_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid

        ! Local variables
        INTEGER(intk) :: ip, len
        
        CALL get_ip1_target(ip, n_grid)
        CALL get_len_jj_target(len, n_grid)
        grid_ptr(1:len) => arr_ptr(ip:ip+len-1)
    END SUBROUTINE ptr_to_grid_y

    SUBROUTINE ptr_to_grid_z(arr_ptr, n_grid, grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(inout) :: grid_ptr(:)
        INTEGER(intk), INTENT(in) :: n_grid
        ! Result variables
        ! Local variables
        INTEGER(intk) :: ip, len
        
        CALL get_ip1_target(ip, n_grid)
        CALL get_len_kk_target(len, n_grid)
        grid_ptr(1:len) => arr_ptr(ip:ip+len-1)
    END SUBROUTINE ptr_to_grid_z

    SUBROUTINE ptr_to_grid3(arr_ptr, n_grid, grid_ptr)
        !$omp declare target
        ! Function arguments
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: arr_ptr(:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(inout) :: grid_ptr(:, :, :)
        INTEGER(intk), INTENT(in) :: n_grid
        ! Result variables
        ! Local variables
        INTEGER(intk) :: ip, ii, jj, kk
        
        CALL get_mgdims_target(kk, jj, ii, n_grid)
        CALL get_ip3_target(ip, n_grid)
        grid_ptr(1:kk, 1:jj, 1:ii) => arr_ptr(ip:ip+kk*jj*ii-1)
    END SUBROUTINE ptr_to_grid3
END MODULE offload_helper_mod
