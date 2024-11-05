MODULE offload_helper_mod
    USE precision_mod, ONLY: intk, realk
    USE pointers_mod, ONLY: ip3d, ip1d
    USE grids_mod, ONLY: nmygrids, get_mgdims, get_mgbasb, nboconds, get_bc_ctyp
    USE fields_mod
    USE realfield_mod
    USE scacore_mod
    USE flowcore_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Module constants
    INTEGER(intk), PARAMETER :: N_DIMS = 3
    INTEGER(intk), PARAMETER :: N_BASB = 6
    INTEGER(intk), PARAMETER :: N_FACES = 6
    INTEGER(intk), PARAMETER :: N_BC_RANGE = 2
    INTEGER(intk), PARAMETER :: ISCA_FIELD = 1
    INTEGER(intk), PARAMETER :: ISCA_LEVEL = 1

    ! ┌────────────────────────────────────────────────────────────────────────────┐
    ! | Keeps a pointer to the data that is required on the target device          |
    ! | WHY POINTERS?                                                              |
    ! |     - Keep omp directives central for the sake of code readability         |
    ! |     - Prevent any unwanted intereference with the core flow implementation |
    ! |     - Allows to directly map field data without omp directives in fields   |
    ! └────────────────────────────────────────────────────────────────────────────┘
    ! ----- Pointers to fields -----
    ! Grid parameters
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:) :: ip3d_offload, ip1d_offload
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:, :) :: nboconds_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: rdx_offload, rdy_offload, rdz_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: rddx_offload, rddy_offload, rddz_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: dx_offload, dy_offload, dz_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: ddx_offload, ddy_offload, ddz_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: bt_offload
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: g_offload

    ! Flow/Scalar fields
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: u_offload, v_offload, w_offload, t_offload
    
    ! ----- Newly encoded or global arrays -----
    REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:) :: qtu_offload, qtv_offload, qtw_offload, qtt_offload
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:) :: mgdims_offload, mgbasb_offload
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:) :: encoded_ctyp_offload
    INTEGER(intk), POINTER, CONTIGUOUS, DIMENSION(:, :) :: bc_indexing

    ! Make all data available on the target device
    !$omp declare target(ip3d_offload, ip1d_offload, mgdims_offload, mgbasb_offload)
    !$omp declare target(encoded_ctyp_offload, nboconds_offload, bc_indexing)
    !$omp declare target(rdx_offload, rdy_offload, rdz_offload, rddx_offload, rddy_offload, rddz_offload)
    !$omp declare target(dx_offload, dy_offload, dz_offload, ddx_offload, ddy_offload, ddz_offload)
    !$omp declare target(bt_offload)
    !$omp declare target(g_offload)
    !$omp declare target(u_offload, v_offload, w_offload, t_offload)
    !$omp declare target(qtu_offload, qtv_offload, qtw_offload, qtt_offload)


    ! Public subroutines for host
    PUBLIC :: offload_constants, finish_offload_constants

    ! Public variables for host
    PUBLIC :: ISCA_FIELD, ISCA_LEVEL

    ! Public subroutines for device
    PUBLIC :: ptr_to_grid_x, ptr_to_grid_y, ptr_to_grid_z, ptr_to_grid3, get_mgdims_target, get_mgbasb_target, &
        get_encoded_ctyp_offload

    ! Public variables for device
    PUBLIC :: rdx_offload, rdy_offload, rdz_offload, rddx_offload, rddy_offload, rddz_offload, &
        dx_offload, dy_offload, dz_offload, ddx_offload, ddy_offload, ddz_offload, bt_offload, g_offload, &
        u_offload, v_offload, w_offload, t_offload, qtu_offload, qtv_offload, qtw_offload, qtt_offload, nboconds_offload

CONTAINS
    SUBROUTINE offload_constants()        
        CALL map_grid_data()
        CALL map_constant_grid_fields()
        CALL map_bc_data()
        CALL map_bc_encoding()
        CALL map_flow_sca()
    END SUBROUTINE offload_constants

    SUBROUTINE map_grid_data()
        ! Local variables
        INTEGER(intk) :: igrid, i, mgdims_arr_size, kk, jj, ii

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

        ! Create pointers to pointers_mod fields just to have all omp directives to map data in this file
        ip3d_offload => ip3d
        ip1d_offload => ip1d
        
        !$omp target enter data map(to: ip3d_offload, ip1d_offload, mgdims_offload)
    END SUBROUTINE

    SUBROUTINE map_bc_data()
        ! Local variables
        INTEGER(intk) :: igrid, i, nfro, nbac, nrgt, nlft, nbot, ntop

        ! Fill new structure to map mgbasb
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

        ! Create pointers to nboconds just to have all omp directives to map data in this file
        nboconds_offload => nboconds

        !$omp target enter data map(to: nboconds_offload, mgbasb_offload)
    END SUBROUTINE

    SUBROUTINE map_constant_grid_fields()
        ! Local variables
        TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f, dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f, bt_f, g_f

        ! Create copy for grid constants
        CALL get_field(rdx_f, "RDX")
        CALL get_field(rdy_f, "RDY")
        CALL get_field(rdz_f, "RDZ")
        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")
        CALL get_field(bt_f, "BT")
        CALL get_field(g_f, "G")


        
        rdx_offload => rdx_f%arr
        rdy_offload => rdy_f%arr
        rdz_offload => rdz_f%arr
        rddx_offload => rddx_f%arr
        rddy_offload => rddy_f%arr
        rddz_offload => rddz_f%arr
        dx_offload => dx_f%arr
        dy_offload => dy_f%arr
        dz_offload => dz_f%arr
        ddx_offload => ddx_f%arr
        ddy_offload => ddy_f%arr
        ddz_offload => ddz_f%arr
        bt_offload => bt_f%arr
        g_offload => g_f%arr




        !$omp target enter data map(to: rdx_offload, rdy_offload, rdz_offload, rddx_offload, rddy_offload, rddz_offload, &
        !$omp& dx_offload, dy_offload, dz_offload, ddx_offload, ddy_offload, ddz_offload, bt_offload, g_offload)
    END SUBROUTINE

    SUBROUTINE map_flow_sca()
        ! Local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f, sca_f

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(sca_f, scalar(ISCA_FIELD)%name)

        u_offload => u_f%arr
        v_offload => v_f%arr
        w_offload => w_f%arr
        t_offload => sca_f%arr

        ALLOCATE(qtu_offload, source=u_f%arr)
        ALLOCATE(qtv_offload, source=v_f%arr)
        ALLOCATE(qtw_offload, source=w_f%arr)
        ALLOCATE(qtt_offload, source=u_f%arr)

        !$omp target enter data map(to: u_offload, v_offload, w_offload, t_offload, qtu_offload, qtv_offload, qtw_offload, qtt_offload)
    END SUBROUTINE

    SUBROUTINE map_bc_encoding()
        ! Local variables
        INTEGER(intk) :: igrid, iface, ibocd, n_bo_conds, num_bcs, bc_counter, bctypid
        CHARACTER(len=8) :: ctyp

        ! Allocate fields based on the number of boundary conditions
        CALL count_num_bc(num_bcs)
        ALLOCATE(bc_indexing(N_FACES, nmygrids))
        ALLOCATE(encoded_ctyp_offload(num_bcs))

        ! Encode boundary conditions based on grid, face and type
        ! Each face may have multiple boundary conditions
        bc_indexing = 0
        bc_counter = 1
        DO igrid = 1, nmygrids
            DO iface = 1, N_FACES
                n_bo_conds = nboconds(iface, igrid)
                bc_indexing(iface, igrid) = bc_counter

                DO ibocd = 1, n_bo_conds
                    CALL get_bc_ctyp(ctyp, ibocd, iface, igrid)
                    CALL encode_bc_ctyp(bctypid, ctyp)

                    encoded_ctyp_offload(bc_counter) = bctypid
                    bc_counter = bc_counter + 1
                END DO
            END DO
        END DO

        !$omp target enter data map(to: encoded_ctyp_offload, bc_indexing)
    END SUBROUTINE

    SUBROUTINE get_encoded_ctyp_offload(ctyp_encoded, ibocd, iface, igrid)
        !$omp declare target
        INTEGER(intk), INTENT(out) :: ctyp_encoded
        INTEGER(intk), INTENT(in) :: ibocd
        INTEGER(intk), INTENT(in) :: iface
        INTEGER(intk), INTENT(in) :: igrid

        INTEGER(intk) :: istart

        istart = bc_indexing(iface, igrid)

        ctyp_encoded = encoded_ctyp_offload(istart + ibocd - 1)
    END SUBROUTINE

    SUBROUTINE encode_bc_ctyp(bctypid, ctyp)
        INTEGER(intk), INTENT(OUT) :: bctypid
        CHARACTER(len=8), INTENT(IN) :: ctyp

        ! Encodes a ctyp boundary condition type to an integer
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
        CASE ("PAR")
            bctypid = 8
        CASE DEFAULT
            print *, "Could not map BC:", ctyp
        END SELECT
    END SUBROUTINE encode_bc_ctyp

    SUBROUTINE count_num_bc(num_bc)
        ! Subroutine arguments
        INTEGER(intk), INTENT(OUT) :: num_bc

        ! Local variables
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
        !$omp target exit data map(delete: mgdims_offload, ip3d_offload, ip1d_offload)
        !$omp target exit data map(delete: nboconds_offload, mgbasb_offload)
        !$omp target exit data map(delete: rdx_offload, rdy_offload, rdz_offload, rddx_offload, rddy_offload, rddz_offload, &
        !$omp& dx_offload, dy_offload, dz_offload, ddx_offload, ddy_offload, ddz_offload, bt_offload, g_offload)
        !$omp target exit data map(delete: u_offload, v_offload, w_offload, t_offload, qtu_offload, qtv_offload, qtw_offload, qtt_offload)
        !$omp target exit data map(delete: bc_indexing, encoded_ctyp_offload)

        DEALLOCATE(mgdims_offload)
        DEALLOCATE(mgbasb_offload)
        DEALLOCATE(qtu_offload)
        DEALLOCATE(qtv_offload)
        DEALLOCATE(qtw_offload)
        DEALLOCATE(qtt_offload)
        DEALLOCATE(bc_indexing)
        DEALLOCATE(encoded_ctyp_offload)
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
