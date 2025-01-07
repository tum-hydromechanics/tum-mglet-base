MODULE particle_runtimestat_mod

    USE MPI_f08
    USE comms_mod
    USE err_mod
    USE precision_mod
    USE timer_mod
    USE utils_mod

    USE particle_config_mod

    IMPLICIT NONE

    ! RUNTIME STATISTICS / NOTIFICATIONS

    ! Particle Timeintegration
    REAL(realk) :: psim_max_adv_dx = 0.0
    REAL(realk) :: psim_max_adv_dy = 0.0
    REAL(realk) :: psim_max_adv_dz = 0.0

    REAL(realk) :: psim_max_dif_dx = 0.0
    REAL(realk) :: psim_max_dif_dy = 0.0
    REAL(realk) :: psim_max_dif_dz = 0.0

    REAL(realk) :: psim_max_dx = 0.0
    REAL(realk) :: psim_max_dy = 0.0
    REAL(realk) :: psim_max_dz = 0.0

    REAL(realk) :: psim_max_disp = 0.0

    ! Particle Boundaries/ Particle Motion
    INTEGER(intk) :: psim_n_replaced_tot = 0
    INTEGER(intk) :: psim_n_bcerr = 0
    REAL(realk) :: psim_max_bcerr = 0.0

    ! Particle Exchange
    INTEGER(intk) :: psim_n_sent = 0

CONTAINS

    SUBROUTINE itinfo_particles()

        IF (TRIM(particle_terminal) == "none") THEN
            RETURN
        END IF

        CALL reduce_runtimestat()

        IF (myid == 0) THEN
            WRITE(*, '("--- Particle Runtime Statistics ---")')
            WRITE(*, '("MAX ADV DX/DY/DZ (decoupled):   ", 3F14.9)') psim_max_adv_dx, psim_max_adv_dy, psim_max_adv_dz
            WRITE(*, '("MAX DIF DX/DY/DZ (decoupled):   ", 3F14.9)') psim_max_dif_dx, psim_max_dif_dy, psim_max_dif_dz

            WRITE(*, '("MAX DISPLACEMENT:               ", 1F14.9, "   (", 3F14.9, " )")') &
             psim_max_disp, psim_max_dx, psim_max_dy, psim_max_dz

            WRITE(*, '("MAX OBST ERR:                      ", F17.15)') psim_max_bcerr
            WRITE(*, '("NUM PART IN OBST.:              ", I10)') psim_n_bcerr
            WRITE(*, '("NUM PART REPLACED:              ", I10)') psim_n_replaced_tot

            WRITE(*, '("NUM PART SENT:                  ", I10)') psim_n_sent
            WRITE(*, '("-----------------------------------")')
        END IF

        CALL reset_runtimestat()

    END SUBROUTINE itinfo_particles

    SUBROUTINE reduce_runtimestat()

        ! local variables
        INTEGER(intk) :: iSend = 0
        INTEGER(intk) :: source_proc
        INTEGER(intk) :: integer_temp
        REAL(realk) :: real_temp

        ! Particle Timeintegration
        CALL MPI_Reduce(psim_max_adv_dx, real_temp, 1, mglet_mpi_real, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_max_adv_dx = real_temp

        CALL MPI_Reduce(psim_max_adv_dy, real_temp, 1, mglet_mpi_real, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_max_adv_dy = real_temp

        CALL MPI_Reduce(psim_max_adv_dz, real_temp, 1, mglet_mpi_real, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_max_adv_dz = real_temp

        CALL MPI_Reduce(psim_max_dif_dx, real_temp, 1, mglet_mpi_real, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_max_dif_dx = real_temp

        CALL MPI_Reduce(psim_max_dif_dy, real_temp, 1, mglet_mpi_real, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_max_dif_dy = real_temp

        CALL MPI_Reduce(psim_max_dif_dz, real_temp, 1, mglet_mpi_real, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_max_dif_dz = real_temp

        CALL MPI_Allreduce(psim_max_disp, real_temp, 1, mglet_mpi_real, &
         MPI_MAX, MPI_COMM_WORLD)

        IF (reals_are_equal(psim_max_disp, real_temp, EPSILON(psim_max_disp))) THEN
            iSend = myid
        END IF

        CALL MPI_Allreduce(iSend, source_proc, 1, mglet_mpi_int, &
         MPI_MAX, MPI_COMM_WORLD)

        psim_max_disp = real_temp

        CALL MPI_Bcast(psim_max_dx, 1, mglet_mpi_real, &
         source_proc, MPI_COMM_WORLD)

        CALL MPI_Bcast(psim_max_dy, 1, mglet_mpi_real, &
         source_proc, MPI_COMM_WORLD)

        CALL MPI_Bcast(psim_max_dz, 1, mglet_mpi_real, &
         source_proc, MPI_COMM_WORLD)

        ! Particle Boundaries/ Particle Motion
        CALL MPI_Reduce(psim_n_replaced_tot, integer_temp, 1, mglet_mpi_int, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_n_replaced_tot = integer_temp

        CALL MPI_Reduce(psim_n_bcerr, integer_temp, 1, mglet_mpi_int, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_n_bcerr = integer_temp

        CALL MPI_Reduce(psim_max_bcerr, real_temp, 1, mglet_mpi_real, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_max_bcerr = real_temp

        ! Particle Exchange
        CALL MPI_Reduce(psim_n_sent, integer_temp, 1, mglet_mpi_int, &
         MPI_MAX, 0, MPI_COMM_WORLD)
        psim_n_sent = integer_temp

    END SUBROUTINE reduce_runtimestat

    SUBROUTINE reset_runtimestat()

        ! Particle Timeintegration
        psim_max_adv_dx = 0.0
        psim_max_adv_dy = 0.0
        psim_max_adv_dz = 0.0

        psim_max_dif_dx = 0.0
        psim_max_dif_dy = 0.0
        psim_max_dif_dz = 0.0

        psim_max_dx = 0.0
        psim_max_dy = 0.0
        psim_max_dz = 0.0

        psim_max_disp = 0.0

        ! Particle Boundaries/ Particle Motion
        psim_n_replaced_tot = 0
        psim_n_bcerr = 0
        psim_max_bcerr = 0.0

        ! Particle Exchange
        psim_n_sent = 0

    END SUBROUTINE reset_runtimestat

END MODULE particle_runtimestat_mod