MODULE particle_loadbalance_mod

    USE MPI_f08
    USE comms_mod
    USE grids_mod

    USE particle_runtimestat_mod, ONLY: psim_n_sent
    USE particle_list_mod
    USE particle_statistics_mod
    USE particle_utils_mod

        PRIVATE

        ! configuration parameters
        ! TODO: make some of these config parameters
        INTEGER(intk) :: loadbalance_step
        INTEGER(intk) :: exess_tolerance_nodes_abs, exess_tolerance_procs_abs, min_particles_to_offload, max_iterations 
        REAL(realk) :: exess_tolerance_nodes_rel, exess_tolerance_procs_rel

        ! Arrays that store information on the gridwise!!! (cpu to cpu) offloading connections of this proc
        !
        !   Dim 1: Information about a specific connection
        !   Dim 2: The different connections
        !   ...

        ! ... in my_helper_conns, the information in the first dimension is sorted as follows:
        !   Field 1: rank in mpi_world_comm of "helpers" (procs to which particles of this proc are offloaded)
        !   Field 2: corresponding grid
        !   Field 3: number of particles
        INTEGER(intk), ALLOCATABLE :: my_helper_conns(:,:)

        ! ... in my_burden_conns, the information in the first dimension is sorted as follows:
        !   Field 1: rank in mpi_world_comm of "burdens" (procs which offloads particles to this proc)
        !   Field 2: corresponding grid
        !   Field 3: number of particles
        INTEGER(intk), ALLOCATABLE :: my_burden_conns(:,:) 

        LOGICAL :: i_am_helper, i_am_burden
        
        INTEGER(intk) :: max_burden_conns, max_helper_conns, n_burden_conns, n_helper_conns

        PUBLIC :: loadbalance_step
        PUBLIC :: init_particle_loadbalance, set_loadbalance_connections

    CONTAINS

        SUBROUTINE init_particle_loadbalance()

            ! TODO: find a more reasonable measure 
            max_helper_conns = ngrid * numprocs +100
            ALLOCATE(my_helper_conns(3, max_helper_conns))

            ! TODO: find a more reasonable measure 
            max_burden_conns = ngrid * numprocs +100
            ALLOCATE(my_burden_conns(3, max_burden_conns))

            max_iterations = 10
            loadbalance_step = 10
            min_particles_to_offload = 10 ! TODO: make this a function of the number of cells per grid
            exess_tolerance_nodes_rel = 0.001
            exess_tolerance_procs_rel = 0.001
            
        END SUBROUTINE init_particle_loadbalance

        SUBROUTINE set_loadbalance_connections()

            ! TODO: optimize MPI communication ! ! !

            INTEGER(intk) :: i, j, counter, pgrid
            INTEGER(intk) :: np_per_proc, np_per_node, np_requested
            INTEGER(kind=int32) :: ierr
            LOGICAL :: finished_among_nodes, finished_on_nodes
            
            ! number of grids per process on shmcomm (index corresponds to the rank on shmcomm)
            INTEGER(intk), ALLOCATABLE :: shm_ngrids(:), displ1(:), displ3(:), recvcount_lb(:)

            ! array with that connects shmid to myid
            INTEGER(intk) :: shm_members(0:shmprocs-1), grank_to_shmid(0:numprocs-1)

            ! array that holds the number of particles and corresponding rank for each grid on the shmcomm
            INTEGER(intk), ALLOCATABLE :: shm_particle_gridinfo(:, :)
            
            ! absolute exess of particles per rank in mpi_comm_world
            INTEGER(intk) :: particle_exess(0:numprocs-1)
            ! absolute exess of particles per rank in mpi_comm_world
            INTEGER(intk) :: particle_exess_on_nodes(0:shmprocs-1)
            ! absolute exess of particles on shmcomm per node (rank in shm_masters_comm)
            INTEGER(intk) :: particle_exess_among_nodes(0:num_shm_masters-1)
            
            INTEGER(intk) :: conns_buffer(3), burden_buffer(3), helper_buffer(2)

            INTEGER(intk) :: offload_order_on_nodes(3, 0:shmprocs-1)
            INTEGER(intk) :: my_burden_potential_shm(2), max_burden_potential_shm(2)
            INTEGER(intk) :: my_helper_potential_shm(2), max_helper_potential_shm(2)
            
            INTEGER(intk) :: burden_info(3)
            INTEGER(intk) :: helper_info(2)

            INTEGER(intk) :: offload_order_among_nodes(3, 0:num_shm_masters-1) 
            INTEGER(int32) :: my_burden_potential_shmmasters(2)
            INTEGER(int32) :: max_burden_potential_shmmasters(2)

            INTEGER(intk), ALLOCATABLE :: grids_np_temp(:)

            INTEGER(intk) :: offload_potential(0:shmprocs-1)
            
            shm_members(shmid) = myid
            CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, shm_members, 1, &
             mglet_mpi_int, shmcomm)

            grank_to_shmid(myid) = shmid
            CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, grank_to_shmid, 1, &
             mglet_mpi_int, MPI_COMM_WORLD)

            ! >>> COMPUTE NUMBER OF PARTICLES PER GRID (LOCAL) <<< 
            ! TODO: replace by particle counting routine 
            grids_np = 0
            DO i = 1, my_particle_list%ifinal
                pgrid = my_particle_list%particles(i)%igrid
                grids_np(particle_grid_ptr(pgrid)) = grids_np(particle_grid_ptr(pgrid)) + 1
            END DO

            !guest_grids_np = 0
            !DO i = 1, guest_particle_list%ifinal
            !    pgrid = guest_particle_list%particles(i)%igrid
            !    guest_grids_np(guest_particle_grid_ptr(pgrid)) = guest_grids_np(guest_particle_grid_ptr(pgrid)) + 1
            !END DO
            
            ALLOCATE(grids_np_temp(SIZE(grids_np)))
            grids_np_temp = grids_np

            ! >>> COMPUTE AND DISTRIBUTE PARTICLE EXCESS <<< 
            local_np = my_particle_list%active_np !+ guest_particle_list%active_np
            
            CALL MPI_Allreduce(local_np, global_np, 1, mglet_mpi_int, MPI_SUM, MPI_COMM_WORLD)
            CALL MPI_Allreduce(local_np, node_np, 1, mglet_mpi_int, MPI_SUM, shmcomm)
            
            np_per_node = NINT(REAL(global_np / num_shm_masters))
            np_per_proc = NINT(REAL(global_np / numprocs))

            ! TODO: multiply np_per_proc and np_per_node with node/rank local factor to account for differences in performance
            ! this factor can be updated each timeintegration step (the sum of factors on all ranks should be 1)

            ! >>> DETERMINE PARTICLE EXESS <<<
            particle_exess(myid) = local_np - np_per_proc
            CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess, 1, &
             mglet_mpi_int, MPI_COMM_WORLD)
            
            particle_exess_on_nodes(shmid) = particle_exess(myid)
            CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess_on_nodes, 1, &
                mglet_mpi_int, shmcomm)
 
            IF (shmid == 0) THEN
                particle_exess_among_nodes(shm_masters_id) = node_np - np_per_node
                CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess_among_nodes, 1, &
                mglet_mpi_int, SHM_MASTERS_COMM)
            END IF

            CALL MPI_Bcast(particle_exess_among_nodes, num_shm_masters, mglet_mpi_int, 0, shmcomm)
            
            ! small sanity check
            IF (shm_masters_id_bcast == -1) CALL errr(__FILE__,__LINE__)

            CALL print_particle_exess(shm_masters_id_bcast, SUM(grids_np_temp), particle_exess(myid), particle_exess_among_nodes(shm_masters_id_bcast))

            ! set tolerance 
            exess_tolerance_nodes_abs = np_per_node * exess_tolerance_nodes_rel
            exess_tolerance_procs_abs = np_per_proc * exess_tolerance_procs_rel

            ! >>> REWIND PREVIOUS DISTRIBUTION <<<
            ! TODO: rewind guest grids if possible

            ! TODO: remove this temporary nullyfication
            my_burden_conns = 0
            my_helper_conns = 0
            n_burden_conns = 0
            n_helper_conns = 0

            ! >>> COLLECT PARTICLE GRIDINFO ON SHM <<<
            ALLOCATE(shm_ngrids(0:shmprocs-1))
            CALL MPI_Allgather(nmy_particle_grids, 1, mglet_mpi_int, shm_ngrids, 1, mglet_mpi_int, shmcomm)
            
            ! sanity check
            IF (.NOT. shm_ngrids(shmid) == nmy_particle_grids) CALL errr(__FILE__, __LINE__)

            ALLOCATE(shm_particle_gridinfo(3, SUM(shm_ngrids)))
            
            ALLOCATE(displ1(0:shmprocs-1))
            ALLOCATE(displ3(0:shmprocs-1))
            ALLOCATE(recvcount_lb(0:shmprocs-1))
            displ1 = 0
            displ3 = 0
            recvcount_lb(0) = shm_ngrids(0) * 3_intk
            DO i = 1, shmprocs - 1
                displ1(i) = displ1(i) + shm_ngrids(i-1)
                displ3(i) = displ3(i-1) + recvcount_lb(i-1)
                recvcount_lb(i) = shm_ngrids(i) * 3_intk
            END DO

            DO i = 1, shm_ngrids(shmid)
                shm_particle_gridinfo(1, i + displ1(shmid)) = myid
                shm_particle_gridinfo(2, i + displ1(shmid)) = my_particle_grids(i)
                shm_particle_gridinfo(3, i + displ1(shmid)) = grids_np(i)
            END DO 

            CALL MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
             shm_particle_gridinfo, recvcount_lb, displ3, mglet_mpi_int, shmcomm)

            offload_order_among_nodes = 0
            ! >>> B: DISTRIBUTE AMONG NODES <<<
            IF (num_shm_masters > 1) THEN
                finished_among_nodes = .FALSE.

                ! B.1: sort procs on shm_masters_comm according to their number of exess particles (ascending) 
                ! on all members of shm_masters_comm
                IF (shmid == 0) THEN
                    offload_order_among_nodes(1, shm_masters_id) = myid
                    offload_order_among_nodes(2, shm_masters_id) = shm_masters_id
                    offload_order_among_nodes(3, shm_masters_id) = particle_exess_among_nodes(shm_masters_id)
                    CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, offload_order_among_nodes, 3, &
                    mglet_mpi_int, SHM_MASTERS_COMM)
                    ! sort offload_order_among_nodes from lowest to highest number of particle_exess
                    CALL sort_conns_unique(offload_order_among_nodes, 3, .FALSE.)
                END IF

                CALL MPI_Bcast(offload_order_among_nodes, 3 * num_shm_masters, mglet_mpi_int, 0, shmcomm)

                ! B.2: iterate over members of shm_masters_comm (start with lowet particle exess (= highest capacity to take particles))
                DO i = 0, num_shm_masters-1

                    IF ((-1_intk * offload_order_among_nodes(3, i)) < min_particles_to_offload) finished_among_nodes = .TRUE.
                    IF (finished_among_nodes) EXIT
                    
                    counter = 1
                    ! B.2.1 While the (negative) particle exess is smaller than the threshold,
                    ! find procs that can send particles to this proc 
                    DO WHILE (offload_order_among_nodes(3, i) < SIGN(exess_tolerance_nodes_abs, -1_intk))

                        ! obsolete?
                        IF ((-1_intk * offload_order_among_nodes(3, i)) < min_particles_to_offload) EXIT 

                        burden_info(1) = -1 ! = burden_proc
                        burden_info(2) = -1 ! = burden_grid
                        burden_info(3) = 0  ! = burden_np
                        
                        IF (shmid == 0) THEN
                            ! obsolete?
                            !np_requested = ABS(offload_order_among_nodes(3, i))
                            !CALL MPI_Bcast(np_requested, 1, mglet_mpi_int, offload_order_among_nodes(2, i), SHM_MASTERS_COMM)

                            ! find the proc on node (shmcomm) with the best offloading potential for each shm_master repectively
                            my_burden_potential_shmmasters(1) = 0
                            my_burden_potential_shmmasters(2) = shm_masters_id
                            IF (.NOT. (shm_masters_id == offload_order_among_nodes(2, i))) THEN
                                DO j = 1, SIZE(shm_particle_gridinfo, 2)
                                    IF (my_burden_potential_shmmasters(1) < MIN(particle_exess(shm_particle_gridinfo(1, j)), particle_exess_among_nodes(shm_masters_id), shm_particle_gridinfo(3, j))) THEN
                                        ! TODO: dont use the maximum number of particles on grid but the highest number of particles per cell on that grid instead ?
                                        my_burden_potential_shmmasters(1) = MIN(particle_exess(shm_particle_gridinfo(1, j)), particle_exess_among_nodes(shm_masters_id), shm_particle_gridinfo(3, j))
                                        burden_info(1) = shm_particle_gridinfo(1, j)
                                        burden_info(2) = shm_particle_gridinfo(2, j)
                                        burden_info(3) = my_burden_potential_shmmasters(1)
                                    END IF
                                END DO
                            END IF

                            ! distribute grid (and corresponding rank) best suited to send particles
                            CALL MPI_Allreduce(my_burden_potential_shmmasters, max_burden_potential_shmmasters, 1, MPI_2INTEGER, MPI_MAXLOC, SHM_MASTERS_COMM)
                            
                            CALL MPI_Bcast(burden_info, 3, mglet_mpi_int, max_burden_potential_shmmasters(2), SHM_MASTERS_COMM)

                            CALL MPI_Bcast(burden_info, 3, mglet_mpi_int, 0, shmcomm)
                        ELSE 
                            CALL MPI_Bcast(burden_info, 3, mglet_mpi_int, 0, shmcomm)
                        END IF

                        CALL MPI_Bcast(max_burden_potential_shmmasters, 2, MPI_INTEGER, 0, shmcomm)

                        IF (burden_info(3) < min_particles_to_offload) THEN
                            finished_among_nodes = .TRUE.
                            EXIT
                        END IF  

                        ! sanity check
                        IF (burden_info(1) < 0 .OR. burden_info(1) > numprocs .OR. burden_info(2) < 1 .OR. burden_info(2) > ngrid) THEN
                            CALL errr(__FILE__, __LINE__)
                        END IF

                        IF (burden_info(1) == offload_order_among_nodes(1, i)) THEN
                            WRITE(*,*) "WARNING: unexpected exit criterion met!"
                            finished_among_nodes = .TRUE.
                            EXIT
                        END IF

                        helper_info(1) = 0 ! = helper_proc
                        helper_info(2) = 0 ! = helper_np
                        ! search for the rank on the current target node best suited to accept particles
                        IF (shm_masters_id_bcast == offload_order_among_nodes(2, i)) THEN
                            IF (shmid == 0) THEN
                                DO j = 0, shmprocs-1
                                    IF (SIGN(helper_info(2), -1_intk) > particle_exess_on_nodes(j)) THEN
                                        helper_info(1) = shm_members(j)
                                        helper_info(2) = ABS(particle_exess_on_nodes(j))
                                    END IF 
                                END DO
                                CALL MPI_Bcast(helper_info, 2, mglet_mpi_int, offload_order_among_nodes(2, i), shm_masters_comm)
                            END IF
                            CALL MPI_Bcast(helper_info, 2, mglet_mpi_int, 0, shmcomm)
                        ELSEIF (shmid == 0) THEN
                            CALL MPI_Bcast(helper_info, 2, mglet_mpi_int, offload_order_among_nodes(2, i), shm_masters_comm)
                            CALL MPI_Bcast(helper_info, 2, mglet_mpi_int, 0, shmcomm)
                        ELSE
                            CALL MPI_Bcast(helper_info, 2, mglet_mpi_int, 0, shmcomm)
                        END IF

                        IF (helper_info(2) < min_particles_to_offload) THEN
                                finished_among_nodes = .TRUE.
                                EXIT
                        END IF  

                        IF (myid == helper_info(1)) THEN
                            n_burden_conns = n_burden_conns + 1

                            my_burden_conns(1, n_burden_conns) = burden_info(1)
                            my_burden_conns(2, n_burden_conns) = burden_info(2)
                            my_burden_conns(3, n_burden_conns) = MIN(helper_info(2), burden_info(3))

                            particle_exess(myid) = particle_exess(myid) + my_burden_conns(3, n_burden_conns)
                            particle_exess_on_nodes(shmid) = particle_exess_on_nodes(shmid) + my_burden_conns(3, n_burden_conns)
                            particle_exess_among_nodes(shm_masters_id_bcast) = particle_exess_among_nodes(shm_masters_id_bcast) + my_burden_conns(3, n_burden_conns)

                            local_np = local_np + my_burden_conns(3, n_burden_conns)

                        END IF 
                        
                        IF (myid == burden_info(1)) THEN
                            n_helper_conns = n_helper_conns + 1
                            
                            my_helper_conns(1, n_helper_conns) = helper_info(1)
                            my_helper_conns(2, n_helper_conns) = burden_info(2)
                            my_helper_conns(3, n_helper_conns) = MIN(helper_info(2), burden_info(3))
                            
                            particle_exess(myid) = particle_exess(myid) - my_helper_conns(3, n_helper_conns)
                            particle_exess_on_nodes(shmid) = particle_exess_on_nodes(shmid) - my_helper_conns(3, n_helper_conns)
                            particle_exess_among_nodes(shm_masters_id_bcast) = particle_exess_among_nodes(shm_masters_id_bcast) - my_helper_conns(3, n_helper_conns)

                            grids_np_temp(particle_grid_ptr(burden_info(2))) = grids_np_temp(particle_grid_ptr(burden_info(2))) - my_helper_conns(3, n_helper_conns)

                            local_np = local_np - my_helper_conns(3, n_helper_conns)

                            WRITE(*, '("New CPU-CPU offloading connection: Node ", I3, "  :  Proc ", I3, "  :  Grid ", I3, "  -->  ", I12," particles  -->  Node", I3, "  :  Proc", I3)') & 
                             shm_masters_id_bcast, myid, burden_info(2), my_helper_conns(3, n_helper_conns), offload_order_among_nodes(2, i), helper_info(1)
                        END IF

                        ! exchange new particle exess and shm_particle_gridinfo
                        CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess, 1, &
                         mglet_mpi_int, MPI_COMM_WORLD)

                        CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess_on_nodes, 1, &
                         mglet_mpi_int, shmcomm)

                        IF (INT(shm_masters_id_bcast, int32) == max_burden_potential_shmmasters(2)) THEN
                            CALL MPI_Bcast(particle_exess_among_nodes, 1 * num_shm_masters, mglet_mpi_int, grank_to_shmid(burden_info(1)), shmcomm)
                        END IF 

                        IF (shm_masters_id_bcast == offload_order_among_nodes(2, i)) THEN
                            CALL MPI_Bcast(particle_exess_among_nodes, 1 * num_shm_masters, mglet_mpi_int, grank_to_shmid(helper_info(1)), shmcomm)
                        END IF 

                        IF (shmid == 0) THEN
                            CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess_among_nodes, 1, &
                                mglet_mpi_int, SHM_MASTERS_COMM)
                        END IF

                        CALL MPI_Bcast(particle_exess_among_nodes, 1 * num_shm_masters, mglet_mpi_int, 0, shmcomm)

                        DO j = 1, shm_ngrids(shmid)
                            shm_particle_gridinfo(1, j + displ1(shmid)) = myid
                            shm_particle_gridinfo(2, j + displ1(shmid)) = my_particle_grids(j)
                            shm_particle_gridinfo(3, j + displ1(shmid)) = grids_np_temp(j)
                        END DO 

                        CALL MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        shm_particle_gridinfo, recvcount_lb, displ3, mglet_mpi_int, shmcomm)

                        DO j = 0, num_shm_masters-1
                            offload_order_among_nodes(3, j) = particle_exess_among_nodes(offload_order_among_nodes(2, j))
                        END DO

                        counter = counter + 1
                        IF (counter == max_iterations) THEN
                            WRITE(*, '("Warning in particle_loadbalance_mod.F90: maximum number of interations in while loop reached!")')
                            EXIT
                        END IF

                    END DO

                END DO
            END IF

            CALL MPI_Barrier(MPI_COMM_WORLD)

            CALL print_particle_exess(shm_masters_id_bcast, local_np, particle_exess(myid), particle_exess_among_nodes(shm_masters_id_bcast))

            ! >>> DISTRIBUTE ON NODES <<<
            CALL MPI_Allreduce(local_np, node_np, 1, mglet_mpi_int, MPI_SUM, shmcomm)
            particle_exess_on_nodes(shmid) = local_np - node_np / shmprocs
            CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess_on_nodes, 1, &
                mglet_mpi_int, shmcomm)

            offload_order_on_nodes = 0
            offload_order_on_nodes(1, shmid) = myid
            offload_order_on_nodes(2, shmid) = shmid
            offload_order_on_nodes(3, shmid) = particle_exess_on_nodes(shmid)

            CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, offload_order_on_nodes, 3, &
             mglet_mpi_int, shmcomm)

            ! sort offload_order_on_nodes from lowest to highest number of particle_exess
            CALL sort_conns_unique(offload_order_on_nodes, 3, .FALSE.)

            IF (shmprocs > 1) THEN
            finished_on_nodes = .FALSE.

                DO i = 0, shmprocs-1

                    IF ((-1_intk * offload_order_on_nodes(3, i)) < min_particles_to_offload) finished_on_nodes = .TRUE.
                    IF (finished_on_nodes) EXIT
                    
                    counter = 1

                    DO WHILE (offload_order_on_nodes(3, i) < SIGN(exess_tolerance_procs_abs, -1_intk))

                        !obsolete?
                        IF ((-1_intk * offload_order_on_nodes(3, i)) < min_particles_to_offload) EXIT 

                        helper_info(1) = offload_order_on_nodes(1, i)! = helper_proc
                        helper_info(2) = ABS(offload_order_on_nodes(3, i))! = helper_np

                        !burden_info(1) = -1 ! = burden_proc
                        !burden_info(2) = -1 ! = burden_grid
                        !burden_info(3) = 0  ! = burden_np

                        my_burden_potential_shm(1) = 0
                        my_burden_potential_shm(2) = shmid
                        IF (.NOT. shmid == i) THEN
                            DO j = 1, nmy_particle_grids
                                IF (my_burden_potential_shm(1) < MIN(particle_exess_on_nodes(shmid), grids_np_temp(j))) THEN
                                    ! TODO: dont use the maximum number of particles on grid but the highest number of particles per cell on that grid instead
                                    my_burden_potential_shm(1) = MIN(particle_exess_on_nodes(shmid), grids_np_temp(j))
                                    !burden_info(1) = myid
                                    burden_info(2) = my_particle_grids(j)
                                    !burden_info(3) = my_burden_potential_shm(1)
                                END IF
                            END DO
                        END IF 

                        CALL MPI_Allreduce(my_burden_potential_shm, max_burden_potential_shm, 1, MPI_2INTEGER, MPI_MAXLOC, shmcomm)

                        IF (max_burden_potential_shm(1) < min_particles_to_offload) THEN
                            finished_on_nodes = .TRUE.
                            EXIT
                        END IF  

                        IF (max_burden_potential_shm(2) == offload_order_on_nodes(2, i)) THEN
                            WRITE(*,*) "WARNING: unexpected exit criterion met!"
                            finished_on_nodes = .TRUE.
                            EXIT
                        END IF

                        IF (shmid == i) THEN
                            n_burden_conns = n_burden_conns + 1

                            ! Blocking receive
                            CALL MPI_Recv(conns_buffer, 3, mglet_mpi_int, &
                            max_burden_potential_shm(2), 123, shmcomm, MPI_STATUS_IGNORE)
                            
                            my_burden_conns(:, n_burden_conns) = conns_buffer

                            particle_exess_on_nodes(shmid) = particle_exess_on_nodes(shmid) + my_burden_conns(3, n_burden_conns)
                            ! TODO: remove update of particle-exess here (this is only for temporary debugging)
                            particle_exess(myid) = particle_exess(myid) + my_burden_conns(3, n_burden_conns)

                            local_np = local_np + my_burden_conns(3, n_burden_conns)
                        END IF 

                        IF (shmid == max_burden_potential_shm(2)) THEN
                            n_helper_conns = n_helper_conns + 1

                            my_helper_conns(1, n_helper_conns) = offload_order_on_nodes(1, i)
                            my_helper_conns(2, n_helper_conns) = burden_info(2)
                            my_helper_conns(3, n_helper_conns) = MIN(max_burden_potential_shm(1), helper_info(2))

                            conns_buffer = my_helper_conns(:, n_helper_conns)
                            conns_buffer(1) = myid
                            
                            ! blocking send
                            CALL MPI_Send(conns_buffer, 3, mglet_mpi_int, &
                            offload_order_on_nodes(2, i), 123, shmcomm, ierr)

                            grids_np_temp(particle_grid_ptr(burden_info(2))) = grids_np_temp(particle_grid_ptr(burden_info(2))) - my_helper_conns(3, n_helper_conns)
                            particle_exess_on_nodes(shmid) = particle_exess_on_nodes(shmid) - my_helper_conns(3, n_helper_conns)
                            ! TODO: remove update of particle-exess here (this is only for temporary debugging)
                            particle_exess(myid) = particle_exess(myid) - my_helper_conns(3, n_helper_conns)

                            local_np = local_np - my_helper_conns(3, n_helper_conns)

                            WRITE(*, '("New CPU-CPU offloading connection: Node ", I3, "  :  Proc ", I3, "  :  Grid ", I3, "  -->  ", I12," particles  -->  Node", I3, "  :  Proc", I3)') & 
                            shm_masters_id_bcast, myid, burden_info(2), my_helper_conns(3, n_helper_conns), shm_masters_id_bcast, offload_order_on_nodes(1, i)
                                
                        END IF
                        ! exchange new particle exess and shm_particle_gridinfo
                        CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess_on_nodes, 1, &
                         mglet_mpi_int, shmcomm)

                        DO j = 0, shmprocs-1
                            offload_order_on_nodes(3, j) = particle_exess_on_nodes(offload_order_on_nodes(2, j))
                        END DO

                        counter = counter + 1
                        IF (counter == max_iterations) THEN
                            WRITE(*, '("Warning in particle_loadbalance_mod.F90: maximum number of interations in while loop reached!")')
                            EXIT
                        END IF
                    
                    END DO

                END DO
                
                ! TODO: remove update of particle-exess here (this is only for temporary debugging)
                CALL MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, particle_exess, 1, &
                 mglet_mpi_int, MPI_COMM_WORLD)
            END IF

            CALL print_particle_exess(shm_masters_id_bcast, local_np, particle_exess(myid), particle_exess_among_nodes(shm_masters_id_bcast))

            DEALLOCATE(grids_np_temp)
            DEALLOCATE(shm_ngrids)
            DEALLOCATE(displ1)
            DEALLOCATE(displ3)
            DEALLOCATE(recvcount_lb)
            DEALLOCATE(shm_particle_gridinfo)

        END SUBROUTINE set_loadbalance_connections

        SUBROUTINE finish_loadbalance()

            IF (ALLOCATED(my_helper_conns)) DEALLOCATE(my_helper_conns)
            IF (ALLOCATED(my_burden_conns)) DEALLOCATE(my_burden_conns)

        END SUBROUTINE finish_loadbalance

        SUBROUTINE print_particle_exess(node, proc_np, proc_exess, node_exess)

            ! subroutine arguments
            INTEGER(intk), INTENT(in) :: node, proc_np, proc_exess, node_exess

            ! local variables
            INTEGER(intk) :: i
            INTEGER(intk) :: node_arr(0:numprocs-1), proc_np_arr(0:numprocs-1), proc_exess_arr(0:numprocs-1), node_exess_arr(0:numprocs-1)
            
            CALL MPI_Gather(node, 1, mglet_mpi_int, &
             node_arr, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)

            CALL MPI_Gather(proc_np, 1, mglet_mpi_int, &
             proc_np_arr, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)

            CALL MPI_Gather(proc_exess, 1, mglet_mpi_int, &
             proc_exess_arr, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)

            CALL MPI_Gather(node_exess, 1, mglet_mpi_int, &
             node_exess_arr, 1, mglet_mpi_int, 0, MPI_COMM_WORLD)

            IF (myid == 0) THEN
                WRITE(*, '("   myid   |   node   |   npart   |   exess   |   node_ex  ")')
                DO i = 0, numprocs-1
                    WRITE(*, '(I10, "|", I10, "|", I11, "|", I11, "|", I11)') i, node_arr(i), proc_np_arr(i), proc_exess_arr(i), node_exess_arr(i) 
                END DO
            END IF

        END SUBROUTINE

END MODULE