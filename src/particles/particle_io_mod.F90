MODULE particle_io_mod

    USE HDF5
    USE MPI_f08
    USE comms_mod
    USE core_mod
    USE stencilio_mod

    USE particle_list_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    ! TO DO: Consider renaming the quite generic types...

    ! int_stencils_t = INTEGER(intk), ALLOCATABLE :: arr(:) with destructor
    ! real_stencils_t = REAL(realk), ALLOCATABLE :: arr(:) with destructor

    INTEGER(intk), PARAMETER :: particle_schema_version = 1_intk

    TYPE(int_stencils_t), ALLOCATABLE :: states_lists(:)

    TYPE(int_stencils_t), ALLOCATABLE :: ipart_lists(:)
    TYPE(int_stencils_t), ALLOCATABLE :: igrid_lists(:)

    TYPE(real_stencils_t), ALLOCATABLE :: x_lists(:)
    TYPE(real_stencils_t), ALLOCATABLE :: y_lists(:)
    TYPE(real_stencils_t), ALLOCATABLE :: z_lists(:)

#ifdef _MGLET_OPENMP_
    TYPE(int_stencils_t), ALLOCATABLE :: seed_lists(:)
#endif

    INTEGER(intk), ALLOCATABLE :: nparticle(:)


    PUBLIC :: write_particles_h5, read_particles_h5

CONTAINS

    SUBROUTINE read_particles_h5(filename)

        ! Subroutine arguments
        CHARACTER(*), INTENT(in) :: filename

        ! Local variables
        INTEGER(hid_t) :: file_id

        IF (myid == 0) THEN
            IF (TRIM(particle_terminal) == "normal" .OR. TRIM(particle_terminal) == "verbose") THEN
                WRITE(*,'("Reading Particles from particles.h5 file.")')
                WRITE(*, '()')
            END IF
        END IF

        ! Function body
        CALL hdf5common_open(filename, 'r', file_id)
        CALL read_particles_list(file_id, my_particle_list)
        CALL hdf5common_close(file_id)

        ! global_np is the number of particles amongst all processes, held by the particle list module
        CALL MPI_Allreduce(my_particle_list%active_np, global_np, 1, mglet_mpi_int, MPI_SUM, MPI_COMM_WORLD)

    END SUBROUTINE read_particles_h5



    SUBROUTINE write_particles_h5(filename)

        ! Subroutine arguments
        CHARACTER(*), INTENT(in) :: filename

        ! Local variables
        INTEGER(hid_t) :: file_id

        ! Function body
        CALL hdf5common_open(filename, 'w', file_id)
        CALL hdf5common_attr_write('PARTICLE_SCHEMA_VERSION', particle_schema_version, file_id)
        CALL write_particles_list(file_id, my_particle_list)
        CALL hdf5common_close(file_id)

    END SUBROUTINE write_particles_h5



    SUBROUTINE write_particles_list(file_id, plist)

        ! Subroutine arguments
        INTEGER(hid_t), INTENT(in) :: file_id
        TYPE(particle_list_t), INTENT(inout) :: plist

        ! Local variables
        INTEGER(intk) :: ig, ip, ic, igrid, npart
        INTEGER(intk), ALLOCATABLE :: icount(:)

        ! Function body

        ! Allocating one array for each grid on process
        IF (.NOT. ALLOCATED(nparticle)) THEN
           ALLOCATE(nparticle(nmygrids))
           nparticle = 0
        END IF

        ALLOCATE(states_lists(nmygrids))

        ALLOCATE(ipart_lists(nmygrids))
        ALLOCATE(igrid_lists(nmygrids))

        ALLOCATE(x_lists(nmygrids))
        ALLOCATE(y_lists(nmygrids))
        ALLOCATE(z_lists(nmygrids))

#ifdef _MGLET_OPENMP_
        ALLOCATE(seed_lists(nmygrids))
#endif

        ! Counting the particles per grid
        CALL defragment(plist)

        DO ip = 1, plist%ifinal
            IF ( plist%particles(ip)%state > 0 ) THEN
                DO ig = 1, nmygrids
                    igrid = mygrids(ig)
                    IF ( plist%particles(ip)%igrid == igrid ) THEN
                        nparticle(ig) = nparticle(ig) + 1
                        EXIT
                    END IF
                END DO
            ELSE
                WRITE(*,*) "Inactive particle in defragmented list"
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        IF (SUM(nparticle) /= plist%ifinal) THEN
            WRITE(*,*) "Not all active particles belong to a local grid"
            CALL errr(__FILE__, __LINE__)
        END IF


        ! Allocating space for the particles on each grid
        DO ig = 1, nmygrids
            npart = nparticle(ig)
            ALLOCATE(states_lists(ig)%arr(npart))

            ALLOCATE(ipart_lists(ig)%arr(npart))
            ALLOCATE(igrid_lists(ig)%arr(npart))

            ALLOCATE(x_lists(ig)%arr(npart))
            ALLOCATE(y_lists(ig)%arr(npart))
            ALLOCATE(z_lists(ig)%arr(npart))

#ifdef _MGLET_OPENMP_
            ALLOCATE(seed_lists(ig)%arr(npart))
#endif
        END DO

        ! Inserting the particle data
        ALLOCATE(icount(nmygrids))
        icount = 0

        DO ip = 1, plist%ifinal
            DO ig = 1, nmygrids
                igrid = mygrids(ig)
                IF ( plist%particles(ip)%igrid == igrid ) THEN

                    icount(ig) = icount(ig) + 1
                    ic = icount(ig)

                    states_lists(ig)%arr(ic) = plist%particles(ip)%state

                    ipart_lists(ig)%arr(ic) = plist%particles(ip)%ipart
                    igrid_lists(ig)%arr(ic) = plist%particles(ip)%igrid

                    x_lists(ig)%arr(ic) = plist%particles(ip)%x
                    y_lists(ig)%arr(ic) = plist%particles(ip)%y
                    z_lists(ig)%arr(ic) = plist%particles(ip)%z

#ifdef _MGLET_OPENMP_
                    seed_lists(ig)%arr(ic) = plist%particles(ip)%seed
#endif
                    EXIT
                END IF
            END DO
        END DO

        DO ig = 1, nmygrids
            IF ( nparticle(ig) /= icount(ig) ) THEN
                igrid = mygrids(ig)
                WRITE(*,*) "Not all slots filled in igrid = ", igrid
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        DEALLOCATE(icount)

        ! Using stencils infrastructure for parallel I/O
        ! (functions manage all grids of process)
        CALL stencilio_write(file_id, 'state', states_lists)

        CALL stencilio_write(file_id, 'ipart', ipart_lists)
        CALL stencilio_write(file_id, 'igrid', igrid_lists)

        CALL stencilio_write(file_id, 'x', x_lists)
        CALL stencilio_write(file_id, 'y', y_lists)
        CALL stencilio_write(file_id, 'z', z_lists)

#ifdef _MGLET_OPENMP_
        CALL stencilio_write(file_id, 'seed', seed_lists)
#endif

        ! Deallocate all allocated attribute arrays
        DEALLOCATE(nparticle)

        DEALLOCATE(states_lists)

        DEALLOCATE(ipart_lists)
        DEALLOCATE(igrid_lists)

        DEALLOCATE(x_lists)
        DEALLOCATE(y_lists)
        DEALLOCATE(z_lists)

#ifdef _MGLET_OPENMP_
        DEALLOCATE(seed_lists)
#endif
    END SUBROUTINE write_particles_list



    SUBROUTINE read_particles_list(file_id, plist)

        ! Subroutine arguments
        INTEGER(hid_t), INTENT(in) :: file_id
        TYPE(particle_list_t), INTENT(inout) :: plist

        ! Local variables
        INTEGER(intk) :: ig, igrid, npart, n, addlen, cpart, i
        LOGICAL :: has_state, has_igrid
#ifdef _MGLET_OPENMP_
        LOGICAL :: has_seed
#endif

        ! Function body
        ALLOCATE(ipart_lists(nmygrids))
        ALLOCATE(x_lists(nmygrids))
        ALLOCATE(y_lists(nmygrids))
        ALLOCATE(z_lists(nmygrids))

        ! Using stencils infrastructure for parallel I/O
        ! (functions manage all grids of process)
        CALL hdf5common_dataset_exists('state', file_id, has_state)
        CALL hdf5common_dataset_exists('igrid', file_id, has_igrid)

        CALL stencilio_read(file_id, 'ipart', ipart_lists)
        CALL stencilio_read(file_id, 'x', x_lists)
        CALL stencilio_read(file_id, 'y', y_lists)
        CALL stencilio_read(file_id, 'z', z_lists)

        IF (has_state) THEN
            ALLOCATE(states_lists(nmygrids))
            CALL stencilio_read(file_id, 'state', states_lists)
        END IF
        IF (has_igrid) THEN
            ALLOCATE(igrid_lists(nmygrids))
            CALL stencilio_read(file_id, 'igrid', igrid_lists)
        END IF
#ifdef _MGLET_OPENMP_
        CALL hdf5common_dataset_exists('seed', file_id, has_seed)
        IF (has_seed) THEN
            ALLOCATE(seed_lists(nmygrids))
            CALL stencilio_read(file_id, 'seed', seed_lists)
        ELSEIF (myid == 0) THEN
            WRITE(*,*) "WARNING: Particle restart has no RNG seeds; seeds will be reinitialized from particle IDs."
        END IF
#endif

        ! Determine the number of particles
        npart = 0
        DO ig = 1, nmygrids
            IF (ALLOCATED(ipart_lists(ig)%arr)) THEN
                n = SIZE(ipart_lists(ig)%arr)
                npart = npart + n
            END IF
        END DO

        ! Extend list if necessary
        IF (list_limit .AND. npart > plist_len) THEN
            WRITE(*, '("ERROR on Process ", I0, ": Number of particles to be read from particles.h5 exceeds the given list limit!")') myid
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (npart > plist%max_np) THEN
            addlen = npart - plist%max_np
            CALL reallocate_particle_list(plist, addlen)
        END IF

        cpart = 0
        DO ig = 1, nmygrids

            IF (.NOT. ALLOCATED(ipart_lists(ig)%arr)) THEN
                CYCLE
            END IF

            DO i = 1, SIZE(ipart_lists(ig)%arr)

                ! Checking consistency
                igrid = mygrids(ig)
                IF (has_state) THEN
                    IF (states_lists(ig)%arr(i) < 1) THEN
                        WRITE(*,*) "Inactive particle found in restart file"
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END IF
                IF (has_igrid) THEN
                    IF (igrid_lists(ig)%arr(i) /= igrid) THEN
                        WRITE(*,*) "Particle for wrong grid read in"
                        CALL errr(__FILE__, __LINE__)
                    END IF
                END IF

                ! Incrementing the particle counter
                cpart = cpart + 1

                ! The stencil bucket identifies the grid. set_particle marks
                ! the particle active, sets its owner and reconstructs ijkcell.
                CALL set_particle(plist%particles(cpart), ipart_lists(ig)%arr(i), &
                    x_lists(ig)%arr(i), y_lists(ig)%arr(i), z_lists(ig)%arr(i), &
                    iproc=myid, igrid=igrid)
#ifdef _MGLET_OPENMP_
                IF (has_seed) plist%particles(cpart)%seed = seed_lists(ig)%arr(i)
#endif

            END DO

        END DO

        plist%ifinal = cpart
        plist%active_np = cpart

        IF (cpart /= npart) THEN
            WRITE(*,*) "Counter unequal number of particles"
            CALL errr(__FILE__, __LINE__)
        END IF

        CALL sort_by_grid(plist)
        CALL check_plist(plist, .TRUE.)

        ! Deallocate all allocated attribute arrays
        IF (ALLOCATED(states_lists)) DEALLOCATE(states_lists)

        DEALLOCATE(ipart_lists)
        IF (ALLOCATED(igrid_lists)) DEALLOCATE(igrid_lists)

        DEALLOCATE(x_lists)
        DEALLOCATE(y_lists)
        DEALLOCATE(z_lists)

#ifdef _MGLET_OPENMP_
        IF (ALLOCATED(seed_lists)) DEALLOCATE(seed_lists)
#endif
    END SUBROUTINE read_particles_list

END MODULE particle_io_mod
