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

    TYPE(int_stencils_t), ALLOCATABLE :: states_lists(:)

    TYPE(int_stencils_t), ALLOCATABLE :: ipart_lists(:)
    TYPE(int_stencils_t), ALLOCATABLE :: igrid_lists(:)
    TYPE(int_stencils_t), ALLOCATABLE :: islice_lists(:)

    TYPE(int_stencils_t), ALLOCATABLE :: gitstep_lists(:)
    TYPE(int_stencils_t), ALLOCATABLE :: sitstep_lists(:)

    TYPE(real_stencils_t), ALLOCATABLE :: x_lists(:)
    TYPE(real_stencils_t), ALLOCATABLE :: y_lists(:)
    TYPE(real_stencils_t), ALLOCATABLE :: z_lists(:)

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
        ALLOCATE(islice_lists(nmygrids))

        ALLOCATE(gitstep_lists(nmygrids))
        ALLOCATE(sitstep_lists(nmygrids))

        ALLOCATE(x_lists(nmygrids))
        ALLOCATE(y_lists(nmygrids))
        ALLOCATE(z_lists(nmygrids))

        ! Counting the particles per grid
        CALL plist%defragment()

        DO ip = 1, plist%ifinal
            IF ( plist%particles(ip)%state > 0 ) THEN
                DO ig = 1, nmygrids
                    igrid = mygrids(ig)
                    IF ( plist%particles(ip)%igrid == igrid ) THEN
                        nparticle(ig) = nparticle(ig) + 1
                        ! EXIT
                    END IF
                END DO
            ELSE
                WRITE(*,*) "Inactive particle in defragmented list"
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO


        ! Allocating space for the particles on each grid
        DO ig = 1, nmygrids
            npart = nparticle(ig)
            ALLOCATE(states_lists(ig)%arr(npart))

            ALLOCATE(ipart_lists(ig)%arr(npart))
            ALLOCATE(igrid_lists(ig)%arr(npart))
            ALLOCATE(islice_lists(ig)%arr(npart))

            ALLOCATE(gitstep_lists(ig)%arr(npart))
            ALLOCATE(sitstep_lists(ig)%arr(npart))

            ALLOCATE(x_lists(ig)%arr(npart))
            ALLOCATE(y_lists(ig)%arr(npart))
            ALLOCATE(z_lists(ig)%arr(npart))
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
                    islice_lists(ig)%arr(ic) = plist%particles(ip)%islice

                    gitstep_lists(ig)%arr(ic) = plist%particles(ip)%gitstep
                    sitstep_lists(ig)%arr(ic) = plist%particles(ip)%sitstep

                    x_lists(ig)%arr(ic) = plist%particles(ip)%x
                    y_lists(ig)%arr(ic) = plist%particles(ip)%y
                    z_lists(ig)%arr(ic) = plist%particles(ip)%z

                    ! EXIT
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
        CALL stencilio_write(file_id, 'islice', islice_lists)

        CALL stencilio_write(file_id, 'gitstep', gitstep_lists)
        CALL stencilio_write(file_id, 'sitstep', sitstep_lists)

        CALL stencilio_write(file_id, 'x', x_lists)
        CALL stencilio_write(file_id, 'y', y_lists)
        CALL stencilio_write(file_id, 'z', z_lists)

        ! Deallocate all allocated attribute arrays
        DEALLOCATE(nparticle)

        DEALLOCATE(states_lists)

        DEALLOCATE(ipart_lists)
        DEALLOCATE(igrid_lists)
        DEALLOCATE(islice_lists)

        DEALLOCATE(gitstep_lists)
        DEALLOCATE(sitstep_lists)

        DEALLOCATE(x_lists)
        DEALLOCATE(y_lists)
        DEALLOCATE(z_lists)

    END SUBROUTINE write_particles_list



    SUBROUTINE read_particles_list(file_id, plist)

        ! Subroutine arguments
        INTEGER(hid_t), INTENT(in) :: file_id
        TYPE(particle_list_t), INTENT(inout) :: plist

        ! Local variables
        INTEGER(intk) :: ig, igrid, npart, n, addlen, cpart, i

        ! Function body
        ALLOCATE(states_lists(nmygrids))

        ALLOCATE(ipart_lists(nmygrids))
        ALLOCATE(igrid_lists(nmygrids))
        ALLOCATE(islice_lists(nmygrids))

        ALLOCATE(gitstep_lists(nmygrids))
        ALLOCATE(sitstep_lists(nmygrids))

        ALLOCATE(x_lists(nmygrids))
        ALLOCATE(y_lists(nmygrids))
        ALLOCATE(z_lists(nmygrids))

        ! Using stencils infrastructure for parallel I/O
        ! (functions manage all grids of process)
        CALL stencilio_read(file_id, 'state', states_lists)

        CALL stencilio_read(file_id, 'ipart', ipart_lists)
        CALL stencilio_read(file_id, 'igrid', igrid_lists)
        CALL stencilio_read(file_id, 'islice', islice_lists)

        CALL stencilio_read(file_id, 'gitstep', gitstep_lists)
        CALL stencilio_read(file_id, 'sitstep', sitstep_lists)

        CALL stencilio_read(file_id, 'x', x_lists)
        CALL stencilio_read(file_id, 'y', y_lists)
        CALL stencilio_read(file_id, 'z', z_lists)

        ! Determine the number of particles
        npart = 0
        DO ig = 1, nmygrids
            IF (ALLOCATED(states_lists(ig)%arr)) THEN
                n = SIZE(states_lists(ig)%arr)
                npart = npart + n
            END IF
        END DO

        ! Extend list of necessary
        IF (npart > plist%max_np) THEN
            addlen = npart - plist%max_np
            CALL reallocate_particle_list(plist, addlen)
        END IF

        cpart = 0
        DO ig = 1, nmygrids

            IF (.NOT. ALLOCATED(states_lists(ig)%arr)) THEN
                CYCLE
            END IF

            DO i = 1, SIZE(states_lists(ig)%arr)

                ! Checking consistency
                igrid = mygrids(ig)
                IF (igrid_lists(ig)%arr(i) /= igrid) THEN
                    WRITE(*,*) "Particle for wrong grid read in"
                    CALL errr(__FILE__, __LINE__)
                END IF

                ! Incrementing the particle counter
                cpart = cpart + 1

                ! Inserting the particle data

                plist%particles(cpart)%state = states_lists(ig)%arr(i)

                plist%particles(cpart)%ipart = ipart_lists(ig)%arr(i)
                plist%particles(cpart)%iproc = myid
                plist%particles(cpart)%igrid = igrid_lists(ig)%arr(i)
                plist%particles(cpart)%islice = islice_lists(ig)%arr(i)

                plist%particles(cpart)%gitstep = gitstep_lists(ig)%arr(i)
                plist%particles(cpart)%sitstep = sitstep_lists(ig)%arr(i)

                plist%particles(cpart)%x = x_lists(ig)%arr(i)
                plist%particles(cpart)%y = y_lists(ig)%arr(i)
                plist%particles(cpart)%z = z_lists(ig)%arr(i)

                CALL set_particle_cell(plist%particles(cpart))

            END DO

        END DO

        plist%ifinal = cpart
        plist%active_np = cpart

        IF (list_limit .AND. plist%max_np > plist_len) THEN
            WRITE(*,*) "WARNING in read_particle_list (h5): Specified List Limit exceeded while reading particles!"
        END IF

        IF (cpart /= npart) THEN
            WRITE(*,*) "Counter unequal number of particles"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Deallocate all allocated attribute arrays
        DEALLOCATE(states_lists)

        DEALLOCATE(ipart_lists)
        DEALLOCATE(igrid_lists)
        DEALLOCATE(islice_lists)

        DEALLOCATE(gitstep_lists)
        DEALLOCATE(sitstep_lists)

        DEALLOCATE(x_lists)
        DEALLOCATE(y_lists)
        DEALLOCATE(z_lists)

    END SUBROUTINE read_particles_list

END MODULE particle_io_mod
