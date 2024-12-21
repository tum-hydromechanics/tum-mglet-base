MODULE particle_io_mod

    USE HDF5
    USE core_mod
    USE stencilio_mod
    USE particle_list_mod

    IMPLICIT NONE(type, external)

    PRIVATE

    ! TO DO: Consider renaming the quite generic types...

    ! int_stencils_t = INTEGER(intk), ALLOCATABLE :: arr(:) with destructor
    ! real_stencils_t = REAL(realk), ALLOCATABLE :: arr(:) with destructor

    TYPE :: particle_io_t

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

        CHARACTER(len=mglet_filename_max) :: file

    CONTAINS

        PROCEDURE, PUBLIC :: init
        PROCEDURE, PUBLIC :: write_particles_io
        PROCEDURE, PUBLIC :: read_particles_io
        PROCEDURE :: write_particles_list
        PROCEDURE :: read_particles_list

    END TYPE particle_io_t

    PUBLIC :: particle_io_t

CONTAINS

    SUBROUTINE init(this, filename)
        ! Subroutine arguments
        CLASS(particle_io_t), INTENT(inout) :: this
        CHARACTER(*), INTENT(in) :: filename

        ! Function body
        this%file = filename

    END SUBROUTINE init



    SUBROUTINE read_particles_io(this)

        ! Subroutine arguments
        CLASS(particle_io_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(hid_t) :: file_id

        ! Function body
        CALL hdf5common_open(this%file, 'r', file_id)
        CALL this%read_particles_list(file_id, my_particle_list)
        CALL hdf5common_close(file_id)

    END SUBROUTINE read_particles_io



    SUBROUTINE write_particles_io(this)

        ! Subroutine arguments
        CLASS(particle_io_t), INTENT(inout) :: this

        ! Local variables
        INTEGER(hid_t) :: file_id

        ! Function body
        CALL hdf5common_open(this%file, 'w', file_id)
        CALL this%write_particles_list(file_id, my_particle_list)
        CALL hdf5common_close(file_id)

    END SUBROUTINE write_particles_io



    SUBROUTINE write_particles_list(this, file_id, plist)

        ! Subroutine arguments
        CLASS(particle_io_t), INTENT(inout) :: this
        INTEGER(hid_t), INTENT(in) :: file_id
        TYPE(particle_list_t), INTENT(inout) :: plist

        ! Local variables
        INTEGER(intk) :: ig, ip, ic, igrid, npart
        INTEGER(intk), ALLOCATABLE :: icount(:)

        ! Function body

        ! Allocating one array for each grid on process
        IF (.NOT. ALLOCATED(this%nparticle)) THEN
           ALLOCATE(this%nparticle(nmygrids))
           this%nparticle = 0
        END IF

        ALLOCATE(this%states_lists(nmygrids))

        ALLOCATE(this%ipart_lists(nmygrids))
        ALLOCATE(this%igrid_lists(nmygrids))
        ALLOCATE(this%islice_lists(nmygrids))

        ALLOCATE(this%gitstep_lists(nmygrids))
        ALLOCATE(this%sitstep_lists(nmygrids))

        ALLOCATE(this%x_lists(nmygrids))
        ALLOCATE(this%y_lists(nmygrids))
        ALLOCATE(this%z_lists(nmygrids))

        ! Counting the particles per grid
        CALL plist%defragment()

        DO ip = 1, plist%ifinal
            IF ( plist%particles(ip)%state > 0 ) THEN
                DO ig = 1, nmygrids
                    igrid = mygrids(ig)
                    IF ( plist%particles(ip)%igrid == igrid ) THEN
                        this%nparticle(ig) = this%nparticle(ig) + 1
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
            npart = this%nparticle(ig)
            ALLOCATE(this%states_lists(ig)%arr(npart))

            ALLOCATE(this%ipart_lists(ig)%arr(npart))
            ALLOCATE(this%igrid_lists(ig)%arr(npart))
            ALLOCATE(this%islice_lists(ig)%arr(npart))

            ALLOCATE(this%gitstep_lists(ig)%arr(npart))
            ALLOCATE(this%sitstep_lists(ig)%arr(npart))

            ALLOCATE(this%x_lists(ig)%arr(npart))
            ALLOCATE(this%y_lists(ig)%arr(npart))
            ALLOCATE(this%z_lists(ig)%arr(npart))
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

                    this%states_lists(ig)%arr(ic) = plist%particles(ip)%state

                    this%ipart_lists(ig)%arr(ic) = plist%particles(ip)%ipart
                    this%igrid_lists(ig)%arr(ic) = plist%particles(ip)%igrid
                    this%islice_lists(ig)%arr(ic) = plist%particles(ip)%islice

                    this%gitstep_lists(ig)%arr(ic) = plist%particles(ip)%gitstep
                    this%sitstep_lists(ig)%arr(ic) = plist%particles(ip)%sitstep

                    this%x_lists(ig)%arr(ic) = plist%particles(ip)%x
                    this%y_lists(ig)%arr(ic) = plist%particles(ip)%y
                    this%z_lists(ig)%arr(ic) = plist%particles(ip)%z

                    ! EXIT
                END IF
            END DO
        END DO

        DO ig = 1, nmygrids
            IF ( this%nparticle(ig) /= icount(ig) ) THEN
                igrid = mygrids(ig)
                WRITE(*,*) "Not all slots filled in igrid = ", igrid
                CALL errr(__FILE__, __LINE__)
            END IF
        END DO

        DEALLOCATE(icount)

        ! Using stencils infrastructure for parallel I/O
        ! (functions manage all grids of process)
        CALL stencilio_write(file_id, 'state', this%states_lists)

        CALL stencilio_write(file_id, 'ipart', this%ipart_lists)
        CALL stencilio_write(file_id, 'igrid', this%igrid_lists)
        CALL stencilio_write(file_id, 'islice', this%islice_lists)

        CALL stencilio_write(file_id, 'gitstep', this%gitstep_lists)
        CALL stencilio_write(file_id, 'sitstep', this%sitstep_lists)

        CALL stencilio_write(file_id, 'x', this%x_lists)
        CALL stencilio_write(file_id, 'y', this%y_lists)
        CALL stencilio_write(file_id, 'z', this%z_lists)

        ! Deallocate all allocated attribute arrays
        DEALLOCATE(this%nparticle)

        DEALLOCATE(this%states_lists)

        DEALLOCATE(this%ipart_lists)
        DEALLOCATE(this%igrid_lists)
        DEALLOCATE(this%islice_lists)

        DEALLOCATE(this%gitstep_lists)
        DEALLOCATE(this%sitstep_lists)

        DEALLOCATE(this%x_lists)
        DEALLOCATE(this%y_lists)
        DEALLOCATE(this%z_lists)

    END SUBROUTINE write_particles_list



    SUBROUTINE read_particles_list(this, file_id, plist)

        ! Subroutine arguments
        CLASS(particle_io_t), INTENT(inout) :: this
        INTEGER(hid_t), INTENT(in) :: file_id
        TYPE(particle_list_t), INTENT(inout) :: plist

        ! Local variables
        INTEGER(intk) :: ig, igrid, npart, n, addlen, cpart, i

        ! Function body
        ALLOCATE(this%states_lists(nmygrids))

        ALLOCATE(this%ipart_lists(nmygrids))
        ALLOCATE(this%igrid_lists(nmygrids))
        ALLOCATE(this%islice_lists(nmygrids))

        ALLOCATE(this%gitstep_lists(nmygrids))
        ALLOCATE(this%sitstep_lists(nmygrids))

        ALLOCATE(this%x_lists(nmygrids))
        ALLOCATE(this%y_lists(nmygrids))
        ALLOCATE(this%z_lists(nmygrids))

        ! Using stencils infrastructure for parallel I/O
        ! (functions manage all grids of process)
        CALL stencilio_read(file_id, 'state', this%states_lists)

        CALL stencilio_read(file_id, 'ipart', this%ipart_lists)
        CALL stencilio_read(file_id, 'igrid', this%igrid_lists)
        CALL stencilio_read(file_id, 'islice', this%islice_lists)

        CALL stencilio_read(file_id, 'gitstep', this%gitstep_lists)
        CALL stencilio_read(file_id, 'sitstep', this%sitstep_lists)

        CALL stencilio_read(file_id, 'x', this%x_lists)
        CALL stencilio_read(file_id, 'y', this%y_lists)
        CALL stencilio_read(file_id, 'z', this%z_lists)

        ! Determine the number of particles
        npart = 0
        DO ig = 1, nmygrids
            IF (ALLOCATED(this%states_lists(ig)%arr)) THEN
                n = SIZE(this%states_lists(ig)%arr)
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

            IF (.NOT. ALLOCATED(this%states_lists(ig)%arr)) THEN
                CYCLE
            END IF

            DO i = 1, SIZE(this%states_lists(ig)%arr)

                ! Checking consistency
                igrid = mygrids(ig)
                IF (this%igrid_lists(ig)%arr(i) /= igrid) THEN
                    WRITE(*,*) "Particle for wrong grid read in"
                    CALL errr(__FILE__, __LINE__)
                END IF

                ! Incrementing the particle counter
                cpart = cpart + 1

                ! Inserting the particle data
                plist%particles(cpart)%state = this%states_lists(ig)%arr(i)

                plist%particles(cpart)%ipart = this%ipart_lists(ig)%arr(i)
                plist%particles(cpart)%igrid = this%igrid_lists(ig)%arr(i)
                plist%particles(cpart)%islice = this%islice_lists(ig)%arr(i)

                plist%particles(cpart)%gitstep = this%gitstep_lists(ig)%arr(i)
                plist%particles(cpart)%sitstep = this%sitstep_lists(ig)%arr(i)

                plist%particles(cpart)%x = this%x_lists(ig)%arr(i)
                plist%particles(cpart)%y = this%y_lists(ig)%arr(i)
                plist%particles(cpart)%z = this%z_lists(ig)%arr(i)

            END DO

        END DO

        IF (cpart /= npart) THEN
            WRITE(*,*) "Counter unequal number of particles"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Deallocate all allocated attribute arrays
        DEALLOCATE(this%states_lists)

        DEALLOCATE(this%ipart_lists)
        DEALLOCATE(this%igrid_lists)
        DEALLOCATE(this%islice_lists)

        DEALLOCATE(this%gitstep_lists)
        DEALLOCATE(this%sitstep_lists)

        DEALLOCATE(this%x_lists)
        DEALLOCATE(this%y_lists)
        DEALLOCATE(this%z_lists)

    END SUBROUTINE read_particles_list

END MODULE particle_io_mod
