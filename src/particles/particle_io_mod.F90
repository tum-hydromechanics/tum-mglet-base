MODULE particle_io_mod

    USE MPI_f08
    USE HDF5

    USE err_mod, ONLY: errr
    USE comms_mod
    USE hdf5common_mod
    USE precision_mod

    USE particle_core_mod
    USE particle_list_mod
    USE particle_exchange_mod ! for particle mpi type

    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    PRIVATE

    ! H5 variables
    ! TODO: make particle write mode conditional
    CHARACTER(len=1) :: ph5_writemode = "w"

    CHARACTER(len=12) :: ph5_filename = "particles.h5"

    INTEGER(hid_t) :: ph5_file_id = 0
    INTEGER(hid_t) :: ph5_group_id = 0

    TYPE, BIND(C) :: poffset_t
        INTEGER(c_long_long) :: offset
        INTEGER(c_long_long) :: length
    END TYPE poffset_t

    INTEGER(hid_t) :: poffset_h5t

    INTEGER(hid_t) :: particle_h5t

    ! MPI cariables
    TYPE(MPI_Datatype) :: poffset_mpit

    INTEGER(intk), ALLOCATABLE :: piogr_counts(:)
    INTEGER(intk), ALLOCATABLE :: piogr_displs(:)

    INTEGER(intk), ALLOCATABLE :: pio_counts(:)
    INTEGER(intk), ALLOCATABLE :: pio_displs(:)

    TYPE(baseparticle_t), ALLOCATABLE :: sendBufParticle_grio(:)
    TYPE(baseparticle_t), ALLOCATABLE :: recvBufParticle_grio(:)

    TYPE(MPI_Request), ALLOCATABLE :: piogr_sendReqs(:), piogr_recvReqs(:)


CONTAINS

    SUBROUTINE init_particle_io()

        ! here, counts and displs refers to the particles of each proc respectively
        IF (ioproc) THEN
            ALLOCATE(piogr_counts(0:iogrprocs - 1))
            ALLOCATE(piogr_displs(0:iogrprocs - 1))

            ALLOCATE(pio_counts(0:ioprocs - 1))
            ALLOCATE(pio_displs(0:ioprocs - 1))

            ALLOCATE(piogr_recvReqs(1:iogrprocs - 1)) ! no request for ioproc needed
            ALLOCATE(piogr_sendReqs(1))
        ELSE
            ALLOCATE(piogr_counts(0))
            ALLOCATE(piogr_displs(0))

            ALLOCATE(pio_counts(0))
            ALLOCATE(pio_displs(0))

            ALLOCATE(piogr_recvReqs(0))
            ALLOCATE(piogr_sendReqs(1))
        END IF

        IF (.NOT. ioproc) RETURN

        ! Create HDF5 type for offset,count
        BLOCK
            ! Local variables
            TYPE(poffset_t), TARGET :: foo
            INTEGER(int32) :: hdferr
            INTEGER(hid_t) :: int64_h5t

            int64_h5t = h5kind_to_type(c_long_long, H5_INTEGER_KIND)

            CALL h5tcreate_f(H5T_COMPOUND_F, C_SIZEOF(foo), poffset_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(poffset_h5t, "OFFSET", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%offset)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(poffset_h5t, "LENGTH", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%length)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)
        END BLOCK

        ! Create MPI type for offset,count
        BLOCK
            ! Local variables
            TYPE(poffset_t), TARGET :: foo
            INTEGER(int32) :: blocklen(2)
            TYPE(MPI_Datatype) :: types(2)
            INTEGER(mpi_address_kind) :: base, disp(2)

            blocklen = 1

            CALL MPI_Get_address(foo%offset, disp(1))
            CALL MPI_Get_address(foo%length, disp(2))

            base = disp(1)
            DO i = 1, SIZE(disp)
                disp(i) = disp(i) - base
            END DO

            types(1) = MPI_INTEGER8
            types(2) = MPI_INTEGER8

            CALL MPI_Type_create_struct(2, blocklen, disp, types, &
                poffset_mpit)
            CALL MPI_Type_commit(poffset_mpit)
        END BLOCK

        ! Create HDF5 type for particles
        BLOCK
            ! Local variables
            TYPE(baseparticle_t), TARGET :: foo
            INTEGER(int32) :: hdferr
            INTEGER(hsize_t) :: dim1(1)
            INTEGER(hid_t) :: int64_h5t
            INTEGER(hid_t) :: int64x3_h5t
            INTEGER(hid_t) :: real64_h5t


            ! component datatypes
            int64_h5t = h5kind_to_type(c_intk, H5_INTEGER_KIND) ! or c_longlong?

            real64_h5t = h5kind_to_type(c_realk, H5_REAL_KIND) ! or c_double?

            dim1 = 3
            CALL h5tarray_create_f(c_intk, 1, dim1, int64x3_h5t, ierr)
            IF (ierr /= 0) CALL errr(__FILE__, __LINE__)


            ! create particle hid_t
            CALL h5tcreate_f(H5T_COMPOUND_F, C_SIZEOF(foo), particl_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)


            ! insert particle hid_t type attributes
            CALL h5tinsert_f(particle_h5t, "STATE", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%state)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "IPART", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%ipart)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "IPROC", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%iproc)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "IGRID", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%igrid)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "ISLICE", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%islice)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "GITSTEP", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%gitstep)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "SITSTEP", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%sitstep)), int64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "IJKCELL", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%ijkcell)), int64x3_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "X", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%x)), real64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "Y", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%y)), real64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

            CALL h5tinsert_f(particle_h5t, "Z", &
                 H5OFFSETOF(C_LOC(foo), C_LOC(foo%z)), real64_h5t, hdferr)
            IF (hdferr /= 0) CALL errr(__FILE__, __LINE__)

        END BLOCK

    END SUBROUTINE init_particle_io

    SUBROUTINE write_particles()

        ! local variables
        INTEGER(intk) :: i, iproc, tag = 123
        INTEGER(hid_t) :: ph5_data_id, ph5_filespace
        INTEGER(hsize_t) :: shape1(1)

        ! gather particles on io group master from all processes of the iogrcomm communicator
        CALL MPI_Gather(my_particle_list%ifinal, 1, mglet_mpi_int, &
         piogr_counts, 1, mglet_mpi_int, 0, IOGRCOMM)

        IF (iogrid == 0) THEN
            piogr_displs(0) = 1

            DO i = 1, iogrprocs - 1
                piogr_displs(i) = piogr_displs(i-1) + piogr_counts(i-1)
            END DO

            ALLOCATE(recvBufParticle_grio(SUM(piogr_counts)))

            DO i = 1, piogr_counts(0)
                recvBufParticle_grio(i) = my_particle_list%particles(i)
            END DO
        END IF

        IF (iogrid == 0) THEN
            DO i = 1, iogrprocs - 1
                CALL MPI_Irecv(recvBufParticle_grio(piogr_displs(i)), piogr_counts(i), particle_mpitype, &
                 i, tag, IOGRCOMM,_piogr_recvreqs(i))
            END DO
        ELSE
            DO i = 1, iogrprocs - 1
                CALL MPI_Isend(my_particle_list%particles(1:ifinal), piogr_counts(iogrid), mglet_mpi_int, &
                0, tag, IOGRCOMM, piogr_sendreqs)
            END DO
        END IF

        CALL MPI_Waitall(SIZE(piogr_recvreqs), piogr_recvreqs, MPI_STATUSES_IGNORE)

        shape1 = !number of all particles

        ! write particles
        IF (ioproc) THEN
            CALL hdf5common_dataset_create("PARTICLES", shape1, particle_h5t, &
             ph5_file_id, ph5_data_id, ph5_filespace)
        END IF

    END SUBROUTINE write_particles

    SUBROUTINE begin_particle_writing()

        CALL init_particle_io()

        IF (myid == 0) THEN
            WRITE(*, '("Opening file ", A, " with mode ", A)') &
                TRIM(ph5_filename), ph5_writemode
        END IF

        CALL hdf5common_open(ph5_filename, ph5_writemode, ph5_file_id)

        !CALL hdf5common_group_open("PARTICLES", ph5_file_id, ph5_group_id, track_index=.TRUE.)

    END SUBROUTINE begin_particle_writing

END MODULE