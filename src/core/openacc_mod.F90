MODULE openacc_mod
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_openacc

    !INTEGER, PARAMETER :: acc_device_type = acc_device_nvidia
    !INTEGER :: num_devices!, device_num, device_type
CONTAINS
    SUBROUTINE init_openacc()
        USE openacc

        !dev_type = acc_get_device_type()
        CALL acc_init(acc_device_nvidia)
        !num_devices = acc_get_num_devices(acc_device_nvidia)
        !print*, 'NUMBER OF DEVICES ', num_devices
    END SUBROUTINE init_openacc

END MODULE openacc_mod