MODULE particle_config_mod

USE precision_mod

IMPLICIT NONE

LOGICAL :: dsim_particles = .TRUE.

LOGICAL :: dread_particles = .TRUE.

LOGICAL :: dinterp_particles = .FALSE.

LOGICAL :: dwrite_particles = .TRUE.

INTEGER(intk) :: init_npart = 100
INTEGER(intk) :: max_np_per_list = 1000

INTEGER(intk) :: psnapshot_step = 2

REAL(realk), PARAMETER :: D = 1.0_realk / (10.0_realk ** 2.0_realk) ! Diffusion constant in m²/s (homogeneous and isotropic diffusion for now)

END MODULE particle_config_mod