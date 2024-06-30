MODULE particle_config_mod

USE precision_mod

IMPLICIT NONE

LOGICAL :: dsim_particles = .TRUE.

LOGICAL :: dread_particles = .FALSE.

LOGICAL :: dwrite_particles = .TRUE.

INTEGER(intk) :: psnapshot_step = 2

REAL(realk), PARAMETER :: D = 0.0_realk ! (3_realk) / (10_realk ** 9_realk) ! Diffusion constant in m²/s (homogeneous and isotropic diffusion for now)

END MODULE particle_config_mod