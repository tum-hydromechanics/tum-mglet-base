# Overview
This particle module extends the open-source flow solver mglet-base to simulate the dynamics of tracer particles both in (statistically) steady and unsteady flow field.
Currently available under the branch: julius/particle-transport-mpi

# Usage
## What to adjust in compile-time
Set a desired boundary condition mode in bc_coupling_mode, particle_boundaries_mod.F90.
For the time of writing, there are three options: FLOW, SCAL and PART (this might change in future).
For example, if your baseline simulation does *not* include scalar transport, then use FLOW.
After setting an appropriate boundary condition mode, proceed to compile the code.

## What to specify in JSON parameter file
An example of particle parameters specified in parameters.json can be found in the following:

```json
    "particles": {
        "dread_part_h5": false,
        "dread_part_dict": false,
        "dwrite_part_h5": true,
        "terminal": "normal",
        "snapshot_step": 50,
       	"snapshot_npart": 1000,
        "init_npart": 10000,
        "list_length": 10000,
        "dinterp": true,
        "rk_method": "euler",
        "dturb_diff": false,
        "random_walk_mode": "gaussian2",
        "truncation_limit": 2.0,
        "D":[0.005, 0.005, 0.005],
        "dgridstat": true,
	    "rt_ittot_start": 1000,
	    "rt_tstep_max": 9000,
        "slice_direction": "X",
        "nslice_levels": 2,
        "nslices": [1, 5],
        "slice_levels": [0.5, 1.0]
    }
```

| name | type | description |
|------|------|-------------|
| dread_part_h5 | bool | if it reads particles from a h5 file |
| dread_part_dict | bool | if it reads particles from a dict file |
| dread_obst_dict | bool | if it reads obstacles from a dict file |
| dwrite_part_h5 | bool | if it writes out a h5 file at the end |
| dwrite_npc | bool | experimental |
| terminal | string | terminal output mode. normal, none or verbose for debugging |
| snapshot_step | int | every n-th step a VTK snapshot file is written. default none |
| snapshot_npart | int | number of particles to be written out in a snapshot. Initially particles are randomly selected, then they are written out throughout |
| init_npart | int | for initialisation of random particles | 
| init_npart | int | number of particles to be created and distributed if they're not read from a file |
| list_length | int | initial particle array length per MPI rank. if not sufficient, then it's automatically adjusted during runtime |
| dinterp | bool | if field value is interpolated to the particle locations |
| rk_method | string | specify a RK method to advance the particle position. If flow is unsteady, use euler |
| dturb_diff | bool | turbulent diffusivity on/off. irrelevant for time-resolved simulations |
| random_walk_mode | string | set a PDF for random walk. ?, uniform or gaussian2 |
| truncation_limit | real | a multiple of std to truncate and rescale a Gaussian PDF. should be >1.75 |
| D | real vector | molecular diffusivity in each direction |  
| dgridstat | bool | grid-wise particle statistics |
| rt_ittot_start | int | start of particle statistics sampling |
| rt_tstep_max | int | end of particle statistics sampling |
| slice_direction | string |?|
| nslice_levels | int |?|
| nslices | int list |?|
| slice_levels | real list |?|


Note: the code checks auto-generated particles are not distributed outside of fluid cells.
