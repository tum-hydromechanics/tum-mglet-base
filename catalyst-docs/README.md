# Catalyst and MGLET



## Workflow

Running Catalyst alongside a simulation can be archieved by providing a simple configuration in the parameters `JSON` file required by MGLET. The configuration supports the following arguments
- `catalyst`: Path of the directory that contains Catalyst, including `libcatalyst-paraview.so`
- `repr` [Optional]: Boolean whether catalyst shall create a representative dataset
- `script` [Optional]: Path to catalyst pipeline Python script

The following snippet shows an example on how such a configuration could look like:

```json
{
    // Snip: Simulation parameters (e.g. flow, time)

    "catalyst": {
        "path": "/usr/local/lib/catalyst/",
        "repr": false,
        "script": "catalyst_pipeline.py"
    }
}
```

### Preparing a Simulation with Catalyst

Before running the complete simulation with Catalyst, a representative dataset must be created. This dataset serves as a reference for paraview and the user to set up the pipeline script which is used to specify e.g. which fields, areas of the domain shall be investigaged. Furthermore, the pipeline script serves as the only interface on how to view the simulation during execution.

**Warning** Representative datasets will save all fields of the whole domain to disk. Use with caution for large domains. 

### Running a Simulation with Catalyst

### Catalyst Results
