# PlayGround
The objective of this playground is to use CxxOverlappingAMR and adapt it to a case scenario where we could use it for MGLET and try to have a better understanding of Catalyst using AMR. 

We are trying to visualize the mesh with the .print() that it has not been possible due to writing-race between the processors and to `--output name.vthb` so that we could creaty a catalyst_pipeline.py according to the specific mesh we are trying to build, in order to have an image extractor on paraview-catalyst.

Also we are going to try to implement an example by using the 4x3x1 mesh were the lower part will have a refinement.

We were trying to include multimesh to come up with this solution but so far we have not been able to display the multimesh correctly in paraview. This is way we are going to use AMR-example for this 

## Step 1:
Erase mpi dependencies. Goal 13.01.2025

## Step 2:
Create three domains and overwrite the mesh on each of the domains. 20.01.2025

## Step 3:
- Create three domains but erasing any overlapping. 
- Implement velocity (x,y,z), pressure just like AMR. 30.01.2025

## Step 4:
Instead of using structed grid. Try to implement unstructed and randomly select different elements on the domain to have a defined refinement. 











## from previous project
## Build Instructions

Export environment variables for catalyst and the paraview implementation
```
export CATALYST_IMPLEMENTATION_PATHS=<Path to Catalyst Library>
export CATALYST_IMPLEMENTATION_NAME=paraview
```
Build with
```
mkdir build
cd build
cmake ..
cmake --build .
```

## Generating a Representative Dataset

From the overlapping-amr directory root, run the example binary while outputting a representative dataset with
```
./build/bin/pg_mesh1 --output amr-data.vthb
```
A new `amr-data` directory and `amr-data.vthb` file will be created. Import the `amr-data.vthb` into paraview and set `Default Number of Levels` to the desired depth to see the different levels of the tree.

## Running a Simulation with Catalyst

From the overlapping-amr directory root, run the example binary along with catalyst
```
./build/bin/CxxOverlappingAMRExampleV2 catalyst_pipeline.py
```
This will output the simulation details with `print` statements in the `catalyst_execute` function.
